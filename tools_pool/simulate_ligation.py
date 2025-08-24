import json
import os
import uuid
import itertools
from typing import List, Dict, Optional, Union, Tuple
from common_utils.sequence import Fragment, SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json

def simulate_ligation(fragments_json_paths: List[str],
                      allow_circularization: bool = True,
                      sticky_end_tolerance: bool = False,  # 参数保留但不启用“宽松匹配”
                      dephosphorylated_ends: Optional[List[str]] = None,
                      output_path: Optional[str] = None) -> Dict[str, Union[List[SequenceRecord], str]]:
    """
    模拟DNA片段的连接反应，仅输出能够形成环状的连接产物。
    要求粘性末端互补或都是平末端才能连接。
    
    Args:
        fragments_json_paths: 要连接的片段JSON文件路径列表
        allow_circularization: 是否允许环化连接（当前仅生成环化产物）
        sticky_end_tolerance: 是否允许不完全匹配的粘性末端连接（本实现中不启用宽松匹配）
        dephosphorylated_ends: 已去磷酸化的末端列表（防止自连），命名形如 "frag0_3end", "frag1_5end"
        output_path: 输出JSON文件路径(可选)
        
    Returns:
        包含环状连接产物和状态信息的字典
    """
    # 加载所有片段
    fragments = []
    for frag_path in fragments_json_paths:
        data = load_sequence_from_json(frag_path)  # 可能是顶层对象
        # 如果是顶层切割结果，就展开成若干真实片段
        if isinstance(data, dict) and "fragments" in data:
            for idx, frag in enumerate(data["fragments"]):
                frag = frag.copy()
                # 补一个易读名字，便于日志与调试
                frag["name"] = frag.get("name") or frag.get("id") or f"{os.path.basename(frag_path)}::frag{idx}"
                fragments.append(frag)
        else:
            # 已经是单片段结构
            if "sequence" not in data:
                raise ValueError(f"{frag_path} 不是单片段结构，也不含 fragments 列表，无法解析")
            data["name"] = data.get("name") or data.get("id") or os.path.basename(frag_path)
            fragments.append(data)
    
    # 检查去磷酸化状态
    if dephosphorylated_ends is None:
        dephosphorylated_ends = []
    
    circular_products = []  # 只保存环状产物
    messages = []
    
    # 尝试所有可能的片段组合（从1个片段到所有片段）
    for num_fragments in range(1, len(fragments) + 1):
        # 生成所有可能的片段组合
        for frag_indices in itertools.combinations(range(len(fragments)), num_fragments):
            # 尝试所有排列顺序
            for frag_order in itertools.permutations(frag_indices):
                # 检查这个排列是否能形成环状
                product, is_circular = _try_circular_ligation_with_complementary_ends(
                    [fragments[i] for i in frag_order],
                    [f"frag{i}" for i in frag_order],
                    dephosphorylated_ends
                )
                
                if product and is_circular:
                    # 记录使用的片段数量和顺序信息
                    product["_ligation_info"] = {
                        "fragment_count": num_fragments,
                        "fragment_order": frag_order,
                        "fragment_names": [fragments[i].get("name", f"frag{i}") for i in frag_order]
                    }
                    circular_products.append(product)
                    messages.append(f"成功形成环状产物: {product['name']} (使用{num_fragments}个片段)")
    
    # 按连接次数从少到多排序
    circular_products.sort(key=lambda x: x["_ligation_info"]["fragment_count"])
    
    # 保存结果
    if output_path and circular_products:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        # 移除内部信息，准备输出
        output_products = []
        for product in circular_products:
            output_product = {
                "id": product["id"],
                "name": product["name"],
                "sequence": product["sequence"],
                "length": product["length"],
                "circular": product["circular"],
                "features": product["features"]
            }
            output_products.append(output_product)
        
        with open(output_path, 'w') as f:
            json.dump(output_products, f, indent=2)
        messages.append(f"环状连接产物已保存至: {output_path}")
    
    return {
        "products": circular_products,
        "message": "; ".join(messages) if messages else "未产生环状连接产物"
    }

def _try_circular_ligation_with_complementary_ends(fragments: List[Dict], frag_names: List[str], 
                                                   dephosphorylated_ends: List[str]) -> Tuple[Optional[Dict], bool]:
    """尝试将多个片段连接成环状，要求粘性末端互补或都是平末端"""
    if len(fragments) == 0:
        return None, False
    
    # 检查所有连接点的兼容性
    for i in range(len(fragments)):
        current_frag = fragments[i]
        next_frag = fragments[(i + 1) % len(fragments)]
        
        current_name = frag_names[i]
        next_name = frag_names[(i + 1) % len(fragments)]
        
        # 去磷酸化：两侧不能同时缺失 5' 磷酸
        current_dephos = f"{current_name}_3end" in dephosphorylated_ends
        next_dephos = f"{next_name}_5end" in dephosphorylated_ends
        if current_dephos and next_dephos:
            return None, False
        
        # 检查末端兼容性（严格互补）
        current_3end = current_frag.get("overhang_3", {"kind": "blunt", "seq": ""})
        next_5end = next_frag.get("overhang_5", {"kind": "blunt", "seq": ""})
        
        if not _are_compatible_ends(current_3end, next_5end):
            return None, False
    
    # 所有连接都兼容，构建环状产物
    return _build_circular_product(fragments, frag_names), True


# ========= 兼容性判定相关：规范化 + 严格互补 =========

def _normalize_overhang(end: Dict) -> Tuple[str, str]:
    """
    把 overhang 统一为“沿接合方向（左->右）观察到的 5'→3' 单链序列”。
    返回 (kind_norm, seq_norm)，kind_norm ∈ {"blunt", "5_overhang", "3_overhang"}。
    """
    kind = end.get("kind", "blunt")
    seq = end.get("seq", "").upper()

    if kind == "blunt" or not seq:
        return "blunt", ""

    if kind == "5_overhang":
        # 已经是朝接合方向的 5'→3'
        return "5_overhang", seq

    if kind == "3_overhang":
        # 反向互补，统一到“接合方向的 5'→3'”
        complement = {'A':'T','T':'A','C':'G','G':'C'}
        rc = ''.join(complement.get(b, b) for b in reversed(seq))
        return "3_overhang", rc

    # 未知类型兜底
    return kind, seq


def _are_compatible_ends(end1: Dict, end2: Dict) -> bool:
    """两个末端兼容：平/平；或同极性黏性端且严格反向互补"""
    # 平末端直接兼容
    if end1.get("kind") == "blunt" and end2.get("kind") == "blunt":
        return True

    kind1, seq1 = _normalize_overhang(end1)
    kind2, seq2 = _normalize_overhang(end2)

    # 必须同极性（都 5′ 或都 3′）
    if kind1 == "blunt" or kind2 == "blunt" or kind1 != kind2:
        return False

    # 严格互补（不启用宽松匹配）
    if len(seq1) != len(seq2):
        return False

    complement = {'A':'T','T':'A','C':'G','G':'C'}
    rc2 = ''.join(complement.get(b, b) for b in reversed(seq2))
    return seq1 == rc2


# ========= 构建序列：按极性决定从哪一侧裁剪 =========

def _build_circular_product(fragments: List[Dict], frag_names: List[str]) -> Dict:
    """构建环状连接产物（严格根据 5′/3′ 极性裁剪重叠）"""
    try:
        # 起始片段
        full_sequence = fragments[0]["sequence"]
        feature_offsets = [0]  # 第0个片段从0开始
        current_offset = len(fragments[0]["sequence"])

        # 依次把后续片段接到 full_sequence 右侧
        for i in range(1, len(fragments)):
            prev = fragments[i-1]
            curr = fragments[i]

            prev_3 = prev.get("overhang_3", {"kind":"blunt","seq":""})
            kind_norm, _ = _normalize_overhang(prev_3)
            overlap_len = len(prev_3.get("seq","")) if prev_3.get("kind") != "blunt" else 0

            if kind_norm == "5_overhang":
                # 从“当前片段开头”裁掉 overlap_len
                full_sequence += curr["sequence"][overlap_len:]
                feature_offsets.append(current_offset - overlap_len)  # 当前片段特征整体左移 overlap_len
                current_offset += len(curr["sequence"]) - overlap_len
            elif kind_norm == "3_overhang":
                # 从“已累积序列的尾部”裁掉 overlap_len
                if overlap_len > 0:
                    full_sequence = full_sequence[:-overlap_len]
                    current_offset -= overlap_len
                full_sequence += curr["sequence"]
                feature_offsets.append(current_offset)
                current_offset += len(curr["sequence"])
            else:
                # blunt
                full_sequence += curr["sequence"]
                feature_offsets.append(current_offset)
                current_offset += len(curr["sequence"])

        # 闭环：处理最后一个片段 3′ 端与第一个片段 5′ 端
        last_3 = fragments[-1].get("overhang_3", {"kind":"blunt","seq":""})
        kind_norm_last, _ = _normalize_overhang(last_3)
        overlap_len_last = len(last_3.get("seq","")) if last_3.get("kind") != "blunt" else 0

        if overlap_len_last > 0:
            if kind_norm_last == "3_overhang":
                # 从 full_sequence 尾部裁剪
                full_sequence = full_sequence[:-overlap_len_last]
                # current_offset 表示线性构建过程的末端位置，闭环后不再使用，可不改
            elif kind_norm_last == "5_overhang":
                # 从“起点（第一个片段的开头）”裁剪
                full_sequence = full_sequence[overlap_len_last:]
                # 所有特征整体左移 overlap_len_last
                feature_offsets = [ofs - overlap_len_last for ofs in feature_offsets]

        new_id = f"circular_{uuid.uuid4().hex[:8]}"

        # 合并与位移特征
        adjusted_features = []
        for i, frag in enumerate(fragments):
            base_ofs = feature_offsets[i]
            for feat in frag.get("features", []):
                adjusted = feat.copy()
                adjusted["start"] = feat.get("start", 0) + base_ofs
                adjusted["end"] = feat.get("end", 0) + base_ofs
                adjusted_features.append(adjusted)

        name_parts = [frag.get("name", frag_names[i]) for i, frag in enumerate(fragments)]
        product_name = "Circular-" + "-".join(name_parts)

        return {
            "id": new_id,
            "name": product_name,
            "sequence": full_sequence,
            "length": len(full_sequence),
            "circular": True,
            "features": adjusted_features
        }
    except KeyError:
        return None


# ======= 保持向后兼容的辅助函数（逻辑依旧使用严格互补） =======

def _connect_blunt_ends(frag1: Dict, frag2: Dict) -> Optional[Dict]:
    """连接两个平末端片段（单片段环化）"""
    product, is_circular = _try_circular_ligation_with_complementary_ends([frag1, frag2], ["frag1", "frag2"], [])
    return product if is_circular else None

def _connect_sticky_ends(frag1: Dict, frag2: Dict) -> Optional[Dict]:
    """连接两个粘性末端片段（单片段环化）"""
    product, is_circular = _try_circular_ligation_with_complementary_ends([frag1, frag2], ["frag1", "frag2"], [])
    return product if is_circular else None
