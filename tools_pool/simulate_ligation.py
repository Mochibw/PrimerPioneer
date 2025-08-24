import json
import os
import uuid
import itertools
from typing import List, Dict, Optional, Union, Tuple
from common_utils.sequence import Fragment, SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json

def simulate_ligation(fragments_json_paths: List[str],
                      allow_circularization: bool = True,
                      sticky_end_tolerance: bool = False,
                      dephosphorylated_ends: Optional[List[str]] = None,
                      output_path: Optional[str] = None) -> Dict[str, Union[List[SequenceRecord], str]]:
    """
    模拟DNA片段的连接反应，仅输出能够形成环状的连接产物。
    要求粘性末端互补或都是平末端才能连接。
    
    Args:
        fragments_json_paths: 要连接的片段JSON文件路径列表
        allow_circularization: 是否允许环化连接
        sticky_end_tolerance: 是否允许不完全匹配的粘性末端连接
        dephosphorylated_ends: 已去磷酸化的末端列表（防止自连）
        output_path: 输出JSON文件路径(可选)
        
    Returns:
        包含环状连接产物和状态信息的字典
    """
    # 加载所有片段
    fragments = []
    for frag_path in fragments_json_paths:
        frag_data = load_sequence_from_json(frag_path)
        fragments.append(frag_data)
    
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
                    dephosphorylated_ends,
                    sticky_end_tolerance
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
                "features": product["features"],
                "overhang_5": product["overhang_5"],
                "overhang_3": product["overhang_3"]
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
                                                 dephosphorylated_ends: List[str], tolerance: bool) -> Tuple[Optional[Dict], bool]:
    """尝试将多个片段连接成环状，要求粘性末端互补或都是平末端"""
    if len(fragments) == 0:
        return None, False
    
    # 检查所有连接点的兼容性
    for i in range(len(fragments)):
        current_frag = fragments[i]
        next_frag = fragments[(i + 1) % len(fragments)]
        
        current_name = frag_names[i]
        next_name = frag_names[(i + 1) % len(fragments)]
        
        # 检查去磷酸化状态
        current_dephos = f"{current_name}_3end" in dephosphorylated_ends
        next_dephos = f"{next_name}_5end" in dephosphorylated_ends
        
        if current_dephos and next_dephos:
            return None, False
        
        # 检查末端兼容性
        current_3end = current_frag.get("overhang_3", {"kind": "blunt", "seq": ""})
        next_5end = next_frag.get("overhang_5", {"kind": "blunt", "seq": ""})
        
        if not _are_compatible_ends(current_3end, next_5end, tolerance):
            return None, False
    
    # 所有连接都兼容，构建环状产物
    return _build_circular_product(fragments, frag_names), True

def _are_compatible_ends(end1: Dict, end2: Dict, tolerance: bool) -> bool:
    """检查两个末端是否兼容：粘性末端互补或都是平末端"""
    # 都是平末端
    if end1["kind"] == "blunt" and end2["kind"] == "blunt":
        return True
    
    # 都是5'突出末端
    if end1["kind"] == "5_overhang" and end2["kind"] == "5_overhang":
        return _are_complementary_sticky_ends(end1["seq"], end2["seq"], tolerance)
    
    # 都是3'突出末端
    if end1["kind"] == "3_overhang" and end2["kind"] == "3_overhang":
        return _are_complementary_sticky_ends(end1["seq"], end2["seq"], tolerance)
    
    # 类型不匹配
    return False

def _are_complementary_sticky_ends(seq1: str, seq2: str, tolerance: bool) -> bool:
    """检查两个粘性末端序列是否互补"""
    if tolerance:
        # 允许不完全匹配
        return seq1 == seq2 or _are_reverse_complement(seq1, seq2)
    else:
        # 要求完全互补
        return _are_reverse_complement(seq1, seq2)

def _are_reverse_complement(seq1: str, seq2: str) -> bool:
    """检查两个序列是否反向互补"""
    if len(seq1) != len(seq2):
        return False
    
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    rev_comp_seq2 = ''.join(complement.get(base, base) for base in reversed(seq2))
    
    return seq1 == rev_comp_seq2

def _build_circular_product(fragments: List[Dict], frag_names: List[str]) -> Dict:
    """构建环状连接产物"""
    try:
        # 构建完整序列
        full_sequence = ""
        current_offset = 0
        feature_offsets = [0]
        
        for i, frag in enumerate(fragments):
            if i > 0:
                # 移除前一个片段的3'端重叠序列
                prev_3end = fragments[i-1].get("overhang_3", {"kind": "blunt", "seq": ""})
                overlap_len = len(prev_3end.get("seq", "")) if prev_3end["kind"] != "blunt" else 0
                full_sequence = full_sequence[:-overlap_len] if overlap_len > 0 else full_sequence
                current_offset -= overlap_len
            
            full_sequence += frag["sequence"]
            if i < len(fragments) - 1:
                feature_offsets.append(current_offset + len(frag["sequence"]))
            current_offset += len(frag["sequence"])
        
        # 处理最后一个片段与第一个片段的连接
        last_3end = fragments[-1].get("overhang_3", {"kind": "blunt", "seq": ""})
        overlap_len = len(last_3end.get("seq", "")) if last_3end["kind"] != "blunt" else 0
        if overlap_len > 0:
            full_sequence = full_sequence[:-overlap_len]
        
        new_id = f"circular_{uuid.uuid4().hex[:8]}"
        
        # 合并和调整特征
        adjusted_features = []
        
        for i, frag in enumerate(fragments):
            offset = feature_offsets[i]
            for feat in frag.get("features", []):
                adjusted_feat = feat.copy()
                adjusted_feat["start"] = feat.get("start", 0) + offset
                adjusted_feat["end"] = feat.get("end", 0) + offset
                adjusted_features.append(adjusted_feat)
        
        # 构建名称
        name_parts = [frag.get("name", frag_names[i]) for i, frag in enumerate(fragments)]
        product_name = "Circular-" + "-".join(name_parts)
        
        return {
            "id": new_id,
            "name": product_name,
            "sequence": full_sequence,
            "length": len(full_sequence),
            "circular": True,  # 关键：设置为环状
            "features": adjusted_features,
            "overhang_5": {"kind": "blunt", "seq": ""},  # 环状DNA没有突出端
            "overhang_3": {"kind": "blunt", "seq": ""}   # 环状DNA没有突出端
        }
    except KeyError:
        return None

# 保持向后兼容的辅助函数
def _connect_blunt_ends(frag1: Dict, frag2: Dict) -> Optional[Dict]:
    """连接两个平末端片段（单片段环化）"""
    product, is_circular = _try_circular_ligation_with_complementary_ends([frag1, frag2], ["frag1", "frag2"], [], False)
    return product if is_circular else None

def _connect_sticky_ends(frag1: Dict, frag2: Dict) -> Optional[Dict]:
    """连接两个粘性末端片段（单片段环化）"""
    product, is_circular = _try_circular_ligation_with_complementary_ends([frag1, frag2], ["frag1", "frag2"], [], False)
    return product if is_circular else None
