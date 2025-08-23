import json
import os
import uuid
from typing import List, Dict, Optional, Union
from common_utils.sequence import Fragment, SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json

def simulate_ligation(fragments_json_paths: List[str],
                      allow_circularization: bool = True,
                      sticky_end_tolerance: bool = False,
                      dephosphorylated_ends: Optional[List[str]] = None,
                      output_path: Optional[str] = None) -> Dict[str, Union[List[SequenceRecord], str]]:
    """
    模拟DNA片段的连接反应，专门用于将片段克隆到环状载体上。
    
    Args:
        fragments_json_paths: 要连接的片段JSON文件路径列表（第一个应为载体）
        allow_circularization: 必须为True，仅支持环状连接
        sticky_end_tolerance: 是否允许不完全匹配的粘性末端连接
        dephosphorylated_ends: 已去磷酸化的末端列表（防止自连）
        output_path: 输出JSON文件路径(可选)
        
    Returns:
        包含连接产物和状态信息的字典
    """
    # 强制仅允许环状连接
    if not allow_circularization:
        return {
            "products": [],
            "message": "仅支持环状连接模式（类似克隆到质粒）"
        }
    
    # 加载所有片段
    fragments = []
    for frag_path in fragments_json_paths:
        frag_data = load_sequence_from_json(frag_path)
        fragments.append(frag_data)
    
    # 检查去磷酸化状态
    if dephosphorylated_ends is None:
        dephosphorylated_ends = []
    
    products = []
    messages = []
    
    # 确保至少有一个载体和一个插入片段
    if len(fragments) < 2:
        return {
            "products": [],
            "message": "需要至少一个载体和一个插入片段进行克隆连接"
        }
    
    # 假设第一个片段是载体，其他是插入片段
    vector = fragments[0]
    inserts = fragments[1:]
    
    # 检查载体是否为环状
    if not vector.get("circular", False):
        messages.append("载体必须是环状DNA")
    
    # 尝试将每个插入片段连接到载体
    for i, insert in enumerate(inserts):
        # 检查载体和插入片段的末端兼容性
        vector_5end = vector.get("overhang_5", {"kind": "blunt", "seq": ""})
        vector_3end = vector.get("overhang_3", {"kind": "blunt", "seq": ""})
        insert_5end = insert.get("overhang_5", {"kind": "blunt", "seq": ""})
        insert_3end = insert.get("overhang_3", {"kind": "blunt", "seq": ""})
        
        # 检查去磷酸化状态
        vector_dephos = f"vector_ends" in dephosphorylated_ends
        insert_dephos = f"insert_{i}_ends" in dephosphorylated_ends
        
        # 如果两端都已去磷酸化，不能连接
        if vector_dephos and insert_dephos:
            messages.append(f"载体和插入片段{i}都已去磷酸化，无法连接")
            continue
            
        # 检查粘性末端兼容性（模拟限制性内切酶切割位点匹配）
        if (_are_compatible_sticky_ends(vector_5end, insert_5end, sticky_end_tolerance) and
            _are_compatible_sticky_ends(vector_3end, insert_3end, sticky_end_tolerance)):
            
            # 成功连接，创建环状产物
            connected = _create_circular_ligation(vector, insert, i)
            if connected:
                products.append(connected)
                messages.append(f"成功将插入片段{i}连接到载体（形成环状质粒）")
        
        # 检查平末端连接
        elif (vector_5end["kind"] == "blunt" and vector_3end["kind"] == "blunt" and
              insert_5end["kind"] == "blunt" and insert_3end["kind"] == "blunt"):
            
            connected = _create_circular_ligation(vector, insert, i)
            if connected:
                products.append(connected)
                messages.append(f"成功将插入片段{i}连接到载体（平末端连接）")
    
    # 保存结果
    if output_path and products:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        # 确保输出格式与generate_map.py兼容
        output_data = []
        for product in products:
            # 转换为generate_map.py期望的格式
            output_data.append({
                "id": product["id"],
                "name": product["name"],
                "sequence": product["sequence"],
                "length": product["length"],
                "circular": product["circular"],
                "features": product["features"],
                # 清除末端突出，因为环状质粒没有突出端
                "overhang_5": {"kind": "blunt", "seq": ""},
                "overhang_3": {"kind": "blunt", "seq": ""}
            })
        
        with open(output_path, 'w') as f:
            json.dump(output_data, f, indent=2)
        messages.append(f"连接产物已保存至: {output_path}")
    
    return {
        "products": products,
        "message": "; ".join(messages) if messages else "未产生连接产物"
    }

def _are_compatible_sticky_ends(end1: Dict, end2: Dict, tolerance: bool) -> bool:
    """检查两个粘性末端是否兼容"""
    if end1["kind"] != end2["kind"]:
        return False
        
    if end1["kind"] == "5_overhang":
        # 5'突出端需要互补配对
        return _are_complementary(end1["seq"], end2["seq"]) or tolerance
    elif end1["kind"] == "3_overhang":
        # 3'突出端需要相同序列
        return end1["seq"] == end2["seq"] or tolerance
    return False

def _are_complementary(seq1: str, seq2: str) -> bool:
    """检查两个序列是否互补"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    comp_seq2 = ''.join(complement.get(base, base) for base in seq2)
    return seq1 == comp_seq2

def _create_circular_ligation(vector: Dict, insert: Dict, insert_idx: int) -> Optional[SequenceRecord]:
    """创建环状连接产物"""
    try:
        # 创建新的环状序列（载体 + 插入片段）
        new_seq = vector["sequence"] + insert["sequence"]
        new_id = f"circular_ligated_{uuid.uuid4().hex[:8]}"
        
        # 合并特征，调整插入片段特征的位置
        adjusted_features = vector.get("features", [])[:]
        insert_features = insert.get("features", [])
        
        for feat in insert_features:
            # 调整插入片段特征的起始位置（加上载体长度）
            adjusted_feat = feat.copy()
            adjusted_feat["start"] = feat.get("start", 0) + len(vector["sequence"])
            adjusted_feat["end"] = feat.get("end", 0) + len(vector["sequence"])
            adjusted_features.append(adjusted_feat)
        
        return {
            "id": new_id,
            "name": f"Circular {vector.get('name', 'vector')} with {insert.get('name', f'insert_{insert_idx}')}",
            "sequence": new_seq,
            "length": len(new_seq),
            "circular": True,  # 关键：设置为环状
            "features": adjusted_features,
            "overhang_5": {"kind": "blunt", "seq": ""},
            "overhang_3": {"kind": "blunt", "seq": ""}
        }
    except KeyError:
        return None

# 保持向后兼容的辅助函数
def _connect_blunt_ends(frag1: Dict, frag2: Dict) -> Optional[SequenceRecord]:
    """连接两个平末端片段（用于环状连接）"""
    return _create_circular_ligation(frag1, frag2, 0)

def _connect_sticky_ends(frag1: Dict, frag2: Dict) -> Optional[SequenceRecord]:
    """连接两个粘性末端片段（用于环状连接）"""
    return _create_circular_ligation(frag1, frag2, 0)
