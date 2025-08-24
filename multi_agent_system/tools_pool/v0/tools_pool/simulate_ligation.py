import json
import os
import uuid
from typing import List, Dict, Optional, Union
from common_utils.sequence import Fragment, SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json
from langchain.tools import tool

@tool()
def simulate_ligation(fragments_json_paths: List[str],
                      allow_circularization: bool = True,
                      sticky_end_tolerance: bool = False,
                      dephosphorylated_ends: Optional[List[str]] = None,
                      output_path: Optional[str] = None) -> Dict[str, Union[List[SequenceRecord], str]]:
    """
    模拟DNA片段的连接反应,支持粘性末端和平末端连接。
    
    Args:
        fragments_json_paths: 要连接的片段JSON文件路径列表
        allow_circularization: 是否允许环化连接
        sticky_end_tolerance: 是否允许不完全匹配的粘性末端连接
        dephosphorylated_ends: 已去磷酸化的末端列表（防止自连）
        output_path: 输出JSON文件路径(可选)
        
    Returns:
        包含连接产物和状态信息的字典
    """
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
    
    # 尝试所有可能的片段连接组合
    for i, frag1 in enumerate(fragments):
        for j, frag2 in enumerate(fragments):
            if i == j and not allow_circularization:
                continue  # 跳过自连如果不允许环化
                
            # 检查末端兼容性
            frag1_3end = frag1.get("overhang_3", {"kind": "blunt", "seq": ""})
            frag2_5end = frag2.get("overhang_5", {"kind": "blunt", "seq": ""})
            
            # 检查去磷酸化状态
            frag1_dephos = f"frag{i}_3end" in dephosphorylated_ends
            frag2_dephos = f"frag{j}_5end" in dephosphorylated_ends
            
            # 如果两端都已去磷酸化，不能连接
            if frag1_dephos and frag2_dephos:
                messages.append(f"片段{i}的3'端和片段{j}的5'端都已去磷酸化，无法连接")
                continue
                
            # 检查末端兼容性
            if (frag1_3end["kind"] == "blunt" and frag2_5end["kind"] == "blunt"):
                # 平末端连接
                connected = _connect_blunt_ends(frag1, frag2)
                if connected:
                    products.append(connected)
                    messages.append(f"成功连接片段{i}和片段{j}（平末端）")
                    
            elif (frag1_3end["kind"] in ["5_overhang", "3_overhang"] and 
                  frag2_5end["kind"] in ["5_overhang", "3_overhang"]):
                # 粘性末端连接
                if _are_compatible_sticky_ends(frag1_3end, frag2_5end, sticky_end_tolerance):
                    connected = _connect_sticky_ends(frag1, frag2)
                    if connected:
                        products.append(connected)
                        messages.append(f"成功连接片段{i}和片段{j}（粘性末端）")
    
    # 保存结果
    if output_path and products:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, 'w') as f:
            json.dump(products, f, indent=2)
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
        return end1["seq"] == end2["seq"] or tolerance
    else:  # 3_overhang
        return end1["seq"] == end2["seq"] or tolerance

def _connect_blunt_ends(frag1: Dict, frag2: Dict) -> Optional[SequenceRecord]:
    """连接两个平末端片段"""
    try:
        new_seq = frag1["sequence"] + frag2["sequence"]
        new_id = f"ligated_{uuid.uuid4().hex[:8]}"
        
        return {
            "id": new_id,
            "name": f"Ligated {frag1.get('name', 'frag1')}-{frag2.get('name', 'frag2')}",
            "sequence": new_seq,
            "length": len(new_seq),
            "circular": False,
            "features": frag1.get("features", []) + frag2.get("features", []),
            "overhang_5": frag1.get("overhang_5", {"kind": "blunt", "seq": ""}),
            "overhang_3": frag2.get("overhang_3", {"kind": "blunt", "seq": ""})
        }
    except KeyError:
        return None

def _connect_sticky_ends(frag1: Dict, frag2: Dict) -> Optional[SequenceRecord]:
    """连接两个粘性末端片段"""
    try:
        # 移除重叠的粘性末端序列
        overlap_len = len(frag1.get("overhang_3", {}).get("seq", ""))
        frag1_seq = frag1["sequence"][:-overlap_len] if overlap_len > 0 else frag1["sequence"]
        
        new_seq = frag1_seq + frag2["sequence"]
        new_id = f"ligated_{uuid.uuid4().hex[:8]}"
        
        return {
            "id": new_id,
            "name": f"Ligated {frag1.get('name', 'frag1')}-{frag2.get('name', 'frag2')}",
            "sequence": new_seq,
            "length": len(new_seq),
            "circular": False,
            "features": frag1.get("features", []) + frag2.get("features", []),
            "overhang_5": frag1.get("overhang_5", {"kind": "blunt", "seq": ""}),
            "overhang_3": frag2.get("overhang_3", {"kind": "blunt", "seq": ""})
        }
    except KeyError:
        return None
