import json
import os
import uuid
from typing import Dict, List, Optional, Union
from common_utils.sequence import SequenceRecord, Fragment, write_record_to_json, load_sequence_from_json

def simulate_ligation(fragments_paths: List[str],
                     allow_circularization: bool = True,
                     sticky_end_tolerance: bool = False,
                     dephosphorylated_ends: Optional[List[str]] = None,
                     output_path: Optional[str] = None) -> Dict[str, Union[List[SequenceRecord], str]]:
    """
    模拟DNA连接反应，支持平末端和粘性末端连接。
    
    Args:
        fragments_paths: 要连接的DNA片段JSON文件路径列表
        allow_circularization: 是否允许环化连接
        sticky_end_tolerance: 是否允许不完全匹配的粘性末端连接
        dephosphorylated_ends: 已去磷酸化的末端ID列表（防止自连）
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含连接产物和状态信息的字典
    """
    # 加载所有片段
    fragments = []
    for path in fragments_paths:
        fragment = load_sequence_from_json(path)
        fragments.append(fragment)
    
    # 检查去磷酸化状态
    dephosphorylated_ids = dephosphorylated_ends or []
    
    # 模拟连接反应
    products = []
    connected_fragments = set()
    
    # 尝试所有可能的连接组合
    for i, frag1 in enumerate(fragments):
        if frag1["id"] in connected_fragments:
            continue
            
        for j, frag2 in enumerate(fragments):
            if i == j or frag2["id"] in connected_fragments:
                continue
                
            # 检查是否可以连接
            if _can_ligate(frag1, frag2, sticky_end_tolerance, dephosphorylated_ids):
                # 执行连接
                product = _perform_ligation(frag1, frag2)
                products.append(product)
                connected_fragments.add(frag1["id"])
                connected_fragments.add(frag2["id"])
                break
    
    # 处理未连接的片段
    for i, frag in enumerate(fragments):
        if frag["id"] not in connected_fragments:
            products.append(frag)
    
    message = f"连接完成，生成 {len(products)} 个产物"
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        if len(products) == 1:
            write_record_to_json(products[0], output_path)
        else:
            # 保存为产物列表
            with open(output_path, 'w') as f:
                json.dump({"products": products}, f, indent=2)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "products": products,
        "message": message
    }

def _can_ligate(frag1: Dict, frag2: Dict, tolerance: bool, dephosphorylated: List[str]) -> bool:
    """检查两个片段是否可以连接"""
    # 检查去磷酸化状态
    frag1_5p_phosphorylated = frag1.get("metadata", {}).get("phosphorylation", {}).get("5end", True)
    frag1_3p_phosphorylated = frag1.get("metadata", {}).get("phosphorylation", {}).get("3end", True)
    frag2_5p_phosphorylated = frag2.get("metadata", {}).get("phosphorylation", {}).get("5end", True)
    frag2_3p_phosphorylated = frag2.get("metadata", {}).get("phosphorylation", {}).get("3end", True)
    
    # 检查粘性末端匹配
    frag1_3end = frag1.get("overhang_3", {"kind": "blunt", "seq": ""})
    frag2_5end = frag2.get("overhang_5", {"kind": "blunt", "seq": ""})
    
    # 平末端连接
    if frag1_3end["kind"] == "blunt" and frag2_5end["kind"] == "blunt":
        return frag1_3p_phosphorylated and frag2_5p_phosphorylated
    
    # 粘性末端连接
    if frag1_3end["kind"] == "3_overhang" and frag2_5end["kind"] == "5_overhang":
        if tolerance:
            # 允许部分匹配
            min_len = min(len(frag1_3end["seq"]), len(frag2_5end["seq"]))
            match = sum(1 for a, b in zip(frag1_3end["seq"], frag2_5end["seq"]) if a == b)
            return match >= min_len * 0.8  # 80%匹配度
        else:
            # 完全匹配
            return frag1_3end["seq"] == frag2_5end["seq"]
    
    return False

def _perform_ligation(frag1: Dict, frag2: Dict) -> Dict:
    """执行实际的连接操作"""
    # 合并序列
    new_sequence = frag1["sequence"] + frag2["sequence"]
    
    # 创建连接产物
    product = {
        "id": f"ligated_{uuid.uuid4().hex[:8]}",
        "name": f"Ligated {frag1.get('name', 'fragment1')}-{frag2.get('name', 'fragment2')}",
        "sequence": new_sequence,
        "length": len(new_sequence),
        "circular": False,  # 默认为线性
        "features": frag1.get("features", []) + frag2.get("features", []),
        "overhang_5": frag1.get("overhang_5", {"kind": "blunt", "seq": "", "length": 0}),
        "overhang_3": frag2.get("overhang_3", {"kind": "blunt", "seq": "", "length": 0}),
        "metadata": {
            **frag1.get("metadata", {}),
            **frag2.get("metadata", {}),
            "ligation": {
                "parent1": frag1["id"],
                "parent2": frag2["id"],
                "timestamp": str(uuid.uuid4())
            }
        }
    }
    
    return product
