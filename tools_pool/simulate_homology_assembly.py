import json
import os
import uuid
from typing import Dict, List, Optional, Union, Literal
from common_utils.sequence import SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json


def simulate_homology_assembly(fragments_paths: List[str],
                             overlap_min: int = 20,
                             method: Literal["gibson", "isothermal", "SLiCE", "CPEC"] = "gibson",
                             output_path: Optional[str] = None) -> Dict[str, Union[List[SequenceRecord], str]]:
    """
    模拟基于同源臂的无缝组装。
    
    Args:
        fragments_paths: 要组装的DNA片段JSON文件路径列表
        overlap_min: 最小同源臂长度
        method: 组装方法
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含组装产物和状态信息的字典
    """
    # 加载所有片段
    fragments = []
    for path in fragments_paths:
        fragment = load_sequence_from_json(path)
        fragments.append(fragment)
    
    if len(fragments) < 2:
        return {
            "products": fragments,
            "message": "需要至少2个片段进行组装"
        }
    
    # 检查同源臂并组装
    assembly_results = []
    successful_assemblies = 0
    
    for i in range(len(fragments) - 1):
        frag1 = fragments[i]
        frag2 = fragments[i + 1]
        
        # 检查同源臂
        overlap_info = _find_homology_overlap(frag1, frag2, overlap_min)
        
        if overlap_info["has_overlap"]:
            # 执行组装
            assembled = _perform_assembly(frag1, frag2, overlap_info, method)
            assembly_results.append(assembled)
            successful_assemblies += 1
        else:
            assembly_results.append({
                "status": "failed",
                "reason": f"片段 {i} 和 {i+1} 之间缺乏足够同源臂",
                "required_overlap": overlap_min,
                "found_overlap": overlap_info["overlap_length"]
            })
    
    message = f"同源重组组装完成，成功组装 {successful_assemblies} 个连接"
    
    # 保存结果
    if output_path and successful_assemblies > 0:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        if successful_assemblies == 1:
            write_record_to_json(assembly_results[0], output_path)
        else:
            with open(output_path, 'w') as f:
                json.dump({"assembly_results": assembly_results}, f, indent=2)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "assembly_results": assembly_results,
        "message": message
    }

def _find_homology_overlap(frag1: Dict, frag2: Dict, min_overlap: int) -> Dict:
    """查找两个片段之间的同源臂"""
    seq1 = frag1["sequence"]
    seq2 = frag2["sequence"]
    
    # 检查frag1的3'端与frag2的5'端的重叠
    max_overlap = min(len(seq1), len(seq2), 50)  # 限制最大检查长度
    best_overlap = 0
    
    for overlap_len in range(min_overlap, max_overlap + 1):
        end_of_frag1 = seq1[-overlap_len:]
        start_of_frag2 = seq2[:overlap_len]
        
        if end_of_frag1 == start_of_frag2:
            best_overlap = overlap_len
    
    return {
        "has_overlap": best_overlap >= min_overlap,
        "overlap_length": best_overlap,
        "overlap_sequence": seq1[-best_overlap:] if best_overlap > 0 else ""
    }

def _perform_assembly(frag1: Dict, frag2: Dict, overlap_info: Dict, method: str) -> Dict:
    """执行实际的组装操作"""
    seq1 = frag1["sequence"]
    seq2 = frag2["sequence"]
    overlap_len = overlap_info["overlap_length"]
    
    # 合并序列（去除重叠部分）
    assembled_seq = seq1 + seq2[overlap_len:]
    
    # 创建组装产物
    assembled_product = {
        "id": f"assembled_{uuid.uuid4().hex[:8]}",
        "name": f"Assembled {frag1.get('name', 'fragment1')}-{frag2.get('name', 'fragment2')}",
        "sequence": assembled_seq,
        "length": len(assembled_seq),
        "circular": False,
        "features": frag1.get("features", []) + frag2.get("features", []),
        "overhang_5": frag1.get("overhang_5", {"kind": "blunt", "seq": "", "length": 0}),
        "overhang_3": frag2.get("overhang_3", {"kind": "blunt", "seq": "", "length": 0}),
        "metadata": {
            **frag1.get("metadata", {}),
            **frag2.get("metadata", {}),
            "assembly_method": method,
            "overlap_length": overlap_len,
            "overlap_sequence": overlap_info["overlap_sequence"],
            "assembly_timestamp": str(uuid.uuid4())
        }
    }
    
    return assembled_product
