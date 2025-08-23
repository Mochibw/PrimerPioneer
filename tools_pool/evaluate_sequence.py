import json
import os
from typing import Dict
from common_utils.sequence_tools import (
    check_pcr_feasibility, 
    has_high_repetitiveness,
    check_homopolymers,
    predict_secondary_structure,
    find_restriction_sites,
    calculate_molecular_weight
)
from common_utils.sequence import get_sequence
from tools_pool.select_cloning_method import get_method_recommendation

def evaluate_sequence(json_path: str) -> Dict:
    """
    综合评估DNA序列的克隆适宜性。
    
    Args:
        json_path: 包含序列的JSON文件路径
        
    Returns:
        包含评估结果的字典
    """
    sequence, length = get_sequence(json_path, full_length=True)
    
    if not sequence:
        return {"error": "无法读取序列"}
    
    # 执行各种评估
    pcr_feasibility = check_pcr_feasibility(sequence)
    repetitiveness = has_high_repetitiveness(sequence)
    homopolymers = check_homopolymers(sequence)
    secondary_structure = predict_secondary_structure(sequence)
    restriction_sites = find_restriction_sites(sequence)
    molecular_weight = calculate_molecular_weight(sequence)
    
    # 获取方法推荐
    method_recommendation = get_method_recommendation(json_path)
    
    # 总体建议
    recommendation = "适合PCR克隆"
    if not pcr_feasibility["feasible"]:
        if any(issue in pcr_feasibility["recommendation"] for issue in ["GC含量过高", "高度重复", "序列过长"]):
            recommendation = "建议使用DNA化学合成"
        else:
            recommendation = "需要优化PCR条件"
    
    return {
        "sequence_info": {
            "length": length,
            "gc_content": pcr_feasibility["gc_content"],
            "molecular_weight": molecular_weight
        },
        "pcr_feasibility": pcr_feasibility,
        "repetitiveness": {
            "is_highly_repetitive": repetitiveness,
            "entropy_ratio": _calculate_entropy_ratio(sequence)
        },
        "homopolymers": homopolymers,
        "secondary_structure": secondary_structure,
        "restriction_sites": restriction_sites,
        "method_recommendation": method_recommendation,
        "overall_recommendation": recommendation,
        "synthesis_recommended": not pcr_feasibility["feasible"] and 
                               ("化学合成" in pcr_feasibility["recommendation"])
    }

def _calculate_entropy_ratio(sequence: str) -> float:
    """计算序列熵与最大熵的比值"""
    from collections import Counter
    counts = Counter(sequence.upper())
    total = len(sequence)
    
    entropy = 0
    for count in counts.values():
        p = count / total
        entropy -= p * math.log2(p) if p > 0 else 0
    
    max_entropy = math.log2(min(4, len(counts)))
    return entropy / max_entropy if max_entropy > 0 else 0
