import json
import os
from typing import Dict, Literal, List, Tuple, Union
from common_utils.sequence_tools import (
    check_pcr_feasibility,
    has_high_repetitiveness,
    check_homopolymers,
    calculate_gc_content
)
from common_utils.sequence import get_sequence

def select_cloning_method(json_path: str,
                         purpose: Literal["general", "expression", "mutagenesis", "library"] = "general",
                         budget_constraints: bool = False,
                         time_constraints: bool = False) -> Dict:
    """
    根据序列特性选择最适合的克隆方法（PCR或化学合成）。
    
    Args:
        json_path: 包含序列的JSON文件路径
        purpose: 克隆目的
        budget_constraints: 是否有预算限制
        time_constraints: 是否有时间限制
        
    Returns:
        包含方法选择建议和详细分析结果的字典
    """
    # 获取序列
    sequence, length = get_sequence(json_path, full_length=True)
    if not sequence:
        return {"error": "无法读取序列"}
    
    # 执行详细分析
    analysis = _analyze_sequence_for_cloning(sequence, length)
    
    # 根据分析结果选择方法
    recommended_method, confidence, reasoning = _determine_best_method(
        analysis, purpose, budget_constraints, time_constraints
    )
    
    # 生成下一步操作建议
    next_steps = _generate_next_steps(recommended_method, json_path, purpose)
    
    return {
        "sequence_analysis": analysis,
        "recommended_method": recommended_method,
        "confidence": confidence,
        "reasoning": reasoning,
        "next_steps": next_steps,
        "considerations": {
            "purpose": purpose,
            "budget_constraints": budget_constraints,
            "time_constraints": time_constraints
        }
    }

def _analyze_sequence_for_cloning(sequence: str, length: int) -> Dict:
    """对序列进行详细的克隆适宜性分析"""
    # PCR可行性分析
    pcr_analysis = check_pcr_feasibility(sequence)
    
    # 额外分析指标
    gc_content = calculate_gc_content(sequence)
    repetitiveness = has_high_repetitiveness(sequence)
    homopolymers = check_homopolymers(sequence, max_length=6)
    
    # 计算综合评分
    pcr_score = _calculate_pcr_score(pcr_analysis, length, gc_content, repetitiveness)
    synthesis_score = _calculate_synthesis_score(pcr_analysis, length, gc_content)
    
    return {
        "length": length,
        "gc_content": gc_content,
        "pcr_feasibility": pcr_analysis,
        "repetitiveness": repetitiveness,
        "homopolymers": homopolymers,
        "scores": {
            "pcr_suitability": pcr_score,
            "synthesis_suitability": synthesis_score,
            "overall_recommendation": "pcr" if pcr_score > synthesis_score else "synthesis"
        },
        "critical_issues": _identify_critical_issues(pcr_analysis, repetitiveness, homopolymers, length)
    }

def _calculate_pcr_score(analysis: Dict, length: int, gc_content: float, repetitiveness: bool) -> float:
    """计算PCR适宜性评分"""
    score = 0.0
    
    # 长度因素（0-40分）
    if length <= 1000:
        score += 40
    elif length <= 2000:
        score += 30
    elif length <= 3000:
        score += 20
    else:
        score += 10
    
    # GC含量因素（0-30分）
    if 0.4 <= gc_content <= 0.6:
        score += 30
    elif 0.3 <= gc_content <= 0.7:
        score += 20
    else:
        score += 5
    
    # 重复性因素（0-20分）
    if not repetitiveness:
        score += 20
    else:
        score += 5
    
    # 问题数量因素（0-10分）
    issue_count = len(analysis.get("issues", []))
    score += max(0, 10 - issue_count * 2)
    
    return min(100, score)

def _calculate_synthesis_score(analysis: Dict, length: int, gc_content: float) -> float:
    """计算化学合成适宜性评分"""
    score = 0.0
    
    # 长度因素（合成擅长长序列）
    if length > 2000:
        score += 40
    elif length > 1000:
        score += 30
    else:
        score += 10
    
    # GC含量因素（合成可以处理极端GC）
    if gc_content < 0.3 or gc_content > 0.7:
        score += 30
    else:
        score += 20
    
    # 问题数量因素（合成可以解决PCR问题）
    issue_count = len(analysis.get("issues", []))
    score += min(30, issue_count * 5)
    
    return min(100, score)

def _identify_critical_issues(pcr_analysis: Dict, repetitiveness: bool, 
                             homopolymers: Dict, length: int) -> List[str]:
    """识别关键问题"""
    issues = []
    
    # PCR可行性问题
    if not pcr_analysis["feasible"]:
        issues.extend(pcr_analysis["issues"])
    
    # 重复性问题
    if repetitiveness:
        issues.append("高度重复序列")
    
    # 同聚物问题
    if homopolymers["has_long_homopolymer"]:
        issues.append(f"存在长同聚物({homopolymers['longest_homopolymer']}bp)")
    
    # 长度问题
    if length > 3000:
        issues.append("序列过长(>3000bp)")
    
    return issues

def _determine_best_method(analysis: Dict, purpose: str, 
                          budget_constraints: bool, time_constraints: bool) -> tuple:
    """根据分析结果确定最佳方法"""
    pcr_score = analysis["scores"]["pcr_suitability"]
    synthesis_score = analysis["scores"]["synthesis_suitability"]
    
    # 基础决策
    if pcr_score > synthesis_score + 20:
        base_method = "pcr"
        confidence = "high"
        reasoning = "序列非常适合PCR扩增"
    elif synthesis_score > pcr_score + 20:
        base_method = "synthesis"
        confidence = "high"
        reasoning = "序列更适合化学合成"
    else:
        base_method = "pcr" if pcr_score >= synthesis_score else "synthesis"
        confidence = "medium"
        reasoning = "两种方法都可行，基于评分选择"
    
    # 考虑应用目的
    if purpose == "mutagenesis" and len(analysis["critical_issues"]) == 0:
        base_method = "pcr"
        reasoning += "。定点突变通常使用PCR方法"
    elif purpose == "library" and analysis["length"] < 1500:
        base_method = "pcr"
        reasoning += "。文库构建通常使用PCR方法"
    
    # 考虑约束条件
    if budget_constraints and base_method == "synthesis":
        base_method = "pcr"
        reasoning += "。由于预算限制，选择更经济的PCR方法"
        confidence = "low"
    
    if time_constraints and base_method == "synthesis":
        base_method = "pcr"
        reasoning += "。由于时间限制，选择更快速的PCR方法"
        confidence = "medium"
    
    return base_method, confidence, reasoning

def _generate_next_steps(method: str, json_path: str, purpose: str) -> List[Dict]:
    """生成下一步操作建议"""
    if method == "pcr":
        return [
            {
                "action": "design_primers",
                "description": "设计PCR引物",
                "tool": "design_primer_suite",
                "parameters": {
                    "target": json_path,
                    "task": "amplify",
                    "intent": {"amplicon": {"start": 1, "end": "full_length"}}
                }
            },
            {
                "action": "simulate_pcr",
                "description": "模拟PCR扩增",
                "tool": "simulate_pcr",
                "parameters": {
                    "json_path": json_path,
                    "forward": "待设计的前引物",
                    "reverse": "待设计的后引物"
                }
            }
        ]
    else:  # synthesis
        return [
            {
                "action": "simulate_dna_synthesis",
                "description": "模拟DNA化学合成",
                "tool": "simulate_dna_synthesis",
                "parameters": {
                    "sequence": "从文件中提取序列",
                    "circular": False,
                    "provider_constraints": {
                        "max_length": 3000,
                        "min_gc": 0.2,
                        "max_gc": 0.8
                    }
                }
            }
        ]

def get_method_recommendation(json_path: str) -> Dict:
    """
    快速获取方法推荐（简化版接口，供LLM直接调用）
    
    Args:
        json_path: 包含序列的JSON文件路径
        
    Returns:
        简化的方法推荐结果
    """
    full_analysis = select_cloning_method(json_path)
    
    return {
        "recommended_method": full_analysis["recommended_method"],
        "confidence": full_analysis["confidence"],
        "reasoning": full_analysis["reasoning"],
        "next_action": full_analysis["next_steps"][0]["action"] if full_analysis["next_steps"] else "none",
        "critical_issues": full_analysis["sequence_analysis"]["critical_issues"]
    }

# 示例使用
if __name__ == "__main__":
    # 测试示例
    result = select_cloning_method(
        "path/to/sequence.json",
        purpose="expression",
        budget_constraints=False,
        time_constraints=True
    )
    
    print("推荐方法:", result["recommended_method"])
    print("置信度:", result["confidence"])
    print("理由:", result["reasoning"])
    print("下一步:", result["next_steps"][0]["action"])
