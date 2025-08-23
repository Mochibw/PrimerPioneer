import json
import os
import uuid
from typing import Dict, Optional, Union
from common_utils.sequence import SequenceRecord, write_record_to_json
from common_utils.sequence_tools import calculate_gc_content

def simulate_dna_synthesis(sequence: str,
                         circular: bool = False,
                         provider_constraints: Optional[Dict[str, Union[int, float, str, bool]]] = None,
                         output_path: Optional[str] = None) -> Dict[str, Union[SequenceRecord, str]]:
    """
    模拟DNA化学合成过程。
    
    Args:
        sequence: 要合成的DNA序列
        circular: 是否为环状DNA
        provider_constraints: 合成提供商限制条件
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含合成产物和状态信息的字典
    """
    # 标准化序列
    seq = sequence.upper().strip()
    
    # 检查序列是否适合合成
    synthesis_result = _check_synthesis_feasibility(seq, provider_constraints)
    
    if not synthesis_result["feasible"]:
        return {
            "product": None,
            "message": f"DNA合成不可行: {synthesis_result['reason']}",
            "recommendation": synthesis_result.get("recommendation", "")
        }
    
    # 创建合成产物
    synthesized_dna = {
        "id": f"synthesized_{uuid.uuid4().hex[:8]}",
        "name": f"Synthesized DNA ({len(seq)} bp)",
        "sequence": seq,
        "length": len(seq),
        "circular": circular,
        "features": [],
        "overhang_5": {"kind": "blunt", "seq": "", "length": 0},
        "overhang_3": {"kind": "blunt", "seq": "", "length": 0},
        "metadata": {
            "synthesis_method": "chemical",
            "synthesis_status": "success",
            "provider_constraints": provider_constraints or {},
            "gc_content": _calculate_gc_content(seq),
            "sequence_complexity": _calculate_sequence_complexity(seq)
        }
    }
    
    message = "DNA化学合成成功完成"
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(synthesized_dna, output_path)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "product": synthesized_dna,
        "message": message
    }

def _check_synthesis_feasibility(sequence: str, constraints: Optional[Dict] = None) -> Dict:
    """检查DNA序列是否适合化学合成"""
    constraints = constraints or {}
    max_length = constraints.get("max_length", 3000)
    min_gc = constraints.get("min_gc", 0.2)
    max_gc = constraints.get("max_gc", 0.8)
    max_homopolymer = constraints.get("max_homopolymer", 6)
    
    # 检查长度
    if len(sequence) > max_length:
        return {
            "feasible": False,
            "reason": f"序列长度({len(sequence)} bp)超过最大限制({max_length} bp)",
            "recommendation": "考虑分段合成或使用其他方法"
        }
    
    # 检查GC含量
    gc_content = calculate_gc_content(sequence)
    if gc_content < min_gc or gc_content > max_gc:
        return {
            "feasible": False,
            "reason": f"GC含量({gc_content:.2%})超出可接受范围({min_gc:.0%}-{max_gc:.0%})",
            "recommendation": "重新设计序列或选择专业合成服务"
        }
    
    # 检查同聚物长度
    homopolymer_info = _check_homopolymers(sequence, max_homopolymer)
    if homopolymer_info["has_long_homopolymer"]:
        return {
            "feasible": False,
            "reason": f"检测到过长同聚物: {homopolymer_info['longest_homopolymer']}",
            "recommendation": "修改序列以避免长同聚物"
        }
    
    # 检查重复序列
    if _has_high_repetitiveness(sequence):
        return {
            "feasible": False,
            "reason": "序列包含高度重复区域",
            "recommendation": "重新设计序列或选择特殊合成服务"
        }
    
    return {"feasible": True}

# def _calculate_gc_content(seq: str) -> float:
#     """计算GC含量"""
#     gc_count = seq.count('G') + seq.count('C')
#     return gc_count / len(seq) if seq else 0

def _check_homopolymers(seq: str, max_length: int) -> Dict:
    """检查同聚物长度"""
    current_base = None
    current_length = 0
    max_homopolymer = 0
    
    for base in seq:
        if base == current_base:
            current_length += 1
            max_homopolymer = max(max_homopolymer, current_length)
        else:
            current_base = base
            current_length = 1
    
    return {
        "has_long_homopolymer": max_homopolymer > max_length,
        "longest_homopolymer": max_homopolymer
    }

def _calculate_sequence_complexity(seq: str) -> float:
    """计算序列复杂度"""
    if len(seq) <= 1:
        return 1.0
    
    unique_kmers = set()
    k = min(4, len(seq) - 1)  # 使用3-mer或4-mer
    
    for i in range(len(seq) - k + 1):
        unique_kmers.add(seq[i:i+k])
    
    return len(unique_kmers) / (len(seq) - k + 1)

def _has_high_repetitiveness(seq: str, threshold: float = 0.7) -> bool:
    """检查序列是否高度重复"""
    if len(seq) < 20:
        return False
    
    # 简单重复检查：寻找重复模式
    for pattern_length in range(3, min(10, len(seq)//2)):
        for i in range(len(seq) - pattern_length * 2 + 1):
            pattern = seq[i:i+pattern_length]
            count = 1
            j = i + pattern_length
            while j <= len(seq) - pattern_length:
                if seq[j:j+pattern_length] == pattern:
                    count += 1
                    j += pattern_length
                else:
                    break
            
            if count >= 3:  # 重复3次以上
                repetition_density = (count * pattern_length) / len(seq)
                if repetition_density > threshold:
                    return True
    
    return False
