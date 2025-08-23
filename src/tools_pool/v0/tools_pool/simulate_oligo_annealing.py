import json
import os
import uuid
from typing import Dict, Optional, Union
from common_utils.sequence import SequenceRecord, write_record_to_json

def simulate_oligo_annealing(oligo1_seq: str, 
                  oligo2_seq: str, 
                  output_path: Optional[str] = None) -> Dict[str, Union[SequenceRecord, str]]:
    """
    将两条互补的单链寡核苷酸退火成双链DNA片段。
    
    Args:
        oligo1_seq: 第一条寡核苷酸序列
        oligo2_seq: 第二条寡核苷酸序列
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含退火产物和状态信息的字典
    """
    # 标准化序列
    oligo1 = oligo1_seq.upper().strip()
    oligo2 = oligo2_seq.upper().strip()
    
    # 检查互补性
    complement = _get_reverse_complement(oligo2)
    
    if oligo1 == complement:
        # 完全互补，可以退火
        ds_dna = {
            "id": f"annealed_{uuid.uuid4().hex[:8]}",
            "name": "Annealed oligo duplex",
            "sequence": oligo1,  # 使用其中一条链作为代表序列
            "length": len(oligo1),
            "circular": False,
            "features": [],
            "overhang_5": {"kind": "blunt", "seq": "", "length": 0},
            "overhang_3": {"kind": "blunt", "seq": "", "length": 0},
            "metadata": {
                "oligo1": oligo1,
                "oligo2": oligo2,
                "type": "double_stranded_oligo"
            }
        }
        
        message = "寡核苷酸成功退火形成双链DNA"
    else:
        # 不完全互补
        ds_dna = None
        message = f"寡核苷酸不完全互补，无法形成稳定双链。匹配度: {_calculate_complementarity(oligo1, oligo2):.1f}%"
    
    # 保存结果
    if output_path and ds_dna:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(ds_dna, output_path)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "product": ds_dna,
        "message": message
    }

def _get_reverse_complement(seq: str) -> str:
    """获取序列的反向互补序列"""
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N'}
    return ''.join(complement.get(base, 'N') for base in reversed(seq))

def _calculate_complementarity(seq1: str, seq2: str) -> float:
    """计算两条序列的互补性百分比"""
    min_len = min(len(seq1), len(seq2))
    if min_len == 0:
        return 0.0
    
    matches = 0
    for i in range(min_len):
        comp = _get_reverse_complement(seq2[i])
        if seq1[i] == comp:
            matches += 1
    
    return (matches / min_len) * 100
