import json
import os
import uuid
from typing import Dict, Optional, Union
from common_utils.sequence import SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json
from langchain.tools import tool

@tool()
def simulate_a_tailing(json_path: str,
               tail_length: int = 1,
               output_path: Optional[str] = None) -> Dict[str, Union[SequenceRecord, str]]:
    """
    在DNA片段的3'端添加A突出尾。
    
    Args:
        json_path: 输入DNA片段的JSON文件路径
        tail_length: 要添加的A碱基数量
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含加尾后产物和状态信息的字典
    """
    # 加载DNA片段
    fragment = load_sequence_from_json(json_path)
    
    # 添加A尾
    a_tail = "A" * tail_length
    new_sequence = fragment["sequence"] + a_tail
    
    # 创建加尾后的片段
    tailed_fragment = {
        "id": f"a_tailed_{uuid.uuid4().hex[:8]}",
        "name": f"A-tailed {fragment.get('name', 'fragment')}",
        "sequence": new_sequence,
        "length": len(new_sequence),
        "circular": fragment.get("circular", False),
        "features": fragment.get("features", []),
        "overhang_5": fragment.get("overhang_5", {"kind": "blunt", "seq": "", "length": 0}),
        "overhang_3": {
            "kind": "3_overhang",
            "seq": a_tail,
            "length": tail_length
        },
        "metadata": {
            **fragment.get("metadata", {}),
            "a_tailing": {
                "tail_length": tail_length,
                "added_sequence": a_tail
            }
        }
    }
    
    message = f"成功在3'端添加了{tail_length}个A碱基的突出尾"
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(tailed_fragment, output_path)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "product": tailed_fragment,
        "message": message
    }

def simulate_t_overhang(json_path: str,
                   tail_length: int = 1,
                   output_path: Optional[str] = None) -> Dict[str, Union[SequenceRecord, str]]:
    """
    在DNA片段的3'端添加T突出尾（用于TA克隆）。
    
    Args:
        json_path: 输入DNA片段的JSON文件路径
        tail_length: 要添加的T碱基数量
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含加尾后产物和状态信息的字典
    """
    # 加载DNA片段
    fragment = load_sequence_from_json(json_path)
    
    # 添加T尾
    t_tail = "T" * tail_length
    new_sequence = fragment["sequence"] + t_tail
    
    # 创建加尾后的片段
    tailed_fragment = {
        "id": f"t_tailed_{uuid.uuid4().hex[:8]}",
        "name": f"T-tailed {fragment.get('name', 'fragment')}",
        "sequence": new_sequence,
        "length": len(new_sequence),
        "circular": fragment.get("circular", False),
        "features": fragment.get("features", []),
        "overhang_5": fragment.get("overhang_5", {"kind": "blunt", "seq": "", "length": 0}),
        "overhang_3": {
            "kind": "3_overhang",
            "seq": t_tail,
            "length": tail_length
        },
        "metadata": {
            **fragment.get("metadata", {}),
            "t_tailing": {
                "tail_length": tail_length,
                "added_sequence": t_tail
            }
        }
    }
    
    message = f"成功在3'端添加了{tail_length}个T碱基的突出尾"
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(tailed_fragment, output_path)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "product": tailed_fragment,
        "message": message
    }
