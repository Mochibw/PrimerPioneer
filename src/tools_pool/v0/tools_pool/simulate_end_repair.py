import json
import os
import uuid
from typing import Dict, Optional, Union
from common_utils.sequence import SequenceRecord, write_record_to_json, load_sequence_from_json

def simulate_end_repair(json_path: str, 
                output_path: Optional[str] = None) -> Dict[str, Union[SequenceRecord, str]]:
    """
    将DNA片段的粘性末端修复为平末端。
    
    Args:
        json_path: 输入DNA片段的JSON文件路径
        output_path: 输出JSON文件路径(可选)
        
    Returns:
        包含修复后产物和状态信息的字典
    """
    # 加载DNA片段
    fragment = load_sequence_from_json(json_path)
    
    # 检查是否需要修复
    needs_repair = False
    five_end = fragment.get("overhang_5", {"kind": "blunt", "seq": ""})
    three_end = fragment.get("overhang_3", {"kind": "blunt", "seq": ""})
    
    if five_end["kind"] != "blunt" or three_end["kind"] != "blunt":
        needs_repair = True
    
    if needs_repair:
        # 执行末端修复
        repaired_seq = _perform_end_repair(fragment)
        
        repaired_fragment = {
            "id": f"repaired_{uuid.uuid4().hex[:8]}",
            "name": f"End-repaired {fragment.get('name', 'fragment')}",
            "sequence": repaired_seq,
            "length": len(repaired_seq),
            "circular": fragment.get("circular", False),
            "features": fragment.get("features", []),
            "overhang_5": {"kind": "blunt", "seq": "", "length": 0},
            "overhang_3": {"kind": "blunt", "seq": "", "length": 0},
            "metadata": {
                **fragment.get("metadata", {}),
                "end_repair": "completed",
                "original_5overhang": five_end.get("seq", ""),
                "original_3overhang": three_end.get("seq", "")
            }
        }
        
        message = "末端修复完成，粘性末端已转为平末端"
    else:
        repaired_fragment = fragment
        message = "片段已是平末端，无需修复"
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(repaired_fragment, output_path)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "product": repaired_fragment,
        "message": message
    }

def _perform_end_repair(fragment: Dict) -> str:
    """执行实际的末端修复操作"""
    sequence = fragment["sequence"]
    five_end = fragment.get("overhang_5", {"kind": "blunt", "seq": ""})
    three_end = fragment.get("overhang_3", {"kind": "blunt", "seq": ""})
    
    # 处理5'突出端
    if five_end["kind"] == "5_overhang":
        overhang_seq = five_end.get("seq", "")
        # 移除5'突出部分（保留主体序列）
        sequence = sequence[len(overhang_seq):]
    
    # 处理3'突出端
    if three_end["kind"] == "3_overhang":
        overhang_seq = three_end.get("seq", "")
        # 移除3'突出部分
        sequence = sequence[:-len(overhang_seq)]
    
    # 如果是5'凹陷端（3'突出），需要填充
    if five_end["kind"] == "3_overhang":
        overhang_seq = five_end.get("seq", "")
        # 这里简化处理，实际酶切修复会填充碱基
        pass
        
    # 如果是3'凹陷端（5'突出），需要填充
    if three_end["kind"] == "5_overhang":
        overhang_seq = three_end.get("seq", "")
        # 这里简化处理，实际酶切修复会填充碱基
        pass
    
    return sequence
