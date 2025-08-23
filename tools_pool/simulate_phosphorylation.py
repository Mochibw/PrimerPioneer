import json
import os
import uuid
from typing import Dict, List, Optional, Union, Literal
from common_utils.sequence import SequenceRecord
from common_utils.file_operations import write_record_to_json, load_sequence_from_json


def simulate_phosphorylation(json_path: str,
                           action: Literal["phosphorylate", "dephosphorylate"],
                           ends: List[Literal["5end", "3end", "both"]] = ["both"],
                           output_path: Optional[str] = None) -> Dict[str, Union[SequenceRecord, str]]:
    """
    对DNA片段的5'端进行磷酸化或去磷酸化处理。
    
    Args:
        json_path: 输入DNA片段的JSON文件路径
        action: 操作类型，"phosphorylate"或"dephosphorylate"
        ends: 要处理的末端列表，["5end", "3end", "both"]
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含处理后的产物和状态信息的字典
    """
    # 加载DNA片段
    fragment = load_sequence_from_json(json_path)
    
    # 创建修改后的片段
    modified_fragment = fragment.copy()
    
    # 更新metadata记录磷酸化状态
    metadata = modified_fragment.get("metadata", {})
    phosphorylation_status = metadata.get("phosphorylation", {})
    
    if "both" in ends or "5end" in ends:
        phosphorylation_status["5end"] = action == "phosphorylate"
    
    if "both" in ends or "3end" in ends:
        phosphorylation_status["3end"] = action == "phosphorylate"
    
    metadata["phosphorylation"] = phosphorylation_status
    modified_fragment["metadata"] = metadata
    modified_fragment["id"] = f"{action}ed_{uuid.uuid4().hex[:8]}"
    modified_fragment["name"] = f"{action.capitalize()}d {fragment.get('name', 'fragment')}"
    
    # 生成状态消息
    processed_ends = []
    if "both" in ends or "5end" in ends:
        processed_ends.append("5'端")
    if "both" in ends or "3end" in ends:
        processed_ends.append("3'端")
    
    end_str = "和".join(processed_ends)
    message = f"{end_str}已成功{ '磷酸化' if action == 'phosphorylate' else '去磷酸化' }"
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(modified_fragment, output_path)
        message += f"。产物已保存至: {output_path}"
    
    return {
        "product": modified_fragment,
        "message": message
    }
