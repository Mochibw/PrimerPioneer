import sys
import os
import json
from typing import List, Dict, Optional, Union
from common_utils.sequence import Fragment
from common_utils.file_operations import load_sequence_from_json, write_record_to_json

def simulate_gel_purification(json_path: str,
                             selected_indices: List[int],
                             output_path: Optional[str] = None) -> Dict[str, Union[List[str], List[int], List[Fragment], List[str]]]:
    """
    模拟胶回收过程，从酶切产物中选择特定片段。
    
    Args:
        json_path: 输入的酶切产物JSON文件路径（由simulate_restriction_digestion生成）
        selected_indices: 要保留的片段索引列表（0-based）
        output_path: 输出JSON文件路径（可选）
        
    Returns:
        包含选定片段和状态信息的字典
    """
    # 加载酶切产物
    digest_result = load_sequence_from_json(json_path)
    
    # 检查输入数据格式
    if "fragments" not in digest_result:
        raise ValueError("输入文件不是有效的酶切产物格式，缺少'fragments'字段")
    
    all_fragments = digest_result["fragments"]
    
    # 验证索引
    valid_indices = []
    invalid_indices = []
    
    for idx in selected_indices:
        if 0 <= idx < len(all_fragments):
            valid_indices.append(idx)
        else:
            invalid_indices.append(idx)
    
    # 选择指定的片段
    selected_fragments = [all_fragments[i] for i in valid_indices]
    
    # 生成状态消息
    info = digest_result.get("info", [])
    info.append(f"从 {len(all_fragments)} 个片段中选择了 {len(selected_fragments)} 个片段")
    
    if invalid_indices:
        info.append(f"警告：以下索引无效并被忽略: {invalid_indices}")
    
    # 如果有酶切信息，也一并保留
    enzymes = digest_result.get("enzymes", [])
    cuts = digest_result.get("cuts", [])
    
    # 构建输出数据（与DigestResult结构相同）
    output_data = {
        "enzymes": enzymes,
        "cuts": cuts,
        "fragments": selected_fragments,
        "info": info
    }
    
    # 保存结果
    if output_path:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        write_record_to_json(output_data, output_path)
        info.append(f"选定片段已保存至: {output_path}")
    
    return output_data

# Example usage (for testing)
if __name__ == "__main__":
    # 示例：选择第1个和第3个片段（索引0和2）
    result = simulate_gel_purification(
        json_path="data/temp/digest_result.json",
        selected_indices=[0, 2],
        output_path="data/temp/purified_fragments.json"
    )
    # 打印info消息
    for msg in result["info"]:
        print(msg)
