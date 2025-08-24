import sys
import os
from typing import List, Dict, Optional, Tuple
import re

# Add the parent directory to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from common_utils.sequence import *  # 假设这里导入了 Feature, Strand, get_sequence
from tools_pool.get_sequence_info import get_sequence_info

# --- 工具函数 ---
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]

def to_1_based_inclusive_from_fwd_match(m: re.Match) -> Tuple[int, int]:
    # m.start() = 0-based 包含，m.end() = 0-based 排他
    start_1b = m.start() + 1
    end_1b = m.end()  # 已为 1-based 含端点
    return start_1b, end_1b

def to_1_based_inclusive_from_rev_match(m: re.Match, L: int) -> Tuple[int, int]:
    # 在 revcomp 序列中命中区间 [i, j)（0-based、排他端点）
    # 原序列中的 1-based（含端点）坐标映射为：
    # start = L - j + 1, end = L - i
    i, j = m.start(), m.end()
    start_1b = L - j + 1
    end_1b = L - i
    return start_1b, end_1b

def clamp_and_validate(start: int, end: int, L: int) -> Optional[Tuple[int,int]]:
    if start > end:
        start, end = end, start  # 兜底：交换
    if start < 1 or end > L:
        return None
    return start, end


# --- 主函数 ---
def find_features(
    json_path: str,
    scan_builtin: bool = True,
    custom_patterns: Optional[List[Dict[str, str]]] = None
) -> List[Feature]:
    """
    Identifies functional sites (features) in a given DNA sequence based on built-in and custom patterns.

    Args:
        json_path (str): The path to a JSON file containing sequence information.
                         This file is used to retrieve the sequence and any pre-existing features.
        scan_builtin (bool, optional): If True, scan for a predefined set of common restriction sites
                                       and other features (e.g., T7 Promoter, polyA Signal). Defaults to True.
        custom_patterns (Optional[List[Dict[str, str]]], optional): A list of dictionaries, where each dictionary
                                                                     defines a custom pattern to scan for.
                                                                     Each dictionary must contain:
                                                                     - "name" (str): The name of the feature (e.g., "EcoRI_site").
                                                                     - "type" (str): The type of the feature (e.g., "restriction_site", "promoter").
                                                                     - "pattern" (str): The DNA sequence pattern to search for (e.g., "GAATTC").
                                                                     Defaults to None.

    Returns:
        List[Feature]: A list of dictionaries, where each dictionary represents a found feature.
                       Each feature dictionary conforms to the `Feature` TypedDict schema and includes:
                       - "type" (str): The type of the feature.
                       - "start" (int): The 1-based inclusive start coordinate of the feature on the original sequence.
                       - "end" (int): The 1-based inclusive end coordinate of the feature on the original sequence.
                       - "strand" (int): The strand on which the feature was found (1 for forward, -1 for reverse).
                       - "qualifiers" (Dict[str, str]): A dictionary containing additional information,
                                                         including "label" with the feature's name.
                       Pre-existing features from the input JSON are included in the returned list.
    """
    # 先拿已有的（假定已符合 Schema）
    sequence_info = get_sequence_info(json_path)
    found_features: List[Feature] = list(sequence_info.get("features", []))  # 拷贝一份，避免原地修改

    sequence, _ = get_sequence(json_path, full_length=True)
    if not sequence:
        return found_features
    sequence = sequence.upper()
    L = len(sequence)
    rev = revcomp(sequence)

    patterns_to_scan: List[Dict[str, str]] = []
    if scan_builtin:
        # 修正和去重后的内置位点库
        builtin_patterns = [
            {"name": "EcoRI",  "type": "restriction_site", "pattern": "GAATTC"},
            {"name": "BsaI",   "type": "restriction_site", "pattern": "GGTCTC"},
            {"name": "BamHI",  "type": "restriction_site", "pattern": "GGATCC"},
            {"name": "BglII",  "type": "restriction_site", "pattern": "AGATCT"},
            {"name": "HindIII","type": "restriction_site", "pattern": "AAGCTT"},
            {"name": "KpnI",   "type": "restriction_site", "pattern": "GGTACC"},
            {"name": "NcoI",   "type": "restriction_site", "pattern": "CCATGG"},
            {"name": "NdeI",   "type": "restriction_site", "pattern": "CATATG"},
            {"name": "NheI",   "type": "restriction_site", "pattern": "GCTAGC"},  # 修正
            {"name": "NotI",   "type": "restriction_site", "pattern": "GCGGCCGC"},
            {"name": "PacI",   "type": "restriction_site", "pattern": "TTAATTAA"},
            {"name": "PstI",   "type": "restriction_site", "pattern": "CTGCAG"},
            {"name": "SacI",   "type": "restriction_site", "pattern": "GAGCTC"},
            {"name": "SalI",   "type": "restriction_site", "pattern": "GTCGAC"},
            {"name": "SmaI",   "type": "restriction_site", "pattern": "CCCGGG"},
            {"name": "XbaI",   "type": "restriction_site", "pattern": "TCTAGA"},
            {"name": "XhoI",   "type": "restriction_site", "pattern": "CTCGAG"},
            {"name": "AflII",  "type": "restriction_site", "pattern": "CTTAAG"},
            {"name": "T7 Promoter", "type": "promoter",     "pattern": "TAATACGACTCACTATAGGG"},
            {"name": "polyA Signal","type": "polyA_signal", "pattern": "AATAAA"},
        ]
        patterns_to_scan.extend(builtin_patterns)

    if custom_patterns:
        # 允许外部自定义（同样按精准序列匹配，不用正则/IUPAC）
        patterns_to_scan.extend(custom_patterns)

    # 扫描
    for pinfo in patterns_to_scan:
        pat = pinfo["pattern"].upper()

        # 正链
        for m in re.finditer(re.escape(pat), sequence):
            s1, e1 = to_1_based_inclusive_from_fwd_match(m)
            rng = clamp_and_validate(s1, e1, L)
            if not rng:
                continue
            start_1b, end_1b = rng
            f: Feature = {
                "type": pinfo["type"],
                "start": int(start_1b),
                "end": int(end_1b),
                "strand": 1,
                "qualifiers": {"label": pinfo["name"]}
            }
            found_features.append(f)

        # 反链
        for m in re.finditer(re.escape(pat), rev):
            s1, e1 = to_1_based_inclusive_from_rev_match(m, L)
            rng = clamp_and_validate(s1, e1, L)
            if not rng:
                continue
            start_1b, end_1b = rng
            f: Feature = {
                "type": pinfo["type"],
                "start": int(start_1b),
                "end": int(end_1b),
                "strand": -1,
                "qualifiers": {"label": pinfo["name"]}
            }
            found_features.append(f)

    return found_features


if __name__ == "__main__":
    # 测试用例：假设有一个 JSON 文件路径
    json_path = "D:/github_repo/PPversion2/data/temp/pcDNA3.1(-).json"
    features = find_features(json_path)
    print(f"找到 {len(features)} 个功能位点：")
    print(features)
