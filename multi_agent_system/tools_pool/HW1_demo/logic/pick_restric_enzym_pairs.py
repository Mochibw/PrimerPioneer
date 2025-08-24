# pick_restric_enzyme_pairs.py
from typing import List, Tuple, Dict
from snapgene_reader import snapgene_file_to_dict
from Bio.Restriction import RestrictionBatch, CommOnly
from pathlib import Path
import json
from Bio.Seq import Seq

def _norm_dna(s: str) -> str:
    return "".join(ch for ch in s.upper() if ch in "ATGC")

def _get_mcs_range(features: list[dict], name: str = "MCS") -> tuple[int, int]:
    """
    在 features 中精确找到名为 'MCS' 的特征，返回 (start, end) 半开区间（0-based）。
    说明：按你这份 .dna，start/end 就在顶层；不处理跨 origin。
    """
    for f in features:
        if f.get("name") == name:
            # 这里的 start/end 已经是整数（snapgene_reader 输出）
            start = int(f["start"])
            end = int(f["end"])
            # 若你希望用 1-based 的 @range，也可解析："882-1026" → (881, 1026)
            return start, end
    raise ValueError("未找到名为 'MCS' 的 feature，请在 .dna 中精确标注。")

def _scan_unique_sites(seq: str, mcs_start: int, mcs_end: int) -> list[dict]:
    rb = RestrictionBatch(list(CommOnly))
    seq_obj = Seq(seq)              # ← 把 str 包成 Biopython 的 Seq
    hits = rb.search(seq_obj)       # ← 用 Seq 对象搜索
    keep = []
    for enz, pos_list in hits.items():
        if not pos_list or len(pos_list) != 1:
            continue
        pos0 = pos_list[0] - 1
        if mcs_start <= pos0 < mcs_end:
            keep.append({"name": enz.__name__, "site": enz.site, "pos0": pos0})
    keep.sort(key=lambda x: x["pos0"])
    return keep

def _filter_by_insert(insert_seq: str, enzymes: List[Dict]) -> List[Dict]:
    ins = _norm_dna(insert_seq)
    return [e for e in enzymes if _norm_dna(e["site"]) not in ins]

def pick_enzyme_pairs_from_dna(dna_path: str, insert_seq: str) -> List[Dict]:
    d = snapgene_file_to_dict(dna_path)
    full_seq = _norm_dna(d["seq"])
    mcs_start, mcs_end = _get_mcs_range(d["features"], "MCS")
    enzymes = _scan_unique_sites(full_seq, mcs_start, mcs_end)
    enzymes = _filter_by_insert(insert_seq, enzymes)
    return enzymes

# # ---------------- 调试入口 ----------------
# if __name__ == "__main__":
#     # ====== 配置区 ======
#     dna_path = "data/pcDNA3.1(-).dna"
#     insert_path = "data/temp_cds/RTCB_Homo_sapiens_lcl_NM_014306.5_cds_NP_055121.1_1.fasta"
#     # ====================

#     # 读取 insert 文件内容（支持 FASTA / 纯序列）
#     insert_text = Path(insert_path).read_text()
#     insert_seq = "".join(line.strip() for line in insert_text.splitlines() if not line.startswith(">"))

#     enzymes = pick_enzyme_pairs_from_dna(dna_path, insert_seq)
#     print("\n通过的酶（按 MCS 顺序）:")
#     for e in enzymes:
#         print(f"{e['name']} (Site: {e['site']}, Pos: {e['pos0']})")
