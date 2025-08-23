import sys
import os
from typing import List, Dict, Optional, Tuple
import re
import uuid
from typing_extensions import TypedDict, Literal, Union

# Add the parent directory to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from common_utils.sequence import SequenceRecord, Fragment, FragmentEnd, EndType
from common_utils.file_operations import load_sequence_from_json, write_record_to_json

# ---- 本文件内部使用的规格 ----
class EnzymeSpec(TypedDict, total=False):
    name: str            # 酶名，如 "EcoRI"
    site: str            # 识别位点（仅 A/C/G/T/N；本实现精确匹配，不支持 IUPAC 模糊）
    top_cut: int         # 正义链相对位点的切割索引（0..len(site)）
    bottom_cut: int      # 反义链相对位点的切割索引（0..len(site)）

class DigestResult(TypedDict, total=False):
    enzymes: List[str]
    cuts: List[int]             # 1-based、表示切口“左侧碱基”的索引（k 表示在 k 与 k+1 之间断裂；环状时 k 可为 L 表示末尾与首位之间）
    fragments: List[Fragment]
    info: List[str]             # 关于酶切效率的警告信息

# ---- 小工具 ----
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def _revcomp(s: str) -> str:
    return s.translate(_RC)[::-1]

def _include_seq_for_len(L: int, max_len: int = 200) -> bool:
#    return L <= max_len 这个原来的作用是太长的序列就不写入json，后来发现跟simulate_ligation冲突了，所以改成总是写入
    return True
def _normalize_record(seqrec: SequenceRecord) -> Tuple[str, int, bool]:
    s = (seqrec.get("sequence") or "").upper()
    L = len(s)
    circ = bool(seqrec.get("circular", False))
    return s, L, circ

def _calc_overhang(site: str, top_cut: int, bottom_cut: int) -> Tuple[EndType, str]:
    """
    返回（端类型, 5'端片段上的突出端序列seg）。
    约定：当 top_cut < bottom_cut 时为 5' 黏端，其“右侧片段”的 5' 端突出序列为 site[top_cut:bottom_cut]；
         当 top_cut > bottom_cut 时为 3' 黏端，其“左侧片段”的 3' 端突出序列为 site[bottom_cut:top_cut]；
         相等为平端。
    """
    if top_cut == bottom_cut:
        return "blunt", ""
    if top_cut < bottom_cut:
        seg = site[top_cut:bottom_cut]
        return "5_overhang", seg
    else:
        seg = site[bottom_cut:top_cut]
        return "3_overhang", seg

def _cut_to_boundary(start0: int, top_cut: int, L: int) -> int:
    """
    将切割（以 top_cut 的位置作为“边界”定位）映射为 0-based 边界索引 b，表示在碱基 b 与 b+1 之间断裂；
    返回 b 的 0..L-1 （环状取模，线性只会落在 0..L-2）。
    """
    b = (start0 + top_cut - 1 + L) % L
    return b

def _apply_rev_orientation(m_start_rc: int, site_len: int, L: int) -> int:
    """
    RC 序列中的 0-based 起点映射回原序列的 0-based 起点
    原串起点 = L - (m_start_rc + site_len)
    """
    return L - (m_start_rc + site_len)

def _frag_seq(seq: str, a0: int, b0: int, circ: bool) -> str:
    """返回片段序列（0-based 闭区间 a0..b0，环状允许 a0>b0 表示跨零）"""
    if not circ or a0 <= b0:
        return seq[a0:b0+1]
    else:
        return seq[a0:] + seq[:b0+1]

def _end_objects(kind: EndType, seg: str, for_left_right: str) -> Tuple[FragmentEnd, FragmentEnd]:
    """
    根据切口生成左右两侧片段的端对象（左片段的 3' 端；右片段的 5' 端）
    for_left_right: "5_overhang" 时右侧 5' 为 seg，左侧 3' 为 revcomp(seg)
                    "3_overhang" 时左侧 3' 为 seg，右侧 5' 为 revcomp(seg)
                    "blunt"      时二者皆 blunt
    """
    if kind == "blunt":
        fe = {"kind": "blunt", "seq": "", "length": 0}
        return fe, fe

    if for_left_right == "5_overhang":
        left_3 = {"kind": "5_overhang", "seq": _revcomp(seg), "length": len(seg)}
        right_5 = {"kind": "5_overhang", "seq": seg, "length": len(seg)}
        return left_3, right_5
    else:  # "3_overhang"
        left_3 = {"kind": "3_overhang", "seq": seg, "length": len(seg)}
        right_5 = {"kind": "3_overhang", "seq": _revcomp(seg), "length": len(seg)}
        return left_3, right_5


COMMON_ENZYMES: List[EnzymeSpec] = [
    {"name": "EcoRI",   "site": "GAATTC",   "top_cut": 1, "bottom_cut": 5},
    {"name": "BamHI",   "site": "GGATCC",   "top_cut": 1, "bottom_cut": 5},
    {"name": "BglII",   "site": "AGATCT",   "top_cut": 1, "bottom_cut": 5},
    {"name": "HindIII", "site": "AAGCTT",   "top_cut": 1, "bottom_cut": 5},
    {"name": "KpnI",    "site": "GGTACC",   "top_cut": 5, "bottom_cut": 1},
    {"name": "NcoI",    "site": "CCATGG",   "top_cut": 1, "bottom_cut": 5},
    {"name": "NdeI",    "site": "CATATG",   "top_cut": 2, "bottom_cut": 4},
    {"name": "NheI",    "site": "GCTAGC",   "top_cut": 1, "bottom_cut": 5},
    {"name": "NotI",    "site": "GCGGCCGC", "top_cut": 2, "bottom_cut": 6},
    {"name": "PacI",    "site": "TTAATTAA", "top_cut": 5, "bottom_cut": 3},
    {"name": "PstI",    "site": "CTGCAG",   "top_cut": 5, "bottom_cut": 1},
    {"name": "SacI",    "site": "GAGCTC",   "top_cut": 1, "bottom_cut": 5},
    {"name": "SalI",    "site": "GTCGAC",   "top_cut": 1, "bottom_cut": 5},
    {"name": "SmaI",    "site": "CCCGGG",   "top_cut": 3, "bottom_cut": 3},
    {"name": "XbaI",    "site": "TCTAGA",   "top_cut": 1, "bottom_cut": 5},
    {"name": "XhoI",    "site": "CTCGAG",   "top_cut": 1, "bottom_cut": 5},
    {"name": "AflII",   "site": "CTTAAG",   "top_cut": 1, "bottom_cut": 5},
]

def simulate_restriction_digestion(json_path: str,
                                   enzyme_names: List[str],
                                   output_json_path: Optional[str] = None) -> DigestResult:
    """
    Simulates the restriction digestion of a sequence with one or more enzymes simultaneously.

    Args:
        json_path (str): The path to the input JSON file containing the SequenceRecord.
        enzyme_names (List[str]): A list of enzyme names to use for the digestion.
                                  Enzyme details are looked up from an internal database.

    Returns:
        DigestResult: A single result object for the combined digestion, containing all
                      fragments produced by all enzymes cutting simultaneously.

    Raises:
        ValueError: If an enzyme name is not found or the sequence data is invalid.
    """
    seq = load_sequence_from_json(json_path)
    sequence, L, circular = _normalize_record(seq)

    enzyme_map = {e["name"]: e for e in COMMON_ENZYMES}
    enzymes: List[EnzymeSpec] = []
    for name in enzyme_names:
        if name in enzyme_map:
            enzymes.append(enzyme_map[name])
        else:
            raise ValueError(f"Enzyme '{name}' not found in COMMON_ENZYMES list.")

    if L == 0:
        return {"enzymes": enzyme_names, "cuts": [], "fragments": [], "info": []}

    all_cuts_borders: List[int] = []
    all_left3_end: Dict[int, FragmentEnd] = {}
    all_right5_end: Dict[int, FragmentEnd] = {}
    all_info_messages: List[str] = []
    MIN_FLANKING_BASES = 6

    for enz in enzymes:
        name = enz.get("name")
        site = (enz.get("site") or "").upper()
        if not name or not site:
            raise ValueError("EnzymeSpec is missing name or site")
        if re.search(r"[^ACGTN]", site):
            raise ValueError(f"Recognition site can only contain A/C/G/T/N: {name} -> {site}")

        mlen = len(site)
        tcut = int(enz.get("top_cut", 0))
        bcut = int(enz.get("bottom_cut", 0))
        if not (0 <= tcut <= mlen and 0 <= bcut <= mlen):
            raise ValueError(f"{name}: top_cut/bottom_cut must be between 0 and {mlen}")

        # --- Forward strand match ---
        for m in re.finditer(re.escape(site), sequence):
            start0 = m.start()
            if not circular and (start0 < MIN_FLANKING_BASES or (L - (start0 + mlen)) < MIN_FLANKING_BASES):
                msg = f"Site for {name} at {start0 + 1} is too close to an end for efficient cutting. Please consider additional flanking bases when desgining primers."
                if msg not in all_info_messages:
                    all_info_messages.append(msg)
            
            kind, seg = _calc_overhang(site, tcut, bcut)
            border = _cut_to_boundary(start0, tcut, L)
            all_cuts_borders.append(border)
            l3, r5 = _end_objects(kind, seg, "blunt" if kind=="blunt" else ("5_overhang" if kind=="5_overhang" else "3_overhang"))
            all_left3_end[border] = l3
            all_right5_end[border] = r5

        # --- Reverse strand match ---
        rc = _revcomp(sequence)
        for m in re.finditer(re.escape(site), rc):
            start0_rc = m.start()
            start0 = _apply_rev_orientation(start0_rc, mlen, L)
            if not circular and (start0 < MIN_FLANKING_BASES or (L - (start0 + mlen)) < MIN_FLANKING_BASES):
                msg = f"Site for {name} at {start0 + 1} (reverse) is too close to an end for efficient cutting. Please consider additional flanking bases when desgining primers."
                if msg not in all_info_messages:
                    all_info_messages.append(msg)

            tcut_p = mlen - bcut
            bcut_p = mlen - tcut
            kind, seg = _calc_overhang(site, tcut_p, bcut_p)
            border = _cut_to_boundary(start0, tcut_p, L)
            all_cuts_borders.append(border)
            l3, r5 = _end_objects(kind, seg, "blunt" if kind=="blunt" else ("5_overhang" if kind=="5_overhang" else "3_overhang"))
            all_left3_end[border] = l3
            all_right5_end[border] = r5

        # --- Circular wrap-around match ---
        if circular and mlen > 1:
            ext = sequence + sequence[:mlen-1]
            lo = L - (mlen - 1)
            for m in re.finditer(re.escape(site), ext):
                s0 = m.start()
                if s0 < lo or s0 >= L: continue
                kind, seg = _calc_overhang(site, tcut, bcut)
                border = _cut_to_boundary(s0 % L, tcut, L)
                if border not in all_left3_end:
                    all_cuts_borders.append(border)
                    all_left3_end[border], all_right5_end[border] = _end_objects(kind, seg, "blunt" if kind=="blunt" else ("5_overhang" if kind=="5_overhang" else "3_overhang"))
            
            ext_rc = rc + rc[:mlen-1]
            lo_rc = L - (mlen - 1)
            for m in re.finditer(re.escape(site), ext_rc):
                s0 = m.start()
                if s0 < lo_rc or s0 >= L: continue
                start0 = _apply_rev_orientation(s0, mlen, L)
                tcut_p = mlen - bcut
                bcut_p = mlen - tcut
                kind, seg = _calc_overhang(site, tcut_p, bcut_p)
                border = _cut_to_boundary(start0 % L, tcut_p, L)
                if border not in all_left3_end:
                    all_cuts_borders.append(border)
                    all_left3_end[border], all_right5_end[border] = _end_objects(kind, seg, "blunt" if kind=="blunt" else ("5_overhang" if kind=="5_overhang" else "3_overhang"))

    if not all_cuts_borders:
        return {"enzymes": enzyme_names, "cuts": [], "fragments": [], "info": all_info_messages}

    cuts_sorted = sorted(set(all_cuts_borders))
    cuts_1b = [(b + 1 if b < L - 1 else L) for b in cuts_sorted]
    fragments: List[Fragment] = []

    if circular:
        n = len(cuts_sorted)
        if n == 0: # No cuts, the whole plasmid is one fragment
            seq_str = sequence if _include_seq_for_len(L) else None
            frag: Fragment = {
                "id": str(uuid.uuid4()), "start": 1, "end": L, "length": L,
                "strand": 1, "overhang_5": {"kind": "blunt", "seq": "", "length": 0},
                "overhang_3": {"kind": "blunt", "seq": "", "length": 0},
                "sequence": seq_str
            }
            fragments.append(frag)
        else:
            for i in range(n):
                b_left = cuts_sorted[i]
                b_right = cuts_sorted[(i + 1) % n]
                a0 = (b_left + 1) % L
                b0 = b_right % L
                length = (b0 - a0 + 1) if a0 <= b0 else (L - a0) + (b0 + 1)
                seq_str = _frag_seq(sequence, a0, b0, True) if _include_seq_for_len(length) else None
                frag: Fragment = {
                    "id": str(uuid.uuid4()), "start": a0 + 1, "end": b0 + 1, "length": int(length),
                    "strand": 1, "overhang_5": all_right5_end[b_left], "overhang_3": all_left3_end[b_right],
                    "sequence": seq_str
                }
                fragments.append(frag)
    else: # Linear
        cuts_lin = [b for b in cuts_sorted if b <= L - 2]
        points = [-1] + cuts_lin + [L - 1]
        for i in range(len(points) - 1):
            left_b, right_b = points[i], points[i + 1]
            a0, b0 = left_b + 1, right_b
            length = (b0 - a0 + 1) if b0 >= a0 else 0
            if length <= 0: continue
            seq_str = _frag_seq(sequence, a0, b0, False) if _include_seq_for_len(length) else None
            overhang_5 = {"kind": "blunt", "seq": "", "length": 0} if left_b == -1 else all_right5_end[left_b]
            overhang_3 = {"kind": "blunt", "seq": "", "length": 0} if right_b == L - 1 and right_b not in all_left3_end else all_left3_end.get(right_b, {"kind": "blunt", "seq": "", "length": 0})
            frag: Fragment = {
                "id": str(uuid.uuid4()), "start": a0 + 1, "end": b0 + 1, "length": int(length),
                "strand": 1, "overhang_5": overhang_5, "overhang_3": overhang_3,
                "sequence": seq_str
            }
            fragments.append(frag)

    # If output_json_path is provided, save the fragments to a JSON file
    if output_json_path:
        # Create a dictionary with the fragments and other relevant information
        output_data = {
            "enzymes": enzyme_names,
            "cuts": [int(x) for x in cuts_1b],
            "fragments": fragments,
            "info": all_info_messages
        }
        write_record_to_json(output_data, output_json_path)
    
    return {
        "enzymes": enzyme_names,
        "cuts": [int(x) for x in cuts_1b],
        "fragments": fragments,
        "info": all_info_messages
    }
