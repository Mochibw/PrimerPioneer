from typing import List, Dict, Optional, Literal, TypedDict, Tuple, Union
from typing_extensions import TypedDict, Literal, Union
import uuid
import json

Strand = Literal[1, -1]
EndType = Literal["blunt", "5_overhang", "3_overhang"]

class SequenceRecord(TypedDict, total=False):
    id: str                    # unique ID for the record
    name: str                  # human-readable name
    sequence: str              # uppercase A/C/G/T (IUPAC optional: RYKMSWBDHVN)
    length: int                # len(sequence)
    circular: bool             # plasmid/topology flag
    features: List["Feature"]  # optional feature annotations
    metadata: Dict[str, str]   # free-form key-value pairs

class Feature(TypedDict, total=False):
    type: str                  # e.g., 'CDS', 'promoter', 'ori', 'MCS', 'polyA', 'tag'
    start: int                 # 1-based inclusive
    end: int                   # 1-based inclusive
    strand: Strand             # 1 or -1
    qualifiers: Dict[str, Union[str, float, int]]

class FragmentEnd(TypedDict, total=False):
    kind: EndType              # 'blunt' | '5_overhang' | '3_overhang'
    seq: str                   # overhang sequence ('' if blunt)
    length: int                # overhang length (0 if blunt)

class Fragment(TypedDict, total=False):
    id: str
    start: int
    end: int
    length: int
    strand: Strand
    overhang_5: FragmentEnd
    overhang_3: FragmentEnd
    sequence: Optional[str]    # optional for large constructs to save payload

class Primer(TypedDict, total=False):
    name: str
    sequence: str
    tm: float
    gc: float
    start: Optional[int]       # binding start (1-based, on + strand coordinates)
    end: Optional[int]         # binding end (1-based)
    direction: Literal["F", "R"]
    length: int

class AssemblyStep(TypedDict, total=False):
    action: str                # e.g., 'PCR', 'Digest', 'Phosphorylation', 'Ligation', 'Gibson'
    inputs: List[str]          # IDs of fragments/sequences consumed
    outputs: List[str]         # IDs produced
    params: Dict[str, Union[str, int, float, bool]]

class AssemblyPlan(TypedDict, total=False):
    method: str                # 'restriction_ligation' | 'gibson' | 'isothermal' | 'TA' | 'custom'
    steps: List[AssemblyStep]
    rationale: str
    expected_product: SequenceRecord


def get_sequence(json_path: str,
                 head: Optional[int] = 200,
                 tail: Optional[int] = None,
                 slice_1based: Optional[Tuple[int, int]] = None,
                 full_length: bool = False) -> Tuple[str, int]:
    """
    从任何包含 "sequence" 字段的 JSON 文件中读取序列，并返回序列及其总长度。

    默认情况下，为避免返回内容过长，函数仅返回序列的前 200 个碱基。
    可以通过设置 `full_length=True` 来获取完整序列（极不推荐LLM使用）。

    Args:
        json_path (str): 输入 JSON 文件的路径。
        head (Optional[int]): 返回序列的前 N 个碱基。默认为 200。
        tail (Optional[int]): 返回序列的后 N 个碱基。如果设置，则忽略 head。
        slice_1based (Optional[Tuple[int, int]]): 返回 1-based 闭区间的序列片段 [start, end]。如果设置，则忽略 head 和 tail。
        full_length (bool): 如果为 True，则返回完整序列，忽略所有分块参数。

    Returns:
        Tuple[str, int]: 一个元组，包含序列字符串（或其片段）和原始序列的总长度。
                         如果找不到 "sequence" 字段，则返回 ("", 0)。
    """
    with open(json_path, "r", encoding="utf-8") as f:
        record: Dict = json.load(f)
    
    sequence = record.get("sequence")
    if not isinstance(sequence, str):
        return "", 0

    total_len = len(sequence)

    if full_length:
        return sequence, total_len

    if slice_1based:
        start, end = slice_1based
        return sequence[max(0, start - 1):end], total_len
    
    if tail is not None:
        return sequence[-tail:], total_len
    
    if head is not None:
        return sequence[:head], total_len

    # Fallback if all optional args are None (though head has a default)
    return sequence, total_len

