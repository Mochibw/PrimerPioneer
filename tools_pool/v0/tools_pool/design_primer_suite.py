import sys
import os
from dataclasses import dataclass
from typing import Optional, Dict, List, Literal, Any

# 允许从上级目录导入
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from common_utils.sequence import get_sequence, Primer
from common_utils.file_operations import load_sequence
from common_utils.sequence_tools import calculate_gc_content, calculate_tm, reverse_complement
from tools_pool.select_cloning_method import get_method_recommendation
from langchain.tools import tool

# =========================
# 变化摘要（读这几行就够了）
# - 外部接口仅保留: design_primer_suite(target, task="amplify", intent=None)
# - 不再接受 tm/长度/3' clamp 等可调参数；内部用通用默认策略
# - 文档与注释大幅精简；只保留任务必要说明
# - 3' GC clamp 作为“偏好”而非强制
# =========================

# ===== 基本工具 =====
def revcomp(seq: str) -> str:
    return reverse_complement(seq)

def gc_content(seq: str) -> float:
    return calculate_gc_content(seq)

def tm_wallace(seq: str) -> float:
    return calculate_tm(seq, method="wallace")

# ===== 轻量“锚点+载荷”描述 =====
@dataclass
class Anchor:
    on: str
    pos: int           # 0-based; 5' 末端在模板的位置
    strand: int        # +1 (F) / -1 (R)
    overhang: str = "" # 5' 载荷（可空）

@dataclass
class PrimerBlueprint:
    forward: Anchor
    reverse: Anchor
    notes: str = ""

# ===== 任务语义 → 蓝图 =====
def intent_to_blueprint(task: str,
                        target_seq: str,
                        intent: Dict) -> PrimerBlueprint:
    """Translates a task-specific intent into a standardized primer blueprint.

    This function interprets the user's high-level goal (e.g., "amplify a region")
    and converts it into a pair of `Anchor` objects that define the precise
    locations and properties for the forward and reverse primers.

    All coordinates in the `intent` dictionary are 1-based and inclusive,
    and are converted to 0-based internally.

    Supported `intent` structures:
    - for task="amplify":
        intent={"amplicon":{"start":int,"end":int}, "fwd_overhang":"", "rev_overhang":""}
    - for task="mutagenesis":
        - Point mutation: {"edits":[{"type":"point","pos":int,"to":"A|C|G|T"}]}
        - Insertion: {"edits":[{"type":"insertion","pos":int,"insert":"ATGC...","embed":"core|overhang"}]}
        - Deletion: {"edits":[{"type":"deletion","start":int,"end":int}]}

    Args:
        task (str): The primer design task, e.g., "amplify", "mutagenesis".
        target_seq (str): The template DNA sequence.
        intent (Dict): A dictionary containing task-specific parameters.

    Returns:
        PrimerBlueprint: A data object containing the forward and reverse anchors.

    Raises:
        ValueError: If the task is unknown, or if the `intent` dictionary
                    is missing required keys or contains invalid values.
    """
    if intent is None:
        intent = {}

    def need(keys):
        missing = [k for k in keys if k not in intent]
        if missing:
            raise ValueError(f"缺少必需字段: {missing}")

    if task == "amplify":
        need(["amplicon"])
        amp = intent["amplicon"]
        start = amp.get("start"); end = amp.get("end")
        if start is None or end is None or not (1 <= start <= end <= len(target_seq)):
            raise ValueError("amplicon 需要合法的 [start,end] (1-based, 闭区间)")
        F = Anchor(on="A", pos=start - 1,      strand=+1, overhang=intent.get("fwd_overhang", "") or "")
        R = Anchor(on="A", pos=end,            strand=-1, overhang=intent.get("rev_overhang", "") or "")
        return PrimerBlueprint(F, R, notes=f"Amplify A[{start}:{end}]")

    elif task == "mutagenesis":
        need(["edits"])
        edits = intent["edits"]
        if not isinstance(edits, list) or len(edits) != 1:
            raise ValueError("最小实现仅支持单一编辑（edits 长度=1）")
        e = edits[0]
        etype = e.get("type", "point")

        if etype == "point":
            pos = e.get("pos"); to = e.get("to")
            if pos is None or to is None or not (1 <= pos <= len(target_seq)):
                raise ValueError("point 需要合法的 pos 与 to")
            pos0 = pos - 1
            F = Anchor(on="A", pos=pos0,     strand=+1)
            R = Anchor(on="A", pos=pos0 + 1, strand=-1)
            return PrimerBlueprint(F, R, notes=f"Point mutation A[{pos}] -> {to}")

        elif etype == "insertion":
            pos = e.get("pos"); ins = e.get("insert", "")
            if pos is None or not ins or not (1 <= pos <= len(target_seq) + 1):
                raise ValueError("insertion 需要合法的 pos 与 insert")
            pos0 = pos - 1
            embed = e.get("embed", "core")
            if embed == "core":
                F = Anchor(on="A", pos=pos0, strand=+1)
                R = Anchor(on="A", pos=pos0, strand=-1)
                note = "Insertion in cores"
            else:
                F = Anchor(on="A", pos=pos0, strand=+1, overhang=ins)
                R = Anchor(on="A", pos=pos0, strand=-1, overhang=revcomp(ins))
                note = "Insertion via overhangs"
            return PrimerBlueprint(F, R, notes=f"{note} at A[{pos}]")

        elif etype == "deletion":
            ds = e.get("start"); de = e.get("end")
            if ds is None or de is None or not (1 <= ds <= de <= len(target_seq)):
                raise ValueError("deletion 需要合法的 [start,end]")
            F = Anchor(on="A", pos=ds - 1,      strand=+1)
            R = Anchor(on="A", pos=de,          strand=-1)
            return PrimerBlueprint(F, R, notes=f"Deletion A[{ds}:{de}]")

        else:
            raise ValueError(f"不支持的 mutagenesis 类型: {etype}")

    else:
        raise ValueError(f"未知 task: {task}")

# ===== 蓝图 → 引物（核心 + 5'载荷） =====
# 内部默认策略（不对外暴露）：
TM_MIN, TM_MAX = 58.0, 64.0
LEN_MIN, LEN_MAX = 18, 28
PREFER_GC_CLAMP = True  # 偏好 3' 为 G/C，但不强制

def _prefer_gc_3p(seq: str) -> bool:
    return bool(seq) and (seq[-1].upper() in "GC")

def realize_core(anchor: Anchor, seqs: Dict[str, str]) -> str:
    """
    从 anchor.pos 开始，选择一个在常规范围内“更像好引物”的核心序列：
    - Tm 接近 61℃ 更好；长度 18~28；3' 末端为 G/C 给予加分（不强制）
    - 命中范围优先；若完全命不中，返回 Tm 最接近的候选
    """
    name = anchor.on
    tmpl = seqs[name]
    pos = anchor.pos

    best = ""
    best_score = float("inf")
    center_tm = 0.5 * (TM_MIN + TM_MAX)

    def score(core: str) -> float:
        t = tm_wallace(core)
        gap = abs(t - center_tm)
        bonus = 0.5 if (PREFER_GC_CLAMP and _prefer_gc_3p(core)) else 0.0
        # 越小越好；带一点 3' GC 的优待
        return gap - bonus

    # 先尝试在范围内直接返回；否则保留最优
    inrange_best = ""
    inrange_score = float("inf")

    for L in range(LEN_MIN, LEN_MAX + 1):
        if anchor.strand == +1:
            start, end = pos, pos + L
            if end > len(tmpl): 
                continue
            core = tmpl[start:end]
        else:
            start, end = pos - L, pos
            if start < 0:
                continue
            core = revcomp(tmpl[start:end])

        t = tm_wallace(core)
        sc = score(core)

        if TM_MIN <= t <= TM_MAX:
            if sc < inrange_score:
                inrange_score = sc
                inrange_best = core

        if sc < best_score:
            best_score = sc
            best = core

    if inrange_best:
        return inrange_best
    if best:
        return best
    raise ValueError("在模板边界内未找到可用核心；可调整任务坐标后重试")

def make_primer_record(name: str, overhang: str, core: str, anchor: Anchor) -> Primer:
    seq = (overhang or "") + core
    core_len = len(core)
    if anchor.strand == 1:  # F
        start_pos = anchor.pos + 1
        end_pos = anchor.pos + core_len
    else:                   # R
        start_pos = anchor.pos - core_len + 1
        end_pos = anchor.pos
    return {
        "name": name,
        "sequence": seq,
        "tm": tm_wallace(core),
        "gc": gc_content(core),
        "start": start_pos,   # 1-based, core 区间
        "end": end_pos,       # 1-based, core 区间
        "direction": "F" if anchor.strand == 1 else "R",
        "length": len(seq)
    }

# ===== 主函数（极简接口） =====
@tool()
def design_primer_suite(target: str,
                        task: Literal["amplify", "mutagenesis"] = "amplify",
                        intent: Optional[Dict] = None) -> Dict[str, Any]:
    """
    设计一对引物。

    主要入口：输入目标序列 + 任务类型 + 参数字典，返回包含前/后引物的设计结果。

    Args:
        target (str):
            DNA 模板序列，可以是：
            - 直接的 ATGC 字符串
            - 序列文件路径 (.dna, .fasta, .fa, .fna)
            - 已由 其他MCP tools 处理产生的 JSON 路径（推荐）
        task (str):
            引物设计任务类型，必须是 "amplify" 或 "mutagenesis"
        intent (dict):
            任务参数，格式依赖于 task。

            当 task="amplify":
                {
                "amplicon": {"start": int, "end": int},   # 1-based 闭区间
                "fwd_overhang": Optional[str],            # 可选，前引物 5' 载荷
                "rev_overhang": Optional[str]             # 可选，后引物 5' 载荷
                }

            当 task="mutagenesis":
                {
                "edits": [ { ... } ]   # 列表长度必须为 1
                }
                其中 edit 对象可以是：
                - 点突变:
                {"type": "point", "pos": int, "to": "A|C|G|T"}
                - 插入:
                {"type": "insertion", "pos": int, "insert": "ATGC...", "embed": "core|overhang"}
                - 缺失:
                {"type": "deletion", "start": int, "end": int}

    Returns:
        dict:
            成功:
            {
            "status": "success",
            "message": "Primer design successful.",
            "data": {
                "pairs": [
                {
                    "forward": {
                    "name": "FWD",
                    "sequence": str,   # 完整引物序列 (overhang + core)
                    "tm": float,       # core Tm
                    "gc": float,       # core GC 含量
                    "start": int,      # core 起始位点 (1-based)
                    "end": int,        # core 终止位点 (1-based)
                    "direction": "F",
                    "length": int      # 总长度
                    },
                    "reverse": { ... },
                    "notes": str
                }
                ],
                "singles": []
            }
            }

            失败:
            {
            "status": "error",
            "message": "错误描述"
            }
    """

    method_rec = get_method_recommendation(target)
    if method_rec["recommended_method"] == "synthesis":
        return {
            "status": "warning",
            "message": f"此序列更适合化学合成: {method_rec['reasoning']}",
            "recommendation": "考虑使用simulate_dna_synthesis工具",
            "critical_issues": method_rec["critical_issues"]
        }

    try:
        # 读取模板序列；允许传入文件路径或裸序列
        def seq_of(x):
            if x is None:
                return None
            if isinstance(x, str):
                # Check if it's a file path that needs conversion
                if os.path.exists(x) and x.endswith((".dna", ".fasta", ".fa", ".fna")):
                    try:
                        # Convert to JSON and get the new path
                        result = load_sequence(x)
                        json_path = result["sequence_info"]["json_path"]
                        seq, _ = get_sequence(json_path, full_length=True)
                        return seq
                    except Exception as e:
                        raise ValueError(f"Failed to process sequence file {x}: {e}")
                # Try to read as a JSON path directly
                try:
                    seq, _ = get_sequence(x, full_length=True)
                    return seq
                except Exception:
                    return x  # 不是路径则当作裸序列
            return getattr(x, "sequence", None) or getattr(x, "seq", None)

        A_seq = seq_of(target)
        if not A_seq:
            raise ValueError("target 缺少序列（应为文件路径或有效序列字符串）")

        seqs = {"A": A_seq}

        # 1) 任务 → 蓝图
        bp = intent_to_blueprint(task, A_seq, intent or {})

        # 2) 蓝图 → 核心
        f_core = realize_core(bp.forward, seqs)
        r_core = realize_core(bp.reverse, seqs)

        # 3) 拼接 5' 载荷 + 核心
        f_rec = make_primer_record("FWD", bp.forward.overhang, f_core, bp.forward)
        r_rec = make_primer_record("REV", bp.reverse.overhang, r_core, bp.reverse)

        return {
            "status": "success",
            "message": "Primer design successful.",
            "data": {"pairs": [{"forward": f_rec, "reverse": r_rec, "notes": bp.notes}], "singles": []}
        }
    except ValueError as e:
        return {"status": "error", "message": str(e)}
    except Exception as e:
        return {"status": "error", "message": f"An unexpected error occurred: {str(e)}"}

if __name__ == "__main__":
    # 示例：扩增 1..500 区段，带简单 overhang
    res = design_primer_suite(
        "D:/github_repo/PPversion2/data/temp/pcDNA3.1(-).json",
        task="amplify",
        intent={
            "amplicon": {"start": 1, "end": 500},
            "fwd_overhang": "ATGC",
            "rev_overhang": "GGCC",
        },
    )
    print(res)
