# snapgene_to_sequence_record.py
# 极简：把 snapgene_reader.snapgene_file_to_dict(...) 的 dict 转为 SequenceRecord 并打印 JSON

import sys, os, json, uuid
from typing import Dict, Any
from tools_pool.get_sequence_info import *

def _assign_id(prefix: str = "seq") -> str:
    return f"{prefix}_{uuid.uuid4().hex[:8]}"

def _strand_to_int(s):
    if s in (1, "+", "+1", "plus"): return 1
    if s in (-1, "-", "-1", "minus"): return -1
    return 1  # '.', None 等统一按正链

def _first_segment_range(segments):
    # 解析类似 "@range": "882-1026"
    try:
        if not segments: return None
        rng = segments[0].get("@range") or segments[0].get("range")
        a, b = str(rng).replace(",", "").split("-")
        return int(a), int(b)
    except Exception:
        return None

def convert_snapgene_dict_to_sequence_record(d: dict, file_name: str = "") -> dict:
    seq = (d.get("seq") or d.get("sequence") or "").upper()
    if not seq:
        raise ValueError("输入字典中未找到序列（'seq'）")

    dna_meta = d.get("dna") or {}
    topology = str(dna_meta.get("topology", "")).lower()
    circular = topology == "circular"

    feats_in = d.get("features") or []
    feats_out = []
    for f in feats_in:
        start = f.get("start")
        end = f.get("end")
        if start is None or end is None:
            seg = _first_segment_range(f.get("segments") or [])
            if seg:
                start, end = seg
        # SnapGene 内多数是 1-based 区间；若遇到 0-based 可在此加统一+1校正
        start = int(start) if start is not None else 1
        end = int(end) if end is not None else start

        strand = _strand_to_int(f.get("strand"))
        ftype = f.get("type") or "misc_feature"

        # 组装 qualifiers：把除了核心坐标/类型外的字段都塞进来
        qualifiers_raw = {}
        for k, v in f.items():
            if k in {"type", "start", "end", "strand", "segments", "color", "textcolor"}: 
                continue
            qualifiers_raw[k] = v
        # 常见别名放进 qualifiers
        if f.get("name") and "name" not in qualifiers_raw:
            qualifiers_raw["name"] = f["name"]
        if isinstance(f.get("qualifiers"), dict):
            qualifiers_raw.update(f["qualifiers"])

        # 只保留 label
        qualifiers = {}
        if "label" in qualifiers_raw:
            qualifiers["label"] = qualifiers_raw["label"]
        elif "name" in qualifiers_raw:
            qualifiers["label"] = qualifiers_raw["name"]


        feats_out.append({
            "type": ftype,
            "start": start,            # 1-based inclusive
            "end": end,                # 1-based inclusive
            "strand": strand,          # 1 或 -1
            "qualifiers": qualifiers
        })

    # metadata
    metadata = {
        "parser": "snapgene_reader",
        "topology": topology,
        "source_file": os.path.abspath(file_name) if file_name else "",
    }

    record = {
        "id": _assign_id("seq"),
        "name": d.get("name") or os.path.basename(file_name) or "unnamed",
        "sequence": seq,
        "length": len(seq),
        "circular": circular,
        "features": feats_out,
        "metadata": metadata,
    }
    return record

def convert_fasta_record_to_sequence_record(fasta_record: Any, file_name: str = "") -> Dict[str, Any]:
    """
    将 Biopython SeqRecord 对象转换为 SequenceRecord 字典格式。
    """
    seq = str(fasta_record.seq).upper()
    if not seq:
        raise ValueError("FASTA 记录中未找到序列")

    # FASTA 文件通常不包含拓扑信息，默认为线性
    circular = False

    # FASTA 文件通常不包含 features，可以根据需要从 description 或 id 中解析
    feats_out = []
    # 示例：如果需要从 description 中解析 feature，可以在这里添加逻辑
    # if fasta_record.description:
    #     # 尝试从 description 中提取信息创建 feature
    #     pass

    metadata = {
        "parser": "biopython_fasta",
        "source_file": os.path.abspath(file_name) if file_name else "",
        "description": fasta_record.description,
        "annotations": dict(fasta_record.annotations), # 转换为字典以确保可序列化
    }

    record = {
        "id": _assign_id("seq"),
        "name": fasta_record.id or os.path.basename(file_name) or "unnamed",
        "sequence": seq,
        "length": len(seq),
        "circular": circular,
        "features": feats_out,
        "metadata": metadata,
    }
    return record

def load_sequence(input_path: str) -> dict:
    """
    读取 SnapGene 文件（.dna）或 FASTA 文件（.fasta, .fa, .fna）
    并转换为 SequenceRecord JSON，并保存到 data/temp 目录下。

    Args:
        input_path (str): 输入的文件路径。

    Returns:
        dict: 输出的 JSON 文件路径。以及序列基本信息。
    """
    file_extension = os.path.splitext(input_path)[1].lower()
    rec = None

    if file_extension == ".dna":
        try:
            from snapgene_reader import snapgene_file_to_dict
        except ImportError:
            raise ImportError("处理 .dna 文件需要安装 snapgene-reader：pip install snapgene-reader")
        d = snapgene_file_to_dict(input_path)
        rec = convert_snapgene_dict_to_sequence_record(d, file_name=input_path)
    elif file_extension in (".fasta", ".fa", ".fna"):
        try:
            from Bio import SeqIO
        except ImportError:
            raise ImportError("处理 FASTA 文件需要安装 Biopython：pip install biopython")
        
        # Biopython SeqIO.parse 返回一个迭代器，通常 FASTA 文件只包含一个序列
        # 如果文件包含多个序列，这里只处理第一个
        fasta_records = list(SeqIO.parse(input_path, "fasta"))
        if not fasta_records:
            raise ValueError(f"FASTA 文件 '{input_path}' 中未找到序列。")
        
        rec = convert_fasta_record_to_sequence_record(fasta_records[0], file_name=input_path)
    else:
        raise ValueError(f"不支持的文件类型: {file_extension}。目前只支持 .dna 和 .fasta/.fa/.fna 文件。")

    if rec is None:
        raise RuntimeError("序列转换失败，未生成 SequenceRecord。")

    output_dir = os.path.join("data", "temp")
    os.makedirs(output_dir, exist_ok=True)

    base_name = os.path.basename(input_path)
    file_name_without_ext = os.path.splitext(base_name)[0]
    output_path = os.path.join(output_dir, f"{file_name_without_ext}.json")

    with open(output_path, "w", encoding="utf-8") as fh:
        json.dump(rec, fh, ensure_ascii=False, indent=2)
    
    info = get_sequence_info(output_path)
    return {
        "message": f"文件转换成功，已保存至: {output_path}",
        "sequence_info": info
    }


if __name__ == "__main__":
    # 示例：处理 SnapGene 文件
    snapgene_file = "D:\github_repo\PPversion2\data\pcDNA3.1(-).dna"
    print(f"处理 SnapGene 文件: {snapgene_file}")
    load_sequence(snapgene_file)

    print("-" * 30)

    # 示例：处理 FASTA 文件
    fasta_file = "D:\github_repo\PPversion2\data\RTCB_Homo_sapiens_lcl_NM_014306.5_cds_NP_055121.1_1.fasta"
    print(f"处理 FASTA 文件: {fasta_file}")
    load_sequence(fasta_file)
