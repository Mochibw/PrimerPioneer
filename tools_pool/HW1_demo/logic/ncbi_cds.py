# ncbi_cds.py
from Bio import Entrez
import os, time, logging, re
from pathlib import Path

# 配置日志输出到 stderr
logging.basicConfig(level=logging.INFO, format='[%(levelname)s] %(message)s')

# 必须填邮箱（NCBI 要求）
Entrez.email = os.getenv("NCBI_EMAIL", "fym22@mails.tsinghua.edu.cn")
Entrez.api_key = os.getenv("NCBI_API_KEY", None)

# ------- Windows/跨平台安全文件名处理 -------
_INVALID_CHARS = r'[<>:"/\\|?*\x00-\x1F]'
_RESERVED_WIN = {
    "CON","PRN","AUX","NUL",
    *(f"COM{i}" for i in range(1,10)),
    *(f"LPT{i}" for i in range(1,10)),
}
def safe_filename(name: str, replacement: str = "_") -> str:
    # 替换非法字符
    name = re.sub(_INVALID_CHARS, replacement, name)
    # 去掉结尾的点/空格（Windows 不允许）
    name = name.rstrip(" .")
    # 避免保留名
    stem, dot, ext = name.partition(".")
    if stem.upper() in _RESERVED_WIN:
        stem += "_"
    return stem + (dot + ext if dot else "")

def write_text(path: Path, text: str, encoding: str = "utf-8") -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text, encoding=encoding)

def fasta_seq_length(fasta_text: str) -> int:
    lines = fasta_text.splitlines()
    seq = "".join(l.strip() for l in lines[1:] if l and not l.startswith(">"))
    return len(seq)

def get_cds_by_gene_simple(gene_name: str, organism: str = "Homo sapiens") -> str:
    logging.info(f"Searching gene: {gene_name} ({organism})")

    # 1) 基因名 → Gene ID
    term = f"{gene_name}[gene] AND {organism}[orgn]"
    with Entrez.esearch(db="gene", term=term, retmax=1) as h:
        rec = Entrez.read(h)
    if not rec.get("IdList"):
        logging.error(f"No Gene ID found for {gene_name}")
        raise ValueError(f"未找到基因：{gene_name}（{organism}）")
    gene_id = rec["IdList"][0]
    logging.info(f"Found Gene ID: {gene_id}")
    time.sleep(0.34)

    # 2) Gene ID → RefSeq mRNA（取第一个）
    with Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore_refseqrna") as h:
        links = Entrez.read(h)
    if not links or not links[0].get("LinkSetDb"):
        logging.error("No RefSeq mRNA found")
        raise ValueError("未找到 RefSeq mRNA")
    nuccore_id = links[0]["LinkSetDb"][0]["Link"][0]["Id"]
    logging.info(f"Found nuccore ID: {nuccore_id}")
    time.sleep(0.34)

    # 3) 取 CDS FASTA
    with Entrez.efetch(db="nuccore", id=nuccore_id, rettype="fasta_cds_na", retmode="text") as h:
        data = h.read().strip()
    if not data:
        logging.error("No CDS data found for selected transcript")
        raise ValueError("该转录本无 CDS 数据")

    # 更准确的长度（合并所有序列行）
    length_bp = fasta_seq_length(data)
    logging.info(f"Retrieved CDS FASTA length (merged): {length_bp} bp")

    # 创建用于存放中间文件的目录
    # 修正路径，使其指向项目根目录下的 data/temp_cds
    output_dir = Path(__file__).resolve().parent.parent / "data" / "temp_cds"
    output_dir.mkdir(parents=True, exist_ok=True)

    # 从 FASTA 第一行提取 accession；常见为 >lcl|NM_014306.5 ...
    header = data.splitlines()[0]
    accession = header.split(" ")[0][1:]  # 去掉前导的 '>'
    # 生成安全文件名（会把 '|' 等非法字符替换为下划线）
    raw_file_name = f"{gene_name}_{organism.replace(' ', '_')}_{accession}.fasta"
    file_name = safe_filename(raw_file_name)
    file_path = output_dir / file_name

    write_text(file_path, data, encoding="utf-8")
    logging.info(f"CDS sequence saved to: {file_path}")

    # 返回文件路径（字符串）
    return str(file_path)


# if __name__ == "__main__":
#     # 直接调试用
#     gene = "RTCB"
#     org = "Homo sapiens"
#     print(get_cds_by_gene_simple(gene, org))
