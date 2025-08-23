# get_cds.py
from Bio import Entrez
import os, time, logging, re, json, sys
from pathlib import Path
from datetime import datetime
from typing import Dict, Any, Optional

# Add the parent directory to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))
from common_utils.sequence import SequenceRecord
from common_utils.file_operations import write_record_to_json

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

def get_cds_by_gene(gene_name: str, organism: str) -> Dict[str, Any]:
    """
    从NCBI获取指定基因的CDS序列并保存为SequenceRecord JSON格式。

    重要提示：在调用此工具前，如果用户没有明确指定物种，请先询问用户需要哪个物种的基因。

    Args:
        gene_name (str): 基因名称（如：RTCB, TP53, GAPDH）
        organism (str): 物种的学名，必须提供。常用物种包括：
            - "Homo sapiens" (人类)
            - "Mus musculus" (小鼠)
            - "Rattus norvegicus" (大鼠)  
            - "Escherichia coli" (大肠杆菌)
            - "Saccharomyces cerevisiae" (酿酒酵母)
            - "Drosophila melanogaster" (果蝇)
            注意：必须使用标准的拉丁学名格式

    Returns:
        Dict[str, Any]: 包含sequence_record_path和info字段的字典
        
    注意：
    - 文件将自动保存到 data/temp_cds/ 目录下
    - 如果用户只提到基因名而未指定物种，请在调用此工具前先询问用户需要哪个物种的基因序列。
    """
    try:
        # 验证参数
        if not gene_name or not gene_name.strip():
            return {
                "sequence_record_path": None,
                "info": "错误：基因名称不能为空。请提供有效的基因名称。"
            }
        
        if not organism or not organism.strip():
            return {
                "sequence_record_path": None,
                "info": "错误：organism参数不能为空。请使用标准拉丁学名，如'Homo sapiens'。如果用户未指定物种，建议询问用户需要哪个物种的基因序列。"
            }

        logging.info(f"Searching gene: {gene_name} ({organism})")

        # 1) 基因名 → Gene ID
        term = f"{gene_name}[gene] AND {organism}[orgn]"
        with Entrez.esearch(db="gene", term=term, retmax=1) as h:
            rec = Entrez.read(h)
        if not rec.get("IdList"):
            logging.error(f"No Gene ID found for {gene_name}")
            return {
                "sequence_record_path": None,
                "info": f"未找到基因：{gene_name}（{organism}）。请检查基因名称和物种名称是否正确。常见物种格式：'Homo sapiens', 'Mus musculus'等。"
            }
        gene_id = rec["IdList"][0]
        logging.info(f"Found Gene ID: {gene_id}")
        time.sleep(0.34)

        # 2) Gene ID → RefSeq mRNA（取第一个）
        with Entrez.elink(dbfrom="gene", db="nuccore", id=gene_id, linkname="gene_nuccore_refseqrna") as h:
            links = Entrez.read(h)
        if not links or not links[0].get("LinkSetDb"):
            logging.error("No RefSeq mRNA found")
            return {
                "sequence_record_path": None,
                "info": f"未找到基因 {gene_name} 的 RefSeq mRNA 记录。该基因可能没有注释的转录本，或者需要检查基因名称是否正确。"
            }
        nuccore_id = links[0]["LinkSetDb"][0]["Link"][0]["Id"]
        logging.info(f"Found nuccore ID: {nuccore_id}")
        time.sleep(0.34)

        # 3) 取 CDS FASTA
        with Entrez.efetch(db="nuccore", id=nuccore_id, rettype="fasta_cds_na", retmode="text") as h:
            data = h.read().strip()
        if not data:
            logging.error("No CDS data found for selected transcript")
            return {
                "sequence_record_path": None,
                "info": f"该转录本无 CDS 数据。基因 {gene_name} 可能是非编码基因，或者该转录本缺少CDS注释。"
            }

        # 解析FASTA数据
        lines = data.splitlines()
        header = lines[0]
        sequence_lines = [line.strip() for line in lines[1:] if line.strip() and not line.startswith(">")]
        sequence = "".join(sequence_lines).upper()
        
        # 从FASTA头部提取信息
        accession = header.split(" ")[0][1:]  # 去掉前导的 '>'
        description = " ".join(header.split(" ")[1:]) if len(header.split(" ")) > 1 else ""
        
        # 创建SequenceRecord
        organism_safe = organism.replace(" ", "_")
        record_id = f"{gene_name}_{organism_safe}_{accession}"
        
        sequence_record: SequenceRecord = {
            "id": record_id,
            "name": f"{gene_name} CDS from {organism}",
            "sequence": sequence,
            "length": len(sequence),
            "circular": False,  # CDS通常是线性的
            "features": [],  # 初始为空，可后续用find_features工具添加
            "metadata": {
                "source": "NCBI",
                "accession": accession,
                "gene_name": gene_name,
                "organism": organism,
                "description": description,
                "retrieval_date": datetime.now().isoformat(),
                "database": "nuccore",
                "nuccore_id": str(nuccore_id)
            }
        }

        # 固定输出路径到 data/temp_cds/ 目录
        output_dir = Path(__file__).resolve().parent.parent / "data" / "temp_cds"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 生成安全文件名
        raw_file_name = f"{gene_name}_{organism_safe}_{accession}.json"
        file_name = safe_filename(raw_file_name)
        json_path = str(output_dir / file_name)

        # 保存为JSON
        write_record_to_json(sequence_record, json_path)
        logging.info(f"CDS sequence saved as JSON to: {json_path}")

        return {
            "sequence_record_path": json_path,
            "info": f"成功获取基因 {gene_name} ({organism}) 的CDS序列，长度 {len(sequence)} bp，已保存为SequenceRecord JSON格式。可以使用 get_sequence_info 工具查看详细信息，或使用 find_features 工具查找特征位点。"
        }

    except Exception as e:
        logging.error(f"Error in get_cds_by_gene: {str(e)}")
        return {
            "sequence_record_path": None,
            "info": f"获取CDS序列时发生错误：{str(e)}。请检查网络连接、基因名称和物种名称是否正确。"
        }


# if __name__ == "__main__":
#     # 直接调试用
#     result = get_cds_by_gene("RTCB", "Homo sapiens")
#     print(result)
