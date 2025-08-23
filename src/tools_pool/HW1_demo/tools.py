from logic.pick_restric_enzym_pairs import pick_enzyme_pairs_from_dna
from logic.ncbi_cds import get_cds_by_gene_simple
from logic.primer_design import design_primers_logic
from logic.fasta_utils import read_fasta

from langchain.tools import tool

# =========================
# 工具定义（各 worker / evaluator 拥有不同工具池）
# =========================
@tool()
def read_fasta_file(path: str) -> str:
    """
    读取FASTA文件并返回其序列内容。

    Args:
        path: FASTA文件的路径。

    Returns:
        FASTA文件中的序列内容。
    """
    return read_fasta(path)


@tool()
def get_cds_sequence(gene_name: str, organism: str = "Homo sapiens") -> str:
    """
    从 NCBI 获取基因的 CDS 序列，并保存到临时文件，返回文件路径。
    Args:
        gene_name: 基因名称。
        organism: 物种名称，默认 "Homo sapiens"。
    Returns:
        保存CDS序列的FASTA文件路径。
    """
    return get_cds_by_gene_simple(gene_name, organism)

@tool()
def select_restriction_sites(vector_dna_path: str, insert_sequence_path: str) -> dict:
    """
    根据提供的载体和插入序列的文件路径，选择合适的双酶切位点。

    Args:
        vector_dna_path: 载体DNA文件的路径。
        insert_sequence_path: 插入片段FASTA文件的路径。

    Returns:
        包含推荐酶切位点和理由的字典。
    """
    insert_sequence = read_fasta(insert_sequence_path)
    # 调用 pick_enzyme_pairs_from_dna 函数
    enzymes = pick_enzyme_pairs_from_dna(vector_dna_path, insert_sequence)
    return {
        "available_enzymes": [{"name": e["name"], "site": e["site"], "pos0": e["pos0"]} for e in enzymes],
        "reason": "These enzymes are unique cutters within the MCS of the vector and do not cut within the insert."
    }

@tool()
def design_primers(cds_sequence: str, forward_enzyme_sequence: str, reverse_enzyme_sequence: str) -> dict:
    """
    根据CDS序列文件路径设计引物。

    Args:
        cds_sequence: CDS序列FASTA文件的路径。
        forward_enzyme_sequence: 正向引物酶切位点（DNA序列，例如 "GAATTC"）。请提供精确的DNA识别序列，区分大小写。
        reverse_enzyme_sequence: 反向引物酶切位点（DNA序列，例如 "CTCGAG"）。请提供精确的DNA识别序列，区分大小写。

    Returns:
        包含正向和反向引物信息的字典。
    """
    cds_content = read_fasta(cds_sequence) # Renamed variable for clarity
    return design_primers_logic(cds_content, forward_enzyme_sequence, reverse_enzyme_sequence)

@tool()
def plan_pcr(
    template_path: str,
    forward_primer: str,
    reverse_primer: str,
    forward_primer_tm: float,
    reverse_primer_tm: float
) -> str:
    """
    设计高保真PCR扩增方案。

    Args:
        template_path: cDNA模板FASTA文件的路径。
        forward_primer: 正向引物。
        reverse_primer: 反向引物。
        forward_primer_tm: 正向引物的Tm值。
        reverse_primer_tm: 反向引物的Tm值。

    Returns:
        PCR扩增方案的文本描述。
    """
    template_sequence = read_fasta(template_path)
    product_length_bp = len(template_sequence)

    # Calculate Annealing Temperature (Ta)
    # Ta = min(Tm_f, Tm_r) - 5°C
    annealing_temp = min(forward_primer_tm, reverse_primer_tm) - 5.0

    # Calculate Extension Time (assuming 30s/kb for high-fidelity polymerase)
    extension_time_seconds = (product_length_bp / 1000.0) * 30.0
    # Round up to nearest 5 seconds for practical use
    extension_time_seconds = (round(extension_time_seconds / 5.0) * 5.0) if extension_time_seconds > 0 else 30.0
    
    # Format extension time for display
    if extension_time_seconds < 60:
        ext_time_str = f"{int(extension_time_seconds)}s"
    else:
        minutes = int(extension_time_seconds // 60)
        seconds = int(extension_time_seconds % 60)
        ext_time_str = f"{minutes}min {seconds}s" if seconds > 0 else f"{minutes}min"


    protocol_parts = [
        f"--- 高保真PCR扩增方案 ---",
        f"模板来源: {template_path} (长度: {product_length_bp} bp)",
        f"正向引物: {forward_primer} (Tm: {forward_primer_tm}°C)",
        f"反向引物: {reverse_primer} (Tm: {reverse_primer_tm}°C)",
        "",
        "**反应体系 (50 µL):**",
        "  - Phanta Max Master Mix: 25 µL",
        "  - 正向引物 (10 µM): 2 µL",
        "  - 反向引物 (10 µM): 2 µL",
        "  - 模板DNA (约50 ng/µL): 1 µL",
        "  - ddH2O: 补足至 50 µL",
        "",
        "**PCR程序:**",
        "1. 预变性: 95°C for 30s",
        "2. 循环 (30-35 cycles):",
        f"   - 变性: 95°C for 15s",
        f"   - 退火: {round(annealing_temp, 1)}°C for 15s (基于引物Tm动态计算)",
        f"   - 延伸: 72°C for {ext_time_str} (基于产物长度动态计算)",
        "3. 最终延伸: 72°C for 5min",
        "4. 4°C hold",
        "",
        "**注意事项:**",
        "  - 请根据实际情况调整引物和模板浓度。",
        "  - 建议进行梯度PCR优化退火温度。"
    ]
    return "\n".join(protocol_parts)

@tool()
def plan_ligation(vector_name: str, vector_sequence_path: str, insert_sequence_path: str) -> str:
    """
    制定连接（Ligation）方案。

    Args:
        vector_name: 载体名称。
        vector_sequence_path: 载体FASTA文件的路径。
        insert_sequence_path: 插入片段FASTA文件的路径。

    Returns:
        连接方案的文本描述。
    """
    vector_sequence = read_fasta(vector_sequence_path)
    insert_sequence = read_fasta(insert_sequence_path)
    
    vector_size_bp = len(vector_sequence)
    insert_size_bp = len(insert_sequence)

    # Calculate molar ratio for ligation (e.g., 1:3 vector:insert)
    # Assuming average_bp_mass = 660 g/mol/bp
    # mass_insert = (mass_vector * molar_ratio_insert * insert_size_bp) / (molar_ratio_vector * vector_size_bp)
    
    # Let's assume a target vector mass, e.g., 50 ng for 1:3 molar ratio
    target_vector_mass_ng = 50.0
    recommended_insert_mass_ng = (target_vector_mass_ng * 3 * insert_size_bp) / (1 * vector_size_bp)
    
    protocol_parts = [
        f"--- 连接（Ligation）方案 ---",
        f"载体: {vector_name} (来自 {vector_sequence_path}, 长度: {vector_size_bp} bp)",
        f"插入片段来源: {insert_sequence_path} (长度: {insert_size_bp} bp)",
        "",
        "**推荐摩尔比:** 载体 : 插入片段 = 1 : 3",
        f"**推荐用量:**",
        f"  - 载体 ({vector_name}): {round(target_vector_mass_ng, 2)} ng",
        f"  - 插入片段 (来自 {insert_sequence_path}): {round(recommended_insert_mass_ng, 2)} ng",
        "",
        "**连接反应体系 (20 µL):**",
        "  - T4 DNA Ligase Buffer (10X): 2 µL",
        "  - 载体DNA: 适量 (根据计算用量)",
        "  - 插入片段DNA: 适量 (根据计算用量)",
        "  - T4 DNA Ligase (1 U/µL): 1 µL",
        "  - ddH2O: 补足至 20 µL",
        "",
        "**反应条件:**",
        "  - 16°C for 1-16 hours (过夜连接效果更佳)",
        "  - 或 25°C for 10-30 minutes (快速连接)",
        "",
        "**注意事项:**",
        "  - 确保载体和插入片段都经过磷酸化和酶切处理。",
        "  - 连接产物可直接用于转化。"
    ]
    return "\n".join(protocol_parts)

@tool()
def generate_plasmid_map_and_protocol(plasmid_name: str, vector_name: str, insert_gene: str) -> dict:
    """
    生成质粒图谱和完整构建Protocol。

    Args:
        plasmid_name: 最终质粒名称。
        vector_name: 载体名称。
        insert_gene: 插入的基因名称。

    Returns:
        包含质粒图谱（图片）和构建Protocol（文本）的字典。
    """
    llm_prompt = f"""
    请以分子生物学专家身份，基于本次对话历史与以下参数，撰写一份简明、步骤明确的 {plasmid_name} 质粒构建实验方案（目标：将 {insert_gene} 插入 {vector_name} 构建 {plasmid_name}）。

    **务必整合对话历史中以下工具输出：**
    1. 引物设计：正/反向引物序列、结合区序列、长度、GC%、Tm值。
    2. PCR条件：“本次实验”的退火温度、延伸时间。
    3. 连接条件：“本次实验”的载体/插入片段用量（摩尔比计算）、连接体系与条件。

    **Protocol结构：**
    1. 实验目标（简述）。
    2. 实验步骤（仅保留必要环节，避免冗余说明）：
    - PCR扩增
    - 产物纯化
    - 酶切及产物纯化
    - 连接
    - 转化
    - 筛选与鉴定
    3. 在对应步骤直接嵌入相关数值与条件，不另起长段解释。

    **要求：**
    - 全程使用专业简洁的语言。务必精准、简明。
    - 仅输出实验方案，不包含材料列表或设备说明。
    """

    
    return {
        "plasmid_map_image": "这是一个质粒图谱的占位符。实时或动态生成质粒图谱的功能将在未来版本中实现。当前版本仅提供文本Protocol。",
        "protocol": llm_prompt
    }

@tool
def validate_construct(sequence: str) -> bool:
    """
    验证构建的正确性。

    Args:
        sequence: 要验证的DNA序列。

    Returns:
        是否通过验证。
    """
    # 简单的长度检查
    if len(sequence) < 100:
        return False
    return True

@tool
def recommend_cloning_strategy(cds_sequence: str) -> str:
    """
    推荐克隆策略。

    Args:
        cds_sequence: CDS序列。

    Returns:
        推荐的克隆策略。
    """
    
    return "推荐使用PCR扩增后直接克隆。"
    # elif insert_size < 5000:
    #     return "推荐使用酶切后连接的方法。"
    # else:
    #     return "推荐使用Gibson组装的方法。"