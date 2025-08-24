# main.py
from mcp.server.fastmcp import FastMCP
from tools_pool.get_sequence_info import get_sequence_info
from tools_pool.find_features import find_features
from common_utils.file_operations import load_sequence
from tools_pool.design_primer_suite import design_primer_suite
from tools_pool.simulate_restriction_digestion import simulate_restriction_digestion
from tools_pool.simulate_pcr import simulate_pcr
from common_utils.file_operations import write_record_to_json, load_sequence_from_json
from common_utils.sequence import get_sequence
from tools_pool.simulate_ligation import simulate_ligation
from tools_pool.simulate_end_repair import simulate_end_repair
from tools_pool.simulate_a_tailing import simulate_a_tailing, simulate_t_overhang
from tools_pool.simulate_oligo_annealing import simulate_oligo_annealing
from tools_pool.simulate_phosphorylation import simulate_phosphorylation
from tools_pool.simulate_dna_synthesis import simulate_dna_synthesis
from tools_pool.simulate_homology_assembly import simulate_homology_assembly
from tools_pool.generate_map import generate_map
from common_utils.file_operations import list_data
from tools_pool.strategy_query import recommend_cloning_strategy
from tools_pool.get_cds import get_cds_by_gene
from tools_pool.simulate_gel_purification import simulate_gel_purification

# 初始化 FastMCP 服务器（注意：stdio 模式下不要向 stdout 打印任意文本）
mcp = FastMCP("my-mcp-server")

# 从Tools Pool注册工具
mcp.tool()(load_sequence) # 从.dna和.fasta等文件加载序列到/tools_pool/Schema.md定义的json格式
mcp.tool()(get_sequence_info) # 获取序列的基本信息，如长度、feature等，不包含序列字符串本身
mcp.tool()(find_features) # 查找序列中的特征（features，包括酶切位点和自定义features），返回符合条件的特征列表
mcp.tool()(design_primer_suite) # 设计引物套件（Primer Suite），返回引物列表
mcp.tool()(simulate_restriction_digestion) # 模拟限制性内切酶消化，返回消化后的片段列表
mcp.tool()(simulate_pcr) # 模拟PCR扩增，返回扩增后的片段列表
#mcp.tool()(write_record_to_json) # 将记录写入JSON文件
#mcp.tool()(load_sequence_from_json) # 从JSON文件加载SequenceRecord等符合Schema的记录
#mcp.tool()(get_sequence) # 从包含 "sequence" 字段的 JSON 文件中读取序列，并返回序列字符串本身
mcp.tool()(simulate_ligation) # 模拟DNA连接，返回连接后的片段
mcp.tool()(simulate_end_repair) # 末端修复，将粘性末端修复为平末端
mcp.tool()(simulate_a_tailing) # 加A尾
mcp.tool()(simulate_t_overhang) # 加T
mcp.tool()(simulate_oligo_annealing) # 模拟寡核苷酸退火
mcp.tool()(simulate_phosphorylation) # DNA末端加磷酸或去磷酸
mcp.tool()(simulate_dna_synthesis) # DNA化学合成
mcp.tool()(simulate_homology_assembly) # 同源重组组装
mcp.tool()(list_data) # 列出 /data/ 目录下的文件树结构
mcp.tool()(recommend_cloning_strategy) # 基于RAG的克隆策略推荐
mcp.tool()(generate_map) # 生成序列的可视化图谱
mcp.tool()(get_cds_by_gene) # 从NCBI获取基因CDS序列，需要明确指定物种
mcp.tool()(simulate_gel_purification) # 模拟胶回收，从酶切产物中选择特定片段

def main():
    """启动 MCP 服务器"""
    # 显式指定 stdio 传输
    mcp.run(transport="stdio")

if __name__ == "__main__":
    main()
