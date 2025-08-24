# rag_lines.py
# Chroma + Alibaba Cloud text-embedding | 一行一块 的最小RAG（函数版，无API）
# 用法示例见文件末尾 __main__ 区域

from __future__ import annotations
import os
from typing import List, Tuple, Dict, Any, Optional
import chromadb
from http import HTTPStatus
import dashscope
from openai import OpenAI
from langchain.tools import tool

# Set up Alibaba Cloud API key
dashscope.api_key = os.getenv("DASHSCOPE_API_KEY") or "sk-e5070c7389b347f0b99730031e61fcfd"
DB_FILE_PATH = r"D:\workspace\python\primepioneer\src\tools_pool\data\rag_db"
DEMO_FILE_PATH = r"D:\workspace\python\primepioneer\src\tools_pool\v0\tools_pool\strategy.md"

# Set up OpenAI client
openai_client = OpenAI(
    base_url="https://openrouter.ai/api/v1",  # Default for gpt-oss
    api_key=os.getenv("OPENROUTER_API_KEY"),
)

def ali_embed(texts: list[str]) -> list[list[float]]:
    # 阿里云一次可批量，texts 为字符串列表
    # 限制批次大小为10
    embeddings = []
    batch_size = 10
    for i in range(0, len(texts), batch_size):
        batch = texts[i:i + batch_size]
        rsp = dashscope.TextEmbedding.call(
            model="text-embedding-v4",
            input=batch
        )
        if rsp.status_code == HTTPStatus.OK:
            embeddings.extend([item["embedding"] for item in rsp.output["embeddings"]])
        else:
            raise RuntimeError(f"Aliyun embedding error: {rsp.code} {rsp.message}")
    return embeddings

def optimize_search_query(query: str) -> str:
    """使用LLM将自然语言查询转换为优化的搜索词"""
    prompt = f"""
    请将以下关于分子克隆策略的自然语言查询转换为适合向量检索的优化搜索词。
    只需提供优化后的搜索词，不需要其他解释或格式。
    
    原始查询: {query}
    
    优化后的搜索词:"""
    
    try:
        response = openai_client.chat.completions.create(
            model="openai/gpt-5-mini",
            messages=[
                {"role": "system", "content": "你是一个专业的分子生物学研究助手。"},
                {"role": "user", "content": prompt}
            ],
            temperature=0.3,
            max_tokens=100
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        # 如果LLM调用失败，返回原始查询
        return query

def summarize_results(query: str, results: List[Dict[str, Any]]) -> str:
    """使用LLM根据检索结果生成自然语言回答"""
    # 准备检索到的文本内容
    context = "\n".join([f"{i+1}. {result['text']}" for i, result in enumerate(results)])
    
    prompt = f"""
    请根据以下分子克隆策略相关的检索结果，回答用户的问题。回答应简洁明了，专业准确。
    
    用户问题: {query}
    
    检索结果:
    {context}
    
    请根据上述检索结果回答用户问题。请注意，只给出解决问题的思路，不要包含任何具体的序列、酶等例子，避免误导用户。你只可以推荐如下注册过的工具，
    (load_sequence) # 从.dna和.fasta等文件加载序列到/tools_pool/Schema.md定义的json格式
    (get_sequence_info) # 获取序列的基本信息，如长度、feature等，不包含序列字符串本身
    (find_features) # 查找序列中的特征（features，包括酶切位点和自定义features），返回符合条件的特征列表
    (design_primer_suite) # 设计引物套件（Primer Suite），返回引物列表
    (simulate_restriction_digestion) # 模拟限制性内切酶消化，返回消化后的片段列表
    (simulate_pcr) # 模拟PCR扩增，返回扩增后的片段列表
    (simulate_ligation) # 模拟DNA片段连接，返回连接产物列表
    (write_record_to_json) # 将记录写入JSON文件
    (load_sequence_from_json) # 从JSON文件加载SequenceRecord等符合Schema的记录
    (list_data) # 列出 /data/ 目录下的文件树结构
    (recommend_cloning_strategy) # 基于RAG的克隆策略推荐。
    请给出你的建议：
    """
    
    try:
        response = openai_client.chat.completions.create(
            model="openai/gpt-5-mini",
            messages=[
                {"role": "system", "content": "你是一个专业的分子生物学研究助手。"},
                {"role": "user", "content": prompt}
            ],
            temperature=0.7,
            max_tokens=5000
        )
        return response.choices[0].message.content.strip()
    except Exception as e:
        # 如果LLM调用失败，返回基于检索结果的简单总结
        print(e)
        return f"根据检索到的信息，相关的克隆策略建议包括：\n" + "\n".join([f"- {result['text']}" for result in results[:3]])

@tool()
def recommend_cloning_strategy(query: str) -> str:
    """
    将自然语言的克隆相关问题转化为优化检索并输出总结性策略建议。

    Args:
        query (str): 自然语言描述的克隆需求或问题。

    Returns:
        str: 针对输入问题的中文克隆策略建议；若无相关结果则返回提示信息。
    """
    # 初始化RAG系统
    rag = LineRAG()
    
    # 第一次LLM调用：优化查询词
    optimized_query = optimize_search_query(query)
    
    # 执行向量搜索
    search_results = rag.search(optimized_query, k=5)
    
    
    # 第二次LLM调用：总结结果
    answer = summarize_results(query, search_results)
    return answer

# -------------------------
# 初始化 & 全局资源
# -------------------------
class LineRAG:
    def __init__(
        self,
        db_path: str = DB_FILE_PATH,
        collection: str = "doc_lines",
        model_name: str = "text-embedding-v4",   # Alibaba Cloud text-embedding model
        distance: str = "cosine"                           # "cosine" | "l2" | "ip"
    ):
        # No need to initialize embed model for Alibaba Cloud
        self.db = chromadb.PersistentClient(path=db_path)
        self.col = self.db.get_or_create_collection(
            collection, metadata={ "hnsw:space": distance }
        )

    # -------------------------
    # 更新数据库入口（写入/更新）
    # -------------------------
    def upsert_texts(
        self,
        texts: List[str],
        ids: Optional[List[str]] = None,
        metadatas: Optional[List[Dict[str, Any]]] = None
    ) -> None:
        """批量写入/更新文本。ids 不给则自动编号。"""
        texts = [t for t in texts if isinstance(t, str) and t.strip()]
        if not texts:
            return
        if ids is None:
            ids = [f"item_{i}" for i in range(len(texts))]
        if metadatas is None:
            metadatas = [{} for _ in texts]
        embs = ali_embed(texts)
        self.col.upsert(documents=texts, ids=ids, embeddings=embs, metadatas=metadatas)

    def upsert_file_by_lines(
        self,
        file_path: str,
        encoding: str = "utf-8",
        strip_blank: bool = True
    ) -> int:
        """将文件逐行入库；返回写入的行数。ID = {basename}:{line_no}。"""
        if not os.path.exists(file_path):
            raise FileNotFoundError(file_path)
        base = os.path.basename(file_path)
        with open(file_path, encoding=encoding) as f:
            lines = f.readlines()
        if strip_blank:
            lines = [ln.strip() for ln in lines if ln.strip()]
        ids = [f"{base}:{i}" for i in range(len(lines))]
        metas = [ {"file": base, "line_no": i} for i in range(len(lines)) ]
        self.upsert_texts(lines, ids=ids, metadatas=metas)
        return len(lines)

    # -------------------------
    # 查询入口（检索）
    # -------------------------
    def search(
        self,
        query: str,
        k: int = 5,
        include_ids: bool = True,
        include_metas: bool = True
    ) -> List[Dict[str, Any]]:
        """向量检索：返回 top-k 行及相似度分数（1 - 距离）。"""
        q_emb = ali_embed([query])
        include = ["documents", "distances"]
        if include_metas: include.append("metadatas")
        r = self.col.query(query_embeddings=q_emb, n_results=max(1, k), include=include)

        docs   = r.get("documents",  [[]])[0]
        ids    = r.get("ids",        [[]])[0]
        metas  = r.get("metadatas",  [[]])[0]
        dists  = r.get("distances",  [[]])[0]
        out = []
        for idx, doc in enumerate(docs):
            item = {
                "text": doc,
                "score": 1 - dists[idx] if idx < len(dists) else None
            }
            if include_ids and idx < len(ids):
                item["id"] = ids[idx]
            if include_metas and idx < len(metas):
                item["meta"] = metas[idx]
            out.append(item)
        return out

    # -------------------------
    # 维护工具（可选）
    # -------------------------
    def count(self) -> int:
        """估算条目数量（Chroma 无直接 count，这里用 ids 近似）。"""
        r = self.col.get(include=[])
        return len(r.get("ids", []))

    def delete_by_file(self, file_basename: str) -> None:
        """按文件名删除其所有行（基于 ID 前缀匹配）。"""
        # 取出所有 ids 后过滤（大集合上会慢；数据量大时请自行维护映射）
        all_ids = self.col.get(include=[]).get("ids", [])
        to_del = [i for i in all_ids if i.startswith(f"{file_basename}:")]
        if to_del:
            self.col.delete(ids=to_del)

# -------------------------
# 示例用法
# -------------------------
if __name__ == "__main__":
    rag = LineRAG(
        db_path=DB_FILE_PATH,
        collection="doc_lines",
        model_name="text-embedding-v4",   # Alibaba Cloud text-embedding model
        distance="cosine"
    )

    # === 1) 逐行导入：把你的文件路径替换为本地实际路径 ===
    demo_files = [
        DEMO_FILE_PATH
    ]
    for fp in demo_files:
        if os.path.exists(fp):
            n = rag.upsert_file_by_lines(fp)
            print(f"[upsert] {fp} -> {n} lines")
        else:
            print(f"[skip]   {fp} (not found)")

    print(f"[count] collection size ≈ {rag.count()} lines")

    # === 3) 测试推荐克隆策略函数 ===
    test_queries = [
        # "克隆 data 目录下的 RTCB 基因到表达载体 PCDNA3.1 的最佳策略，包括酶切位点选择、引物设计和连接方法"
        "克隆 data 目录下的 RTCB 基因到表达载体 PCDNA3.1 的最佳策略，包括酶切位点选择、引物设计和连接方法"
    ]
    for q in test_queries:
        print(f"\n[recommend] {q}")
        recommendation = recommend_cloning_strategy(q)
        print(f"  {recommendation}")
