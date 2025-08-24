from langchain.tools import tool

@tool
def select_cloning_strategy(cds_sequence: str) -> str:
    """
    推荐克隆策略。

    Args:
        cds_sequence: CDS序列。

    Returns:
        推荐的克隆方案（酶切/Gibson克隆）。
    """
    
    return "推荐使用酶切克隆方案。"