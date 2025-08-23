from __future__ import annotations
import os, json, re, sys
from typing import TypedDict, Optional, Literal, Annotated

from langgraph.graph import StateGraph, START, END
from langgraph.graph.message import AnyMessage, add_messages
from langchain_core.messages import HumanMessage, AIMessage

from langchain_deepseek import ChatDeepSeek
from langchain_community.chat_models import ChatTongyi
from langchain_openai import ChatOpenAI
from langgraph.prebuilt import create_react_agent

sys.path.append(r"D:\workspace\python\hello_MCP\multi_agent_system\src")
sys.path.append(r"D:\workspace\python\hello_MCP\multi_agent_system\src\tools_pool\HW1_demo")
from tools_pool.HW1_demo.tools import get_cds_sequence, read_fasta_file, select_restriction_sites, recommend_cloning_strategy, design_primers, plan_pcr, plan_ligation


# =========================
# Config
# =========================
MAX_ITERS = 3

# llm = ChatOpenAI(
#     model="openai/gpt-5-mini",
#     api_key=os.getenv("OPENROUTER_API_KEY"),
#     base_url="https://openrouter.ai/api/v1"
#     # temperature=0.2,              # 可调
#     # max_output_tokens=1024,       # 可调
#     )

llm = ChatTongyi(
    model="qwen2.5-72b-instruct",   # 常见：qwen2.5-32b-instruct / qwen2.5-7b-instruct / qwen-turbo 等
    # temperature=0.2,
    # max_tokens=1024,
)


# =========================
# 各 Agent（用 create_react_agent 包一层，拥有各自工具池）
# =========================
SequenceAnalysisAgent = create_react_agent(
    model=llm,
    tools=[get_cds_sequence, read_fasta_file, select_restriction_sites],
    prompt=(
        "你是一个善于获取并分析cds序列的agent，可以完成cds序列读取、根据dna载体与cds序列进行双酶切位点选择等任务。"
        "选择好酶切位点之后就可以返回结果了，不要尝试进行后续的克隆策略选取、PCR方案设计、引物设计等工作。"
        "你可以在有帮助的情况下使用工具。"
    ),
)

StrategySelectionAgent = create_react_agent(
    model=llm,
    tools=[recommend_cloning_strategy],
    # prompt=(
    #     "你是一个善于推荐基因克隆策略的agent，可以完成基因克隆策略的推荐任务。"
    #     "选择好克隆策略（PCR/酶切/gibson等）之后即可返回结果，不要尝试后续的PCR方案设计、引物设计等工作。"
    #     "你可以在有帮助的情况下使用工具。"
    # ),
    prompt=(
        "你是一个善于推荐基因克隆策略的agent，可以完成基因克隆策略的推荐任务。"
        "选择好克隆策略（PCR/酶切/gibson等）之后即可返回结果，不要尝试后续的PCR方案设计、引物设计等工作。"
        "你可以在有帮助的情况下使用工具。"
    ),
)

GeneticComponentDesignerAgent = create_react_agent(
    model=llm,
    tools=[design_primers, plan_pcr, plan_ligation],
    prompt=(
        "你是一个善于设计基因组件的agent，可以完成引物设计、PCR扩增方案设计、连接方案设计等任务。"
        "你可以在有帮助的情况下使用工具。"
    ),
)

ConstructValidationAgent = create_react_agent(
    model=llm,
    tools=[],
    prompt=(
        "你是一个善于评估基因克隆策略的agent\n"
        "回顾对话和最新的助手草稿（如果有），如果结果中包括具体的引物设计、PCR扩增方案、连接方案，可以同意返回结果，不需要验证具体方案是否合理\n"
        "在有帮助的情况下使用工具\n"
        "仅返回一个包含以下字段的JSON对象：\n"
        '{\n'
        '  "approved": true|false,\n'
        '  "feedback": "如果你有批评意见或需要修复的地方，请指出",\n'
        '  "final_answer": "如果同意，给出最终回复；否则输出空回复"\n'
        "}\n"
        "不要添加任何除JSON对象以外的其他文本。"
    ),
)

# =========================
# 监督者（Router）：决定下一个 worker 或进入 evaluate
# =========================
SUPERVISOR_SYS = (
    "你是一个监督者。检查任务和当前进展。\n"
    "决定下一步：\n"
    '返回一个仅包含以下内容的JSON对象：{"next": "SequenceAnalysisAgent|StrategySelectionAgent|GeneticComponentDesignerAgent|ConstructValidationAgent", "reason": "简短的选择理由"}\n'
    "建议第一步进行cds序列分析，第二步选择克隆策略（不包括具体的基因组件设计），第三步根据cds序列分析结果与克隆策略进行具体的基因组件（引物、PCR扩增方案、连接方案等）设计，最后交给ConstructValidationAgent生成质粒图谱和完整构建Protocol。**不要**在没有生成具体的引物设计、PCR扩增方案与连接方案时把结果交给ConstructValidationAgent*"
    "选择一个能够最大化进展的工作代理。如果你觉得当前回复已经完成任务要求，请选择 'ConstructValidationAgent'。"
)

def parse_json_like(s: str, fallback: dict) -> dict:
    """从文本中抽取第一个 JSON；不行则返回 fallback。"""
    try:
        return json.loads(s)
    except Exception:
        m = re.search(r"\{.*\}", s, flags=re.S)
        if m:
            try:
                return json.loads(m.group(0))
            except Exception:
                return fallback
    return fallback

def supervisor_decide(messages: list[AnyMessage]) -> dict:
    """调用 LLM 让 supervisor 产生路由 JSON。"""
    resp = llm.invoke([{"role": "system", "content": SUPERVISOR_SYS}, *messages])
    data = parse_json_like(
        resp.content,
        {"next": "ConstructValidationAgent", "reason": "fallback-route"},
    )
    nxt = data.get("next", "ConstructValidationAgent")
    if nxt not in {"SequenceAnalysisAgent", "StrategySelectionAgent", "GeneticComponentDesignerAgent", "ConstructValidationAgent"}:
        nxt = "ConstructValidationAgent"
    return {"next": nxt, "reason": data.get("reason", "")}


# =========================
# LangGraph State 定义
# =========================
class GraphState(TypedDict):
    messages: Annotated[list[AnyMessage], add_messages]
    iterations: int
    route: Optional[Literal["SequenceAnalysisAgent", "StrategySelectionAgent", "GeneticComponentDesignerAgent", "ConstructValidationAgent"]]
    approved: Optional[bool]
    feedback: Optional[str]
    final_answer: Optional[str]


# =========================
# 节点实现
# =========================
def node_supervisor(state: GraphState) -> dict:
    # 如果已达到上限，强制进入评审
    if state["iterations"] >= MAX_ITERS:
        return {"route": "ConstructValidationAgent"}

    decision = supervisor_decide(state["messages"])
    # 把 supervisor 的决定记录到对话里（可选）
    sup_msg = AIMessage(content=f'[Supervisor] next={decision["next"]}; reason={decision["reason"]}')
    return {"messages": [sup_msg], "route": decision["next"]}

def _run_worker(agent, name: str):
    def _node(state: GraphState) -> dict:
        res = agent.invoke({"messages": state["messages"]})
        # create_react_agent 返回 {"messages": [...]}
        msgs = res.get("messages", [])
        return {"messages": msgs}
    _node.__name__ = f"node_{name}"
    return _node

node_SequenceAnalysisAgent = _run_worker(SequenceAnalysisAgent, "SequenceAnalysisAgent")
node_StrategySelectionAgent = _run_worker(StrategySelectionAgent, "StrategySelectionAgent")
node_GeneticComponentDesignerAgent = _run_worker(GeneticComponentDesignerAgent, "GeneticComponentDesignerAgent")

def node_ConstructValidationAgent(state: GraphState) -> dict:
    res = ConstructValidationAgent.invoke({"messages": state["messages"]})
    msgs = res.get("messages", [])
    eval_text = msgs[-1].content if msgs else "{}"
    obj = parse_json_like(eval_text, {"approved": False, "feedback": "Invalid JSON", "final_answer": ""})
    # 把评审 JSON 也加入消息流，便于审计
    eval_msg = AIMessage(content=f'[Evaluate] {json.dumps(obj, ensure_ascii=False)}')
    new_iters = state["iterations"] + 1 if not obj.get("approved", False) else state["iterations"]
    return {
        "messages": [eval_msg],
        "approved": bool(obj.get("approved", False)),
        "feedback": obj.get("feedback", ""),
        "final_answer": obj.get("final_answer", ""),
        "iterations": new_iters,
    }


# =========================
# 构图
# =========================
graph = StateGraph(GraphState)

graph.add_node("supervisor", node_supervisor)
graph.add_node("SequenceAnalysisAgent", node_SequenceAnalysisAgent)
graph.add_node("StrategySelectionAgent", node_StrategySelectionAgent)
graph.add_node("GeneticComponentDesignerAgent", node_GeneticComponentDesignerAgent)
graph.add_node("ConstructValidationAgent", node_ConstructValidationAgent)

# 流程：START -> supervisor -> (worker_x or evaluate)
graph.add_edge(START, "supervisor")

def route_from_supervisor(state: GraphState) -> str:
    return state.get("route") or "ConstructValidationAgent"

graph.add_conditional_edges(
    "supervisor",
    route_from_supervisor,
    {
        "SequenceAnalysisAgent": "SequenceAnalysisAgent",
        "StrategySelectionAgent": "StrategySelectionAgent",
        "GeneticComponentDesignerAgent": "GeneticComponentDesignerAgent",
        "ConstructValidationAgent": "ConstructValidationAgent",
    },
)

# worker 执行完回到 supervisor
graph.add_edge("SequenceAnalysisAgent", "supervisor")
graph.add_edge("StrategySelectionAgent", "supervisor")
graph.add_edge("GeneticComponentDesignerAgent", "supervisor")

def route_from_evaluate(state: GraphState) -> str:
    if state.get("approved"):
        return "approve"
    if state.get("iterations", 0) >= MAX_ITERS:
        return "maxed"
    return "revise"

graph.add_conditional_edges(
    "ConstructValidationAgent",
    route_from_evaluate,
    {
        "approve": END,
        "maxed": END,
        "revise": "supervisor",
    },
)

app = graph.compile()
# Visualize the graph
# For Jupyter or GUI environments:
app.get_graph().draw_mermaid_png()

# To save PNG to file:
png_data = app.get_graph().draw_mermaid_png()
with open("graph.png", "wb") as f:
    f.write(png_data)


# =========================
# 运行示例
# =========================
if __name__ == "__main__":
    user_task = "我需要克隆marco基因，请你给出完整实验流程。dna载体文件路径r'D:\workspace\python\hello_MCP\PrimerPioneerSetup\pcDNA3.1(-).dna'。"
    init_state: GraphState = {
        "messages": [HumanMessage(content=user_task)],
        "iterations": 0,
        "route": None,
        "approved": None,
        "feedback": None,
        "final_answer": None,
    }

    # A) 一次性执行
    # final_state = app.invoke(init_state)
    # print("\n=== Conversation tail ===")
    # for m in final_state["messages"][-6:]:
    #     role = "USER" if m.type == "human" else "AI"
    #     print(f"{role}: {m.content}")

    # if final_state.get("approved"):
    #     print("\n=== FINAL ANSWER ===")
    #     print(final_state.get("final_answer", ""))
    # else:
    #     print("\n=== NOT APPROVED ===")
    #     print("feedback:", final_state.get("feedback"))

    # B) 流式查看进度（可选）
    with open("agent_output.txt", "w", encoding="utf-8") as f:
        for ev in app.stream(init_state, stream_mode=["updates"]):
            f.write(str(ev) + "\n")   # 每个事件转成字符串
            print(">>", ev)
