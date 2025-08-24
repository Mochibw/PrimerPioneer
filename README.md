# PrimerPioneer: 自动化分子克隆实验设计平台

本项目是为清华大学2025年8月黑客马拉松开发的演示项目，旨在构建一个可扩展的自动化分子克隆实验设计平台。项目利用大型语言模型（LLM）和一系列生物信息学工具，能够理解用户需求并自主设计、模拟和验证分子克隆实验流程。

## 核心目标

PrimerPioneer 的核心目标是成为一个智能的“虚拟分子生物学家”，能够：
- **自动化设计**：根据用户输入的目标基因和载体，自动设计引物、选择克隆策略并规划完整的实验步骤。
- **模拟与验证**：在实验设计阶段，通过一系列模拟工具（如PCR、酶切、连接等）对关键步骤进行虚拟验证，提前发现潜在问题。
- **高可扩展性**：项目架构设计灵活，可以轻松扩展以支持更多的分子克隆技术（如 Gibson Assembly, Golden Gate Assembly 等）和更复杂的生物学任务。
- **两种架构探索**：项目同时探索了两种不同的实现路径——单 Agent 架构和多 Agent 架构，为后续开发提供了丰富的参考。

## 架构设计

项目提供了两种截然不同的架构实现，以探索不同模式下的任务处理能力。

### 1. 单 Agent 架构 (主目录)

这是项目的基础实现，采用一个统一的 Agent 来调用所有注册的 MCP (Model Context Protocol) 工具。

- **工作流程**:
  1. **MCP 服务器启动**: `main.py` 文件通过 `FastMCP` 启动一个 MCP 服务器。
  2. **工具注册**: 所有位于 `tools_pool/` 和 `common_utils/` 目录下的功能函数（如 `get_cds_by_gene`, `simulate_pcr`, `design_primer_suite` 等）都被注册为 Agent 可调用的工具。
  3. **任务执行**: 用户通过客户端与 Agent 交互。Agent 在接收到任务（例如“克隆人类 MARCO 基因”）后，会自主规划步骤，并按顺序调用所需工具来完成任务。例如，它会先调用 `get_cds_by_gene` 获取序列，然后调用 `design_primer_suite` 设计引物，最后通过 `simulate_restriction_digestion` 和 `simulate_ligation` 等工具验证设计的合理性。

- **优势**:
  - **简洁直观**: 架构简单，易于理解和快速实现。
  - **开发高效**: 所有工具集中管理，便于调试和维护。

- **局限**:
  - **任务规划复杂**: 对于复杂的、多步骤的任务，单个 Agent 需要具备极强的任务规划和分解能力，容易出错。
  - **扩展性受限**: 随着工具数量的增加，Agent 的决策逻辑会变得越来越复杂，难以维护。

### 2. 多 Agent 架构 (`multi_agent_system/`)

为了解决单 Agent 架构的局限性，`multi_agent_system` 目录下的代码探索了一种基于角色的多 Agent 协作系统。该系统使用 `LangGraph` 框架构建，将复杂的分子克隆任务分解给不同角色的 Agent 协同完成。

- **核心组件**:
  - **Supervisor (监督者)**: 作为任务的“总指挥”，负责理解用户意图，并将任务分解成子任务，然后将子任务路由给最合适的专家 Agent。
  - **专家 Agents**:
    - `SequenceAnalysisAgent`: 专注于序列的获取和分析，如下载 CDS 序列、查找酶切位点等。
    - `StrategySelectionAgent`: 负责根据序列特征和任务目标，选择最合适的克隆策略（如双酶切法、Gibson Assembly 等）。
    - `GeneticComponentDesignerAgent`: 核心设计者，负责设计引物、密码子优化等，并调用模拟工具进行验证。
    - `ConstructValidationAgent`: 最终的“质检员”，负责审查整个流程的完整性和合理性，并向用户提交最终报告。

- **工作流程 (基于 LangGraph 的状态图)**:
  1. **任务接收**: 用户任务首先进入 `supervisor` 节点。
  2. **任务路由**: `supervisor` 根据当前任务状态，决定下一个应该由哪个专家 Agent 接手。
  3. **循环执行**: 专家 Agent 完成其子任务后，控制权返回给 `supervisor`，由其决定下一步行动，直到所有步骤完成。
  4. **最终审核**: 所有设计工作完成后，`supervisor` 将任务交给 `ConstructValidationAgent` 进行最终审核和结果汇总。

- **优势**:
  - **职责清晰**: 每个 Agent 只负责自己领域内的任务，降低了单一模型的复杂度。
  - **高可扩展性**: 可以方便地增加新的专家 Agent 来支持新的功能或克隆策略。
  - **鲁棒性强**: 任务分解使得系统在处理复杂问题时更加稳定和可靠。

## 实现方法

### 工具池 (`tools_pool/`)

这是项目的核心功能库，包含了所有分子生物学操作的实现。每个工具都是一个独立的 Python 函数，并有清晰的输入输出定义（遵循 `Schema.md` 规范）。

- **主要工具类别**:
  - **序列获取与分析**: `get_cds_by_gene`, `get_sequence_info`, `find_features`
  - **克隆策略**: `recommend_cloning_strategy`, `select_cloning_strategy`
  - **引物设计**: `design_primer_suite`
  - **实验模拟**: `simulate_pcr`, `simulate_restriction_digestion`, `simulate_ligation`, `simulate_end_repair` 等。
  - **结果可视化**: `generate_map`

### 数据管理 (`data/` & `common_utils/`)

- `data/` 目录用于存放用户上传的序列文件（如 `.fasta`, `.gb`）以及程序生成的中间文件和结果文件。
- `common_utils/` 包含文件操作和序列处理的通用工具函数，确保数据格式的统一和正确传递。

## 后续开发方向

本项目为未来的扩展奠定了坚实的基础，以下是一些建议的开发方向：

1.  **完善多 Agent 架构**:
    - **丰富专家角色**: 增加更多专家 Agent，例如专门负责密码子优化的 `CodonOptimizerAgent`，或负责 Golden Gate / Gibson Assembly 方案设计的 `AdvancedCloningAgent`。
    - **优化 Supervisor 逻辑**: 提升 Supervisor 对复杂任务的理解和分解能力，使其能够处理更具挑战性的设计任务。

2.  **扩展工具池**:
    - **支持更多克隆方法**: 实现 Golden Gate Assembly、Gibson Assembly、CRISPR 等现代分子克隆技术的相关工具。
    - **整合生物学数据库**: 对接更多外部数据库（如 Addgene, UniProt），以获取更丰富的质粒和蛋白信息。

3.  **增强用户交互**:
    - **图形化界面**: 开发一个前端界面，让用户可以更直观地输入任务、查看进度和最终的可视化结果。
    - **交互式设计**: 允许用户在设计的关键节点进行干预和调整，实现人机协同。

4.  **引入 RAG (Retrieval-Augmented Generation)**:
    - 目前 `recommend_cloning_strategy` 已初步尝试使用 RAG，未来可以构建更全面的分子克隆知识库，让 Agent 能够查询最新的实验方案和文献，从而给出更科学、更前沿的设计建议。

---
*本项目由清华大学2025年8月黑客马拉松 PrimerPioneer 团队开发。*
