from typing import List, Dict, Optional, Literal, TypedDict, Tuple, Union

# ------------------------------
# 共享数据结构 (JSON 安全)
# ------------------------------

# 工具开发规范
# 1. 所有函数都要满足本页的Schema; Schema应该在所有工具开发前统一，如果开发后出现改动，要考虑所有下游联动函数的修改。
# 2. 由于LLM调用是“无状态”的，即不会有保存在内存中的变量，所以函数的输入输出都应该是json文件的路径，读入路径先解析到变量，操作完成时保存到json文件，并return对应的文件路径。
# 3. 建议所有的return除了计划返回的内容之外，加一个info项，可以作为LLM的局部提示词（比如遇到某个错误，返回可能的错误原因，或指导LLM，一般下一步应该调用什么工具。）
# 4. 目前cherryStudio中，工具返回结果只会出现在单轮LLM调用的上下文中，而不会保存历史记录。如果LLM没有“口述”一遍工具调用结果，下一轮对话无法再次看到工具返回值。（不优雅的解决方案之一是，在系统提示词强制让LLM复述info信息）
# 5. 超过200bp的序列尽量不要让LLM调用，序列读取工具，可以默认只开放head()和tool()，除非显式指定full_length，并在DocsString建议LLM如无必要不要读取full_length。如果return序列信息，手动限制返回长度（示例：/common_untils/sequence.py里的get_sequence函数）
# 6. 每个tool单独一个py文件，放在/tools_pool/文件夹中，每个Agent各取所需，通过from * import *注册（当前实例为单Agent，在main.py注册），并添加mcp.too()装饰器
# 7. 为了实现每个tool的主函数的辅助函数，如果只被当前tool用到，就写在该工具自己的py文件中；如果可复用，就放在/common_utils文件夹中，如load_sequence_from_json()和write_record_to_json()等。
Strand = Literal[1, -1]
EndType = Literal["blunt", "5_overhang", "3_overhang"]

class SequenceRecord(TypedDict, total=False):
    id: str                    # 记录的唯一 ID
    name: str                  # 人类可读的名称
    sequence: str              # 大写 A/C/G/T (可选 IUPAC 字符: RYKMSWBDHVN)
    length: int                # sequence 的长度
    circular: bool             # 是否为环状质粒
    features: List["Feature"]  # 可选的特征注释
    metadata: Dict[str, str]   # 自由格式的键值对元数据

class Feature(TypedDict, total=False):
    type: str                  # 例如: 'CDS', 'promoter', 'ori', 'MCS', 'polyA', 'tag'
    start: int                 # 起始位置 (1-based, 包含)
    end: int                   # 结束位置 (1-based, 包含)
    strand: Strand             # 链方向 (1 或 -1)
    qualifiers: Dict[str, Union[str, float, int]] # 限定符

class FragmentEnd(TypedDict, total=False):
    kind: EndType              # 'blunt' (平末端) | '5_overhang' (5' 突出) | '3_overhang' (3' 突出)
    seq: str                   # 突出端的序列 (平末端时为空字符串)
    length: int                # 突出端的长度 (平末端时为 0)

class Fragment(TypedDict, total=False):
    id: str
    start: int
    end: int
    length: int
    strand: Strand
    overhang_5: FragmentEnd
    overhang_3: FragmentEnd
    sequence: Optional[str]    # 对于较大的片段可省略，以节省空间

class EnzymeSpec(TypedDict, total=False):
    name: str                  # 例如: 'EcoRI'
    site: str                  # 识别位点, 例如: 'GAATTC'
    cut_index_top: int         # 正链上的切割位置 (相对于识别位点, 0-based)
    cut_index_bottom: int      # 负链上的切割位置 (相对于识别位点, 0-based)
    star_activity: bool        # 是否有星号活性 (可选)

class DigestResult(TypedDict, total=False):
    enzymes: List[str]
    cuts: List[int]             # 切割位点列表 (1-based, 表示切口左侧碱基的索引)
    fragments: List[Fragment]
    info: List[str]             # 关于酶切效率等的警告信息

class Amplicon(TypedDict, total=False):
    id: str
    start: int
    end: int
    length: int
    sequence: str
    primers: Dict[str, "Primer"]  # {'forward':..., 'reverse':...}

class Primer(TypedDict, total=False):
    name: str
    sequence: str
    tm: float
    gc: float
    start: Optional[int]       # 结合起始位点 (1-based, 基于正链坐标)
    end: Optional[int]         # 结合结束位点 (1-based)
    direction: Literal["F", "R"]

class AssemblyStep(TypedDict, total=False):
    action: str                # 例如: 'PCR', 'Digest', 'Phosphorylation', 'Ligation', 'Gibson'
    inputs: List[str]          # 消耗的片段/序列的 ID
    outputs: List[str]         # 产生的片段/序列的 ID
    params: Dict[str, Union[str, int, float, bool]]

class AssemblyPlan(TypedDict, total=False):
    method: str                # 'restriction_ligation' | 'gibson' | 'isothermal' | 'TA' | 'custom'
    steps: List[AssemblyStep]
    rationale: str
    expected_product: SequenceRecord

class ValidationIssue(TypedDict, total=False):
    severity: Literal["info", "warning", "error"]
    code: str
    message: str
    region: Optional[Tuple[int, int]]  # (1-based, 包含)

class ValidationReport(TypedDict, total=False):
    ok: bool
    issues: List[ValidationIssue]
    summary: str

class MethodRecommendation(TypedDict, total=False):
    recommended_method: Literal["pcr", "synthesis", "other"]
    confidence: Literal["high", "medium", "low"]
    reasoning: str
    next_action: str
    critical_issues: List[str]

class SequenceAnalysis(TypedDict, total=False):
    length: int
    gc_content: float
    pcr_feasibility: Dict
    repetitiveness: bool
    homopolymers: Dict
    scores: Dict
    critical_issues: List[str]
    
# ------------------------------
# 已实现的工具
# ------------------------------

def get_sequence_info(path: str) -> Dict[str, Any]:
    """
    读取 SequenceRecord JSON 文件并返回除序列外的所有元数据。

    实现说明:
    - 这是一个辅助函数，用于读取一个假定包含 SequenceRecord 的 JSON 文件。
    - 在返回之前，它会移除记录中的 'sequence' 字段，以提供一个摘要信息。
    
    Args:
        path (str): 包含 SequenceRecord 数据的 JSON 文件路径。

    Returns:
        Dict[str, Any]: 不含 'sequence' 字段的 SequenceRecord 信息。
    """
    # 实现于: tools_pool/get_sequence_info.py
    pass


def find_features(json_path: str,
                  scan_builtin: bool = True,
                  custom_patterns: Optional[List[Dict[str, str]]] = None) -> List[Feature]:
    """
    在给定序列中识别功能位点。

    实现说明:
    - 输入是 `json_path` (一个 SequenceRecord 文件的路径)，而不是一个 SequenceRecord 对象。
    - 它会保留输入文件中的已有特征。
    - 内置库包含常见的限制性酶切位点 (EcoRI, BamHI 等)、T7 启动子和 polyA 信号。
    - 支持通过 `custom_patterns` 定义用户自定义的模式。
    - 扫描会在正链和反链上同时进行。
    
    Args:
        json_path: 输入的 SequenceRecord JSON 文件路径。
        scan_builtin: 若为 True，则启用内置库进行扫描。
        custom_patterns: 可选的自定义模式列表，例如:
            {
              "name": "MySite",
              "type": "misc_feature",
              "pattern": "GATTACA"
            }
    
    Returns:
        List[Feature]: 所有找到的特征列表 (包括已有的和新发现的)。
    """
    # 实现于: tools_pool/find_features.py
    pass


def design_primer_suite(target: str,
                        task: Literal["amplify", "mutagenesis"] = "amplify",
                        intent: Optional[Dict] = None) -> Dict[str, Any]:
    """
    根据指定的任务设计一对引物。

    实现说明:
    - `target` 参数可以是一个原始序列字符串、一个文件路径 (.dna, .fasta) 或一个 SequenceRecord JSON 文件的路径。
    - 详细的设计参数 (如 Tm, 长度等) 不对用户开放，已在内部固定
      (Tm: 58-64°C, 长度: 18-28 bp, 优先但不强制要求 3' 端为 G/C)。
    - 所有任务相关的细节都通过 `intent` 字典传递。
    - 函数返回一个包含 `status`, `message`, 和 `data` 字段的字典，以便进行健壮的错误处理。
    
    Args:
        target: 模板 DNA 序列 (以字符串或文件路径形式提供)。
        task: 设计任务，可选 "amplify" (扩增) 或 "mutagenesis" (诱变)。
        intent: 一个包含任务特定参数的字典。
            - 对于 "amplify": `{"amplicon": {"start": int, "end": int}, "fwd_overhang": Optional[str], "rev_overhang": Optional[str]}`
            - 对于 "mutagenesis": `{"edits": [{"type": "point", "pos": int, "to": "A|C|G|T"}, ... ]}` (仅支持单个编辑)
    
    Returns:
        一个包含状态信息和设计结果（引物对）的字典。
    """
    # 实现于: tools_pool/design_primer_suite.py
    pass


def simulate_pcr(json_path: str,
                 forward: str,
                 reverse: str,
                 min_anneal_len: int = 15,
                 output_path: Optional[str] = None
                 ) -> Dict[str, Union[List[Fragment], str]]:
    """
    模拟 PCR 扩增。

    实现说明:
    - 模板通过 `json_path` (一个 SequenceRecord 文件的路径) 提供。
    - 引物以原始序列字符串的形式提供，而不是 Primer 对象。
    - 通过 `min_anneal_len` 参数处理带 5' 突出端的引物，该参数指定了用于退火的 3' 末端最小长度。
    - 能正确模拟线性和环状模板上的扩增。
    - 可选地，可以将扩增产物作为一个新的 SequenceRecord JSON 文件保存到 `output_path`。
    - 返回值是一个字典，包含 `amplicons` (作为 `Fragment` 对象的列表) 和一个状态 `message`。
    
    Args:
        json_path: 模板 SequenceRecord JSON 文件的路径。
        forward: 正向引物的序列字符串。
        reverse: 反向引物的序列字符串。
        min_anneal_len: 用于退火的引物 3' 末端最小长度。
        output_path: 可选路径，用于将扩增产物保存为新的 SequenceRecord JSON 文件。
    
    Returns:
        一个包含扩增子列表 (作为 Fragments) 和状态消息的字典。
    """
    # 实现于: tools_pool/simulate_pcr.py
    pass


def simulate_restriction_digestion(json_path: str,
                                   enzyme_names: List[str]) -> DigestResult:
    """
    模拟限制性内切酶协同消化。

    实现说明:
    - 序列通过 `json_path` (一个 SequenceRecord 文件的路径) 提供。
    - 酶通过一个名称列表 (`enzyme_names`) 指定，程序会从内部库中查找酶的规格。不支持自定义酶。
    - 执行所有指定酶的同步消化。
    - 能正确处理线性和环状模板。
    - 会检查切割位点是否离线性模板末端过近，并在结果的 `info` 字段中添加警告。
    - 返回值是代表同步消化结果的单个 `DigestResult` 对象。
    
    Args:
        json_path: 输入的 SequenceRecord JSON 文件路径。
        enzyme_names: 用于消化的酶名称列表 (例如: ["EcoRI", "HindIII"])。
    
    Returns:
        代表同步消化结果的单个 `DigestResult` 对象。
    """
    # 实现于: tools_pool/simulate_restriction_digestion.py
    pass


# ------------------------------
# 未实现的工具 (待开发)
# ------------------------------

def recommend_cloning_strategy(insert: SequenceRecord,
                               vector: SequenceRecord,
                               preferences: Optional[Dict[str, Union[str, int, float, bool]]] = None
                               ) -> AssemblyPlan:
    """
    为插入片段克隆到载体生成推荐策略。
    """
    raise NotImplementedError


def validate_construct(product: SequenceRecord,
                       spec: Dict[str, Union[str, List[str], List[Feature], Dict]]
                       ) -> ValidationReport:
    """
    对构建产物进行验证 (序列、特征、阅读框等)。
    """
    raise NotImplementedError

def synthesize_dna(sequence: str,
                   circular: bool = False,
                   provider_constraints: Optional[Dict[str, Union[int, float, str, bool]]] = None
                   ) -> Fragment:
    """
    DNA 合成 (虚拟)。
    """
    raise NotImplementedError

def simulate_ligation(fragments: List[Fragment],
                      allow_circularization: bool = True,
                      sticky_end_tolerance: bool = False) -> List[SequenceRecord]:
    """
    模拟 T4 连接酶连接 (平末端/粘性末端)。
    """
    raise NotImplementedError

def simulate_homology_based_assembly(fragments: List[Fragment],
                                     overlap_min: int = 20,
                                     method: Literal["gibson", "isothermal", "SLiCE", "CPEC"] = "gibson"
                                     ) -> List[SequenceRecord]:
    """
    模拟基于同源臂的无缝组装 (如 Gibson Assembly)。
    """
    raise NotImplementedError
