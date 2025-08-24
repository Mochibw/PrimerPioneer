import math
import re
from typing import Dict, List, Tuple, Optional, Union

def calculate_gc_content(sequence: str) -> float:
    """
    计算DNA序列的GC含量。
    
    Args:
        sequence: DNA序列
        
    Returns:
        GC含量（0.0到1.0之间的浮点数）
    """
    seq = sequence.upper()
    gc_count = seq.count('G') + seq.count('C')
    total_bases = len(seq)
    return gc_count / total_bases if total_bases > 0 else 0.0

def calculate_tm(sequence: str, method: str = "wallace", salt_conc: float = 0.05, 
                dna_conc: float = 0.0000005) -> float:
    """
    计算DNA序列的熔解温度(Tm)。
    
    Args:
        sequence: DNA序列
        method: 计算方法（"wallace", "santalucia", 或 "basic"）
        salt_conc: 盐浓度(M)
        dna_conc: DNA浓度(M)
        
    Returns:
        熔解温度（摄氏度）
    """
    seq = sequence.upper()
    length = len(seq)
    
    if method == "wallace":
        # Wallace规则：Tm = 2*(A+T) + 4*(G+C)
        a_count = seq.count('A')
        t_count = seq.count('T')
        g_count = seq.count('G')
        c_count = seq.count('C')
        return 2 * (a_count + t_count) + 4 * (g_count + c_count)
    
    elif method == "santalucia":
        # 简化的SantaLucia算法
        if length < 14:
            return calculate_tm(seq, "wallace")
        else:
            gc_content = calculate_gc_content(seq)
            return 81.5 + 0.41 * (gc_content * 100) - (675 / length) + 16.6 * math.log10(salt_conc)
    
    elif method == "basic":
        # 基本公式考虑盐浓度
        gc_content = calculate_gc_content(seq)
        return 64.9 + 41 * (gc_content - 0.16) + 16.6 * math.log10(salt_conc)
    
    else:
        raise ValueError(f"不支持的Tm计算方法: {method}")

def reverse_complement(sequence: str) -> str:
    """
    获取DNA序列的反向互补序列。
    
    Args:
        sequence: DNA序列
        
    Returns:
        反向互补序列
    """
    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 
                 'N': 'N', 'R': 'Y', 'Y': 'R', 'S': 'S', 
                 'W': 'W', 'K': 'M', 'M': 'K', 'B': 'V', 
                 'D': 'H', 'H': 'D', 'V': 'B'}
    
    return ''.join(complement.get(base, 'N') for base in reversed(sequence.upper()))

def check_pcr_feasibility(sequence: str, 
                         min_tm: float = 55.0, 
                         max_tm: float = 65.0,
                         max_gc: float = 0.7,
                         min_gc: float = 0.3,
                         max_length: int = 3000) -> Dict:
    """
    检查序列是否适合PCR扩增。
    
    Args:
        sequence: 要检查的序列
        min_tm: 最小Tm值
        max_tm: 最大Tm值
        max_gc: 最大GC含量
        min_gc: 最小GC含量
        max_length: 最大长度限制
        
    Returns:
        包含可行性检查结果的字典
    """
    gc_content = calculate_gc_content(sequence)
    tm = calculate_tm(sequence)
    length = len(sequence)
    
    issues = []
    
    if length > max_length:
        issues.append(f"序列过长: {length} bp > {max_length} bp")
    
    if gc_content > max_gc:
        issues.append(f"GC含量过高: {gc_content:.1%} > {max_gc:.0%}")
    
    if gc_content < min_gc:
        issues.append(f"GC含量过低: {gc_content:.1%} < {min_gc:.0%}")
    
    if tm > max_tm:
        issues.append(f"Tm值过高: {tm:.1f}°C > {max_tm:.1f}°C")
    
    if tm < min_tm:
        issues.append(f"Tm值过低: {tm:.1f}°C < {min_tm:.1f}°C")
    
    # 检查重复序列
    if has_high_repetitiveness(sequence):
        issues.append("序列包含高度重复区域")
    
    # 检查同聚物
    homopolymer_info = check_homopolymers(sequence, max_length=6)
    if homopolymer_info["has_long_homopolymer"]:
        issues.append(f"检测到过长同聚物: {homopolymer_info['longest_homopolymer']}个碱基")
    
    # 检查二级结构
    secondary_structure = predict_secondary_structure(sequence)
    if secondary_structure["has_strong_secondary_structure"]:
        issues.append("序列可能形成强二级结构")
    
    feasible = len(issues) == 0
    
    recommendation = ""
    if not feasible:
        if gc_content > 0.7 or has_high_repetitiveness(sequence) or length > 2000:
            recommendation = "建议使用DNA化学合成而非PCR"
        else:
            recommendation = "尝试优化PCR条件（如使用GC buffer、添加剂）或重新设计引物"
    
    return {
        "feasible": feasible,
        "length": length,
        "gc_content": gc_content,
        "tm": tm,
        "issues": issues,
        "recommendation": recommendation
    }

def has_high_repetitiveness(sequence: str, threshold: float = 0.6) -> bool:
    """
    检查序列是否高度重复。
    
    Args:
        sequence: 要检查的序列
        threshold: 重复度阈值
        
    Returns:
        是否高度重复
    """
    if len(sequence) < 20:
        return False
    
    # 使用序列熵来评估重复性
    base_counts = {}
    for base in sequence:
        base_counts[base] = base_counts.get(base, 0) + 1
    
    # 计算熵
    entropy = 0
    total = len(sequence)
    for count in base_counts.values():
        probability = count / total
        if probability > 0:
            entropy -= probability * math.log2(probability)
    
    # 最大熵（当所有碱基均匀分布时）
    max_entropy = math.log2(min(4, len(set(sequence))))
    
    # 低熵表示高重复性
    return (entropy / max_entropy) < threshold if max_entropy > 0 else False

def check_homopolymers(sequence: str, max_length: int = 6) -> Dict:
    """
    检查序列中的同聚物。
    
    Args:
        sequence: 要检查的序列
        max_length: 最大允许的同聚物长度
        
    Returns:
        同聚物检查结果
    """
    current_base = None
    current_length = 0
    max_homopolymer = 0
    homopolymer_positions = []
    
    for i, base in enumerate(sequence):
        if base == current_base:
            current_length += 1
            if current_length > max_homopolymer:
                max_homopolymer = current_length
        else:
            if current_length > max_length:
                homopolymer_positions.append({
                    "base": current_base,
                    "length": current_length,
                    "position": i - current_length
                })
            current_base = base
            current_length = 1
    
    # 检查最后一个同聚物
    if current_length > max_length:
        homopolymer_positions.append({
            "base": current_base,
            "length": current_length,
            "position": len(sequence) - current_length
        })
    
    return {
        "has_long_homopolymer": max_homopolymer > max_length,
        "longest_homopolymer": max_homopolymer,
        "problematic_homopolymers": homopolymer_positions
    }

def predict_secondary_structure(sequence: str, min_stem_length: int = 4) -> Dict:
    """
    预测序列的二级结构形成倾向。
    
    Args:
        sequence: DNA序列
        min_stem_length: 最小茎长度
        
    Returns:
        二级结构预测结果
    """
    seq = sequence.upper()
    has_strong_secondary_structure = False
    potential_structures = []
    
    # 简单检查回文序列（可能形成发夹结构）
    for stem_length in range(min_stem_length, min(10, len(seq)//2)):
        for i in range(len(seq) - 2 * stem_length):
            stem = seq[i:i+stem_length]
            loop_start = i + stem_length
            loop_end = loop_start + 2  # 假设最小环大小为2
            complementary = reverse_complement(seq[loop_end:loop_end+stem_length])
            
            if stem == complementary:
                has_strong_secondary_structure = True
                potential_structures.append({
                    "type": "hairpin",
                    "stem_length": stem_length,
                    "position": i,
                    "sequence": seq[i:loop_end+stem_length]
                })
    
    return {
        "has_strong_secondary_structure": has_strong_secondary_structure,
        "potential_structures": potential_structures
    }

def find_restriction_sites(sequence: str, enzyme_list: Optional[List[str]] = None) -> Dict:
    """
    在序列中查找限制性酶切位点。
    
    Args:
        sequence: DNA序列
        enzyme_list: 要查找的酶列表（如果为None，则使用常见酶）
        
    Returns:
        酶切位点信息
    """
    common_enzymes = {
        "EcoRI": "GAATTC",
        "BamHI": "GGATCC",
        "HindIII": "AAGCTT",
        "XhoI": "CTCGAG",
        "NotI": "GCGGCCGC",
        "KpnI": "GGTACC",
        "SacI": "GAGCTC",
        "SalI": "GTCGAC",
        "PstI": "CTGCAG",
        "NcoI": "CCATGG"
    }
    
    enzymes_to_check = enzyme_list or list(common_enzymes.keys())
    restriction_sites = {}
    
    seq = sequence.upper()
    
    for enzyme in enzymes_to_check:
        site = common_enzymes.get(enzyme, "")
        if not site:
            continue
        
        positions = []
        for match in re.finditer(site, seq):
            positions.append(match.start() + 1)  # 1-based position
        
        if positions:
            restriction_sites[enzyme] = {
                "recognition_site": site,
                "positions": positions,
                "count": len(positions)
            }
    
    return restriction_sites

def calculate_molecular_weight(sequence: str, strandedness: str = "double") -> float:
    """
    计算DNA序列的分子量。
    
    Args:
        sequence: DNA序列
        strandedness: "single" 或 "double"
        
    Returns:
        分子量（道尔顿）
    """
    # 碱基分子量（单链）
    base_weights = {
        'A': 313.2, 'T': 304.2, 'C': 289.2, 'G': 329.2,
        'N': 300.0  # 平均值
    }
    
    seq = sequence.upper()
    total_weight = 0.0
    
    for base in seq:
        total_weight += base_weights.get(base, 300.0)
    
    # 减去磷酸二酯键形成时失去的水分子
    total_weight -= (len(seq) - 1) * 18.0
    
    # 双链DNA
    if strandedness == "double":
        total_weight *= 2
        # 减去氢键（近似）
        total_weight -= len(seq) * 2.0
    
    return total_weight

def format_sequence(sequence: str, line_length: int = 80) -> str:
    """
    格式化DNA序列，添加行号和换行。
    
    Args:
        sequence: DNA序列
        line_length: 每行碱基数
        
    Returns:
        格式化后的序列
    """
    formatted = []
    for i in range(0, len(sequence), line_length):
        line = sequence[i:i+line_length]
        line_number = i + 1
        formatted.append(f"{line_number:6d} {line}")
    
    return "\n".join(formatted)

def validate_dna_sequence(sequence: str) -> Dict:
    """
    验证DNA序列的有效性。
    
    Args:
        sequence: 要验证的序列
        
    Returns:
        验证结果
    """
    seq = sequence.upper()
    valid_bases = set('ATCGNRYSWKMBDHV')
    
    invalid_chars = []
    for i, char in enumerate(seq):
        if char not in valid_bases:
            invalid_chars.append({
                "position": i + 1,
                "character": char
            })
    
    is_valid = len(invalid_chars) == 0
    
    return {
        "is_valid": is_valid,
        "length": len(seq),
        "invalid_characters": invalid_chars,
        "gc_content": calculate_gc_content(seq) if is_valid else None
    }
