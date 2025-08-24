from Bio.Seq import Seq
from Bio.SeqUtils import MeltingTemp as mt
from Bio.SeqUtils import gc_fraction # For GC content calculation

# Define target parameters for primer design
MIN_PRIMER_LEN = 18
MAX_PRIMER_LEN = 25
TARGET_GC_MIN = 40.0
TARGET_GC_MAX = 60.0
TARGET_TM_MIN = 58.0
TARGET_TM_MAX = 65.0

# Default concentrations for Tm calculation (in mM and nM)
DEFAULT_NA_CONC = 50.0

# Protective bases to add to the 5' end of the primer, before the restriction site
PROTECTIVE_BASES = "AATT"

def _calculate_tm(sequence: str, Na: float = DEFAULT_NA_CONC) -> float:
    """
    Calculates the melting temperature (Tm) of a DNA sequence using the nearest-neighbor method.
    """
    # Tm_NN requires sequence to be a Biopython Seq object
    seq_obj = Seq(sequence)
    return mt.Tm_NN(seq_obj, Na=Na)

def _find_optimal_primer(
    target_sequence: str,
    is_forward: bool,
    enzyme_site: str
) -> dict:
    """
    Iteratively searches for an optimal primer sequence within the target_sequence.
    Considers length, GC content, and Tm for the binding part.
    Returns metrics for both binding part and full primer.
    """
    best_primer_info = None
    min_deviation = float('inf')

    # Iterate through possible primer lengths for the binding part
    for binding_length in range(MIN_PRIMER_LEN, MAX_PRIMER_LEN + 1):
        if is_forward:
            # Forward primer: take from the beginning of the target_sequence
            if len(target_sequence) < binding_length:
                continue # Not enough sequence for this length
            primer_binding_part = target_sequence[:binding_length]
        else:
            # Reverse primer: take from the end, then reverse complement
            if len(target_sequence) < binding_length:
                continue # Not enough sequence for this length
            primer_binding_part = Seq(target_sequence[-binding_length:]).reverse_complement()

        # Construct the full primer sequence (Protective Bases + Enzyme Site + Binding Part)
        full_primer_seq = PROTECTIVE_BASES + enzyme_site + str(primer_binding_part)

        # Calculate metrics for the binding part
        binding_gc = gc_fraction(str(primer_binding_part)) * 100
        binding_tm = _calculate_tm(str(primer_binding_part))

        # Calculate metrics for the full primer
        full_gc = gc_fraction(full_primer_seq) * 100
        full_tm = _calculate_tm(full_primer_seq)

        # Check if criteria are met (based on binding part)
        is_len_ok = MIN_PRIMER_LEN <= len(primer_binding_part) <= MAX_PRIMER_LEN
        is_gc_ok = TARGET_GC_MIN <= binding_gc <= TARGET_GC_MAX
        is_tm_ok = TARGET_TM_MIN <= binding_tm <= TARGET_TM_MAX

        if is_len_ok and is_gc_ok and is_tm_ok:
            # Found a perfect primer, return it immediately
            return {
                "binding_part_sequence": str(primer_binding_part),
                "full_primer_sequence": full_primer_seq,
                "binding_part_length": len(primer_binding_part),
                "full_primer_length": len(full_primer_seq),
                "binding_part_gc_content": round(binding_gc, 2),
                "binding_part_tm": round(binding_tm, 2),
                "full_primer_gc_content": round(full_gc, 2),
                "full_primer_tm": round(full_tm, 2),
                "notes": "Optimal primer found within specified criteria."
            }
        else:
            # Calculate deviation for "best effort" if no perfect primer is found
            deviation = 0
            if not is_len_ok:
                deviation += min(abs(len(primer_binding_part) - MIN_PRIMER_LEN),
                                 abs(len(primer_binding_part) - MAX_PRIMER_LEN)) * 2 # Higher penalty for length
            if not is_gc_ok:
                deviation += min(abs(binding_gc - TARGET_GC_MIN), abs(binding_gc - TARGET_GC_MAX))
            if not is_tm_ok:
                deviation += min(abs(binding_tm - TARGET_TM_MIN), abs(binding_tm - TARGET_TM_MAX))

            if deviation < min_deviation:
                min_deviation = deviation
                notes_list = ["Best effort primer found, but not all criteria met:"]
                if not is_len_ok:
                    notes_list.append(f"  - Length ({len(primer_binding_part)} bp) not in [{MIN_PRIMER_LEN}-{MAX_PRIMER_LEN}] bp.")
                if not is_gc_ok:
                    notes_list.append(f"  - GC Content ({round(binding_gc, 2)}%) not in [{TARGET_GC_MIN}-{TARGET_GC_MAX}]%.")
                if not is_tm_ok:
                    notes_list.append(f"  - Tm ({round(binding_tm, 2)}°C) not in [{TARGET_TM_MIN}-{TARGET_TM_MAX}]°C.")
                
                best_primer_info = {
                    "binding_part_sequence": str(primer_binding_part),
                    "full_primer_sequence": full_primer_seq,
                    "binding_part_length": len(primer_binding_part),
                    "full_primer_length": len(full_primer_seq),
                    "binding_part_gc_content": round(binding_gc, 2),
                    "binding_part_tm": round(binding_tm, 2),
                    "full_primer_gc_content": round(full_gc, 2),
                    "full_primer_tm": round(full_tm, 2),
                    "notes": "\n".join(notes_list)
                }
    return best_primer_info if best_primer_info else {
        "binding_part_sequence": "", "full_primer_sequence": "",
        "binding_part_length": 0, "full_primer_length": 0,
        "binding_part_gc_content": 0.0, "binding_part_tm": 0.0,
        "full_primer_gc_content": 0.0, "full_primer_tm": 0.0,
        "notes": "Could not find any suitable primer within the given constraints."
    }


def design_primers_logic(cds_sequence: str, forward_enzyme_site: str, reverse_enzyme_site: str) -> dict:
    """
    Designs forward and reverse primers based on CDS sequence and restriction enzyme sites.
    Iteratively optimizes for length, GC content, and Tm.
    """
    if not cds_sequence:
        raise ValueError("Invalid CDS sequence: The provided sequence is empty.")

    # Normalize enzyme sites (ensure uppercase)
    forward_enzyme_site = forward_enzyme_site.upper()
    reverse_enzyme_site = reverse_enzyme_site.upper()

    # Design Forward Primer
    forward_primer_data = _find_optimal_primer(cds_sequence, True, forward_enzyme_site)

    # Design Reverse Primer
    # For reverse primer, the target sequence for binding is the reverse complement of the CDS end
    reverse_primer_data = _find_optimal_primer(cds_sequence, False, reverse_enzyme_site)

    return {
        "forward_primer": f"5'-{forward_primer_data['full_primer_sequence']}-3'",
        "forward_primer_details": {
            "binding_part_sequence": forward_primer_data['binding_part_sequence'],
            "binding_part_length": forward_primer_data['binding_part_length'],
            "binding_part_gc_content": forward_primer_data['binding_part_gc_content'],
            "binding_part_tm": forward_primer_data['binding_part_tm'],
            "full_primer_length": forward_primer_data['full_primer_length'],
            "full_primer_gc_content": forward_primer_data['full_primer_gc_content'],
            "full_primer_tm": forward_primer_data['full_primer_tm'],
            "notes": forward_primer_data['notes']
        },
        "reverse_primer": f"5'-{reverse_primer_data['full_primer_sequence']}-3'",
        "reverse_primer_details": {
            "binding_part_sequence": reverse_primer_data['binding_part_sequence'],
            "binding_part_length": reverse_primer_data['binding_part_length'],
            "binding_part_gc_content": reverse_primer_data['binding_part_gc_content'],
            "binding_part_tm": reverse_primer_data['binding_part_tm'],
            "full_primer_length": reverse_primer_data['full_primer_length'],
            "full_primer_gc_content": reverse_primer_data['full_primer_gc_content'],
            "full_primer_tm": reverse_primer_data['full_primer_tm'],
            "notes": reverse_primer_data['notes']
        },
        "overall_notes": "Primers designed considering length, GC content, and Tm. Tm calculated using Nearest-Neighbor method. Metrics provided for both binding part and full primer."
    }
