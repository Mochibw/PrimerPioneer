import sys
import os
import json
import uuid
from typing import List, Dict, Optional, Union
from langchain.tools import tool

# Add the parent directory to the sys.path to find common_utils
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

from common_utils.sequence import Fragment, SequenceRecord
from common_utils.file_operations import load_sequence_from_json, write_record_to_json

# ---- Helper Functions ----
_RC_MAP = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def _revcomp(s: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    return s.translate(_RC_MAP)[::-1]

@tool()
def simulate_pcr(json_path: str,
                 forward: str,
                 reverse: str,
                 min_anneal_len: int = 15,
                 output_path: Optional[str] = None
                 ) -> Dict[str, Union[List[Fragment], str]]:
    """
    Simulates PCR, saves the resulting amplicon as a SequenceRecord, and returns the result.

    Allows for 5' overhangs on primers by matching the 3' annealing region.
    If an output_path is provided, the generated amplicon is written to a JSON file as a SequenceRecord.

    Args:
        json_path (str): Path to the JSON file with the template SequenceRecord.
        forward (str): The forward primer sequence.
        reverse (str): The reverse primer sequence.
        output_path (Optional[str]): If provided, the path to save the generated amplicon
                                     as a SequenceRecord JSON file.

    Returns:
        Dict[str, Union[List[Fragment], str]]: A dictionary containing the list of
        amplicons (as Fragment objects) and a confirmation message.

    Raises:
        ValueError: If primers are shorter than `min_anneal_len`.
    """

    template_record = load_sequence_from_json(json_path)
    template_seq = (template_record.get("sequence") or "").upper()
    is_circular = template_record.get("circular", False)
    L = len(template_seq)

    if not template_seq or not forward or not reverse:
        raise ValueError("Template and primer sequences must not be empty.")

    forward_primer = forward.upper()
    reverse_primer = reverse.upper()

    if len(forward_primer) < min_anneal_len or len(reverse_primer) < min_anneal_len:
        raise ValueError(f"Primers must be at least {min_anneal_len} bases long.")

    # Get the 3' annealing part of the primers
    forward_anneal = forward_primer[-min_anneal_len:]
    reverse_anneal_rc = _revcomp(reverse_primer[-min_anneal_len:])

    # --- Debugging prints ---
    print(f"Template length: {len(template_seq)}")
    print(f"Forward annealing sequence ({len(forward_anneal)}bp): {forward_anneal}")
    print(f"Reverse annealing sequence (RC, {len(reverse_anneal_rc)}bp): {reverse_anneal_rc}")
    # --- End Debugging prints ---

    # Find binding sites of the annealing parts
    fwd_pos = template_seq.find(forward_anneal)
    rev_pos = template_seq.find(reverse_anneal_rc)

    # --- Debugging prints ---
    print(f"Forward position found: {fwd_pos}")
    print(f"Reverse position found: {rev_pos}")
    # --- End Debugging prints ---

    if fwd_pos == -1 or rev_pos == -1:
        message = "PCR simulation completed. No amplicon produced because one or both primers did not bind."
        if output_path:
            message += " Nothing written to file."
        return {"amplicons": [], "message": message}

    # Template coordinates are based on the annealing part
    fwd_template_start = fwd_pos
    rev_template_end = rev_pos + len(reverse_anneal_rc)

    # 1-based coordinates for the Fragment object
    start_coord = fwd_template_start + 1
    end_coord = rev_template_end

    amplicon_seq = ""
    amplicon_len = 0

    if fwd_template_start < rev_template_end:  # Standard linear or intra-circular
        middle_template_start = fwd_template_start + min_anneal_len
        middle_template_end = rev_pos
        middle_template = template_seq[middle_template_start:middle_template_end]
        
        amplicon_seq = forward_primer + middle_template + _revcomp(reverse_primer)
        amplicon_len = len(amplicon_seq)
    elif is_circular:  # Amplification across the circular origin
        middle_template = template_seq[fwd_pos + min_anneal_len:] + template_seq[:rev_pos]
        
        amplicon_seq = forward_primer + middle_template + _revcomp(reverse_primer)
        amplicon_len = len(amplicon_seq)
    else:
        return {"amplicons": [], "message": "No PCR product formed due to primer binding issues on linear template."}

    amplicon: Fragment = {
        "id": str(uuid.uuid4()),
        "start": start_coord,
        "end": end_coord,
        "length": amplicon_len,
        "strand": 1,
        "sequence": amplicon_seq,
        "overhang_5": {"kind": "blunt", "seq": "", "length": 0},
        "overhang_3": {"kind": "blunt", "seq": "", "length": 0},
    }

    # Create a SequenceRecord for the amplicon to be saved
    amplicon_record: SequenceRecord = {
        "id": amplicon["id"],
        "name": f"PCR Product from {template_record.get('name', 'Unknown Template')}",
        "sequence": amplicon["sequence"],
        "length": amplicon["length"],
        "circular": False, # PCR products are typically linear
        "features": [] # Can add features like primer binding sites here if needed
    }

    amplicons_list_for_return = [amplicon] # This list is for the function's return value

    message = "PCR simulation completed successfully. If you add restriction sites on your primers, the next movement should be simulate_restriction_digestion."

    if output_path:
        abs_path = os.path.abspath(output_path)
        try:
            write_record_to_json(amplicon_record, output_path) # Save the SequenceRecord
            message = f"PCR simulation completed. Amplicon written as SequenceRecord to {abs_path}"
        except Exception as e:
            message = f"PCR simulation completed, but failed to write to {abs_path}. Error: {e}"

    return {"amplicons": amplicons_list_for_return, "message": message}


if __name__ == "__main__":
    # Debug parameters
    forward = "GATCGCTAGCATGGAGGAGGCTGAGCTGGA"
    json_path = "data\\temp\\RTCB_Homo_sapiens_lcl_NM_014306.5_cds_NP_055121.1_1.json"
    min_anneal_len = 18
    output_path = "data\\temp\\RTCB_amplicon.json"
    reverse = "GATCCTCGAGTCCCTCAGAAGATGAGCTGA"
    
    # Run PCR simulation
    result = simulate_pcr(json_path, forward, reverse, min_anneal_len, output_path)
    print(result["message"])
    
    # Print amplicon details if any were produced
    if result["amplicons"]:
        for amplicon in result["amplicons"]:
            print(f"Amplicon ID: {amplicon['id']}")
            print(f"Length: {amplicon['length']}")
            print(f"Sequence: {amplicon['sequence'][:50]}...")  # Show first 50 bases
