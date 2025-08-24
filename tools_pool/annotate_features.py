import sys
import os
import re
from typing import List, Dict, Optional, Tuple, Any
from common_utils.file_operations import load_sequence_from_json, write_record_to_json
from common_utils.sequence import SequenceRecord, Feature

# Add the parent directory to the sys.path
sys.path.append(os.path.join(os.path.dirname(__file__), '..'))

# --- 工具函数 ---
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")

def revcomp(s: str) -> str:
    """Returns the reverse complement of a DNA sequence."""
    return s.translate(_RC)[::-1]

def to_1_based_inclusive_from_fwd_match(m: re.Match) -> Tuple[int, int]:
    """Convert 0-based exclusive match to 1-based inclusive coordinates."""
    start_1b = m.start() + 1
    end_1b = m.end()  # Already 1-based inclusive
    return start_1b, end_1b

def to_1_based_inclusive_from_rev_match(m: re.Match, L: int) -> Tuple[int, int]:
    """Convert reverse complement match to 1-based inclusive coordinates."""
    i, j = m.start(), m.end()
    start_1b = L - j + 1
    end_1b = L - i
    return start_1b, end_1b

def clamp_and_validate(start: int, end: int, L: int) -> Optional[Tuple[int, int]]:
    """Validate and clamp coordinates to sequence bounds."""
    if start > end:
        start, end = end, start  # Swap if needed
    if start < 1 or end > L:
        return None
    return start, end

# --- Common Features Index ---
COMMON_FEATURES = [
    {"name": "T7 Promoter", "type": "promoter", "sequence": "TAATACGACTCACTATAGGG"},
    {"name": "AmpR", "type": "resistance", "sequence": "ATGAGTATTCAACATTTCCGTGTCGCCCTTATTCCCTTTTTTGCGGCATTTTGCCTTCCTGTTTTTGCTCACCCAGAAACGCTGGTGAAAGTAAAAGATGCTGAAGATCAGTTGGGTGCACGAGTGGGTTACATCGAACTGGATCTCAACAGCGGTAAGATCCTTGAGAGTTTTCGCCCCGAAGAACGTTTTCCAATGATGAGCACTTTTAAAGTTCTGCTATGTGGCGCGGTATTATCCCGTATTGACGCCGGGCAAGAGCAACTCGGTCGCCGCATACACTATTCTCAGAATGACTTGGTTGAGTACTCACCAGTCACAGAAAAGCATCTTACGGATGGCATGACAGTAAGAGAATTATGCAGTGCTGCCATAACCATGAGTGATAACACTGCGGCCAACTTACTTCTGACAACGATCGGAGGACCGAAGGAGCTAACCGCTTTTTTGCACAACATGGGGGATCATGTAACTCGCCTTGATCGTTGGGAACCGGAGCTGAATGAAGCCATACCAAACGACGAGCGTGACACCACGATGCCTGTAGCAATGGCAACAACGTTGCGCAAACTATTAACTGGCGAACTACTTACTCTAGCTTCCCGGCAACAATTAATAGACTGGATGGAGGCGGATAAAGTTGCAGGACCACTTCTGCGCTCGGCCCTTCCGGCTGGCTGGTTTATTGCTGATAAATCTGGAGCCGGTGAGCGTGGGTCTCGCGGTATCATTGCAGCACTGGGGCCAGATGGTAAGCCCTCCCGTATCGTAGTTATCTACACGACGGGGAGTCAGGCAACTATGGATGAACGAAATAGACAGATCGCTGAGATAGGTGCCTCACTGATTAAGCATTGG"},
    {"name": "CMV enhancer", "type": "enhancer", "sequence": "gacattgattattgactagttattaatagtaatcaattacggggtcattagttcatagcccatatatggagttccgcgttacataacttacggtaaatggcccgcctggctgaccgcccaacgacccccgcccattgacgtcaataatgacgtatgttcccatagtaacgccaatagggactttccattgacgtcaatgggtggagtatttacggtaaactgcccacttggcagtacatcaagtgtatcatatgccaagtacgccccctattgacgtcaatgacggtaaatggcccgcctggcattatgcccagtacatgaccttatgggactttcctacttggcagtacatctacgtattagtcatcgctattaccatg"},
    {"name": "CMV promoter", "type": "promoter", "sequence": "gtgatgcggttttggcagtacatcaatgggcgtggatagcggtttgactcacggggatttccaagtctccaccccattgacgtcaatgggagtttgttttggcaccaaaatcaacgggactttccaaaatgtcgtaacaactccgccccattgacgcaaatgggcggtaggcgtgtacggtgggaggtctatataagcagagct"},
    {"name": "ori", "type": "origin", "sequence": "tttccataggctccgcccccctgacgagcatcacaaaaatcgacgctcaagtcagaggtggcgaaacccgacaggactataaagataccaggcgtttccccctggaagctccctcgtgcgctctcctgttccgaccctgccgcttaccggatacctgtccgcctttctcccttcgggaagcgtggcgctttctcatagctcacgctgtaggtatctcagttcggtgtaggtcgttcgctccaagctgggctgtgtgcacgaaccccccgttcagcccgaccgctgcgccttatccggtaactatcgtcttgagtccaacccggtaagacacgacttatcgccactggcagcagccactggtaacaggattagcagagcgaggtatgtaggcggtgctacagagttcttgaagtggtggcctaactacggctacactagaagaacagtatttggtatctgcgctctgctgaagccagttaccttcggaaaaagagttggtagctcttgatccggcaaacaaaccaccgctggtagcggtttttttgtttgcaagcagcagattacgcgcagaaaaaaaggatctcaa"},
    {"name": "bGH poly(A) signal", "type": "polyA_signal", "sequence": "ctgtgccttctagttgccagccatctgttgtttgcccctcccccgtgccttccttgaccctggaaggtgccactcccactgtcctttcctaataaaatgaggaaattgcatcgcattgtctgagtaggtgtcattctattctggggggtggggtggggcaggacagcaagggggaggattgggaagacaatagcaggcatgctggggatgcggtgggctctatgg"},
    {"name": "INS CDS", "type": "cds", "sequence": "ATGGCCCTGTGGATGCGCCTCCTGCCCCTGCTGGCGCTGCTGGCCCTCTGGGGACCTGACCCAGCCGCAGCCTTTGTGAACCAACACCTGTGCGGCTCACACCTGGTGGAAGCTCTCTACCTAGTGTGCGGGGAACGAGGCTTCTTCTACACACCCAAGACCCGCCGGGAGGCAGAGGACCTGCAGGTGGGGCAGGTGGAGCTGGGCGGGGGCCCTGGTGCAGGCAGCCTGCAGCCCTTGGCCCTGGAGGGGTCCCTGCAGAAGCGTGGCATTGTGGAACAATGCTGTACCAGCATCTGCTCCCTCTACCAGCTGGAGAACTACTGCAACTAG"},
]

def annotate_features(json_path: str, output_path: Optional[str] = None, record_index: int = 0) -> Dict[str, Any]:
    """
    Adds common features to a SequenceRecord based on sequence matching.
    
    Args:
        json_path (str): Path to the input JSON file containing the SequenceRecord or a list of SequenceRecords.
        output_path (Optional[str]): Path to save the annotated SequenceRecord. 
                                    If not provided, overwrites the input file.
        record_index (int): When the input file contains multiple records, specifies which record to annotate. 
                           Defaults to 0 (the first record).
    
    Returns:
        Dict[str, any]: Dictionary containing:
            - "annotated_file_path": Path to the output JSON file with annotated features.
            - "added_annotations": List of annotations that were successfully added.
    """
    # Load the sequence record
    data = load_sequence_from_json(json_path)
    
    # Handle possible list of records
    if isinstance(data, list):
        if not data:
            raise ValueError("No records found in the input file")
        if record_index >= len(data):
            raise ValueError(f"Record index {record_index} is out of range. The file contains only {len(data)} records.")
        record = data[record_index]
    else:
        record = data
    
    # Get sequence information
    sequence = (record.get("sequence") or "").upper()
    if not sequence:
        raise ValueError("No sequence found in the input file")
    
    L = len(sequence)
    rev_sequence = revcomp(sequence)
    
    # Get existing features
    existing_features = record.get("features", [])
    found_features = list(existing_features)  # Copy existing features
    
    # Create a set to store unique feature identifiers for quick lookups
    existing_feature_signatures = set()
    for f in existing_features:
        label = f.get("qualifiers", {}).get("label", "")
        signature = (label, f.get("start"), f.get("end"), f.get("strand"))
        existing_feature_signatures.add(signature)
        
    # Keep track of added annotations
    added_annotations = []
    
    # Scan for common features
    for feature_info in COMMON_FEATURES:
        name = feature_info["name"]
        feature_type = feature_info["type"]
        pattern = feature_info["sequence"].upper()
        
        # Skip if pattern is empty
        if not pattern:
            continue
            
        # Search on forward strand
        for m in re.finditer(re.escape(pattern), sequence):
            start_1b, end_1b = to_1_based_inclusive_from_fwd_match(m)
            validated_coords = clamp_and_validate(start_1b, end_1b, L)
            if validated_coords:
                start, end = validated_coords
                feature_signature = (name, start, end, 1)
                if feature_signature not in existing_feature_signatures:
                    feature: Feature = {
                        "type": feature_type,
                        "start": int(start),
                        "end": int(end),
                        "strand": 1,
                        "qualifiers": {"label": name}
                    }
                    found_features.append(feature)
                    existing_feature_signatures.add(feature_signature)
                    
                    # Add to added annotations list
                    annotation_desc = f"{name} ({feature_type}) at position {start}-{end} on forward strand"
                    added_annotations.append(annotation_desc)
        
        # Search on reverse strand
        for m in re.finditer(re.escape(pattern), rev_sequence):
            start_1b, end_1b = to_1_based_inclusive_from_rev_match(m, L)
            validated_coords = clamp_and_validate(start_1b, end_1b, L)
            if validated_coords:
                start, end = validated_coords
                feature_signature = (name, start, end, -1)
                if feature_signature not in existing_feature_signatures:
                    feature: Feature = {
                        "type": feature_type,
                        "start": int(start),
                        "end": int(end),
                        "strand": -1,
                        "qualifiers": {"label": name}
                    }
                    found_features.append(feature)
                    existing_feature_signatures.add(feature_signature)
                    
                    # Add to added annotations list
                    annotation_desc = f"{name} ({feature_type}) at position {start}-{end} on reverse strand"
                    added_annotations.append(annotation_desc)
    
    # Update the record with new features
    record["features"] = found_features
    
    # If the input was a list, update the specific record in the list
    if isinstance(data, list):
        data[record_index] = record
        final_data = data
    else:
        final_data = record
    
    # Determine output path
    if output_path is None:
        output_path = json_path
    
    # Save the updated record(s)
    write_record_to_json(final_data, output_path)
    
    return {
        "annotated_file_path": output_path,
        "added_annotations": added_annotations
    }

if __name__ == "__main__":
    # Example usage
    json_path = "data/temp/pcdna3.1(-).json"
    annotated_path = annotate_features(json_path)
    print(f"Features annotated and saved to: {annotated_path}")
