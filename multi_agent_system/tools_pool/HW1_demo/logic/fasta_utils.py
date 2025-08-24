from pathlib import Path

def read_fasta(file_path: str) -> str:
    """
    Reads a FASTA file and returns the concatenated sequence.
    Tries to read with UTF-8, then GBK, then latin-1 encoding.
    Ignores header lines (lines starting with '>').
    """
    sequence = ""
    encodings_to_try = ["utf-8", "gbk", "latin-1"]
    for encoding in encodings_to_try:
        try:
            with Path(file_path).open("r", encoding=encoding) as f:
                for line in f:
                    line = line.strip()
                    if not line.startswith(">") and line:
                        sequence += line
            return sequence.upper()  # Return on success
        except UnicodeDecodeError:
            continue  # Try the next encoding
        except FileNotFoundError:
            raise FileNotFoundError(f"FASTA file not found at: {file_path}")
        except Exception as e:
            raise IOError(f"Error reading FASTA file {file_path}: {e}")
    
    # If all encodings fail
    raise IOError(f"Could not decode file {file_path} with any of the attempted encodings: {encodings_to_try}")
