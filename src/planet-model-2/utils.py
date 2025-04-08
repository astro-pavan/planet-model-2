import re
import pandas as pd

def modify_file_by_lines(filename, new_filename, modification_dict):

    try:
            # Read the file into a list of lines
            with open(filename, 'r') as file:
                lines = file.readlines()
            
            # Apply modifications
            for line_number, new_content in modification_dict.items():
                if 1 <= line_number <= len(lines):
                    lines[line_number - 1] = new_content + '\n'
                else:
                    print(f"Line {line_number} is out of range. Skipping.")
            
            # Write the modified lines back to the file
            with open(new_filename, 'w') as file:
                file.writelines(lines)

            print("File modified successfully.")
        
    except FileNotFoundError:
        print(f"Error: The file at '{filename}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")


def extract_block_lines(filename, header_keyword, block_number=1):
    """Extracts lines under a dashed header block."""
    with open(filename, 'r') as file:
        lines = file.readlines()

    header_pattern = re.compile(rf"-+{re.escape(header_keyword)}-+")
    header_indices = [i for i, line in enumerate(lines) if header_pattern.search(line)]

    if len(header_indices) < block_number or block_number < 1:
        raise ValueError(f"Block number {block_number} for header '{header_keyword}' not found.")

    start_idx = header_indices[block_number - 1] + 1
    block_lines = []

    for line in lines[start_idx:]:
        if header_pattern.search(line) or (line.strip() == ""):
            break
        block_lines.append(line.rstrip('\n'))

    return block_lines


def parse_key_value_block(filename, header_keyword, block_number=1):
    """Parses a key-value-style block into a dictionary."""
    block_lines = extract_block_lines(filename, header_keyword, block_number)

    pattern = re.compile(r"(?P<key>.+?)\s*=\s*(?P<value>[-+eE0-9.]+)")
    data = {}
    for line in block_lines:
        match = pattern.search(line)
        if match:
            key = match.group("key").strip()
            value = float(match.group("value"))
            data[key] = value
    return data


def parse_table_block(filename, header_keyword, block_number=1):
    """Parses a fixed-width or space-delimited table block into a DataFrame."""
    block_lines = extract_block_lines(filename, header_keyword, block_number)

    # Skip blank lines
    block_lines = [line for line in block_lines if line.strip()]

    if len(block_lines) < 2:
        raise ValueError("Table block must contain at least header and one data line.")

    # Try to detect header line(s)
    header_line = block_lines[0]
    if re.search(r"\s{2,}", block_lines[1]):
        # Multi-line header? Merge first two lines
        header_line += " " + block_lines[1]
        data_start_idx = 2
    else:
        data_start_idx = 1

    # Parse using pandas
    from io import StringIO
    table_text = "\n".join([header_line] + block_lines[data_start_idx:])
    df = pd.read_csv(StringIO(table_text), delim_whitespace=True, engine="python", comment='#', skip_blank_lines=True)

    return df