

def modify_file_by_lines(file_path, new_file_path, modification_dict):

    try:
            # Read the file into a list of lines
            with open(file_path, 'r') as file:
                lines = file.readlines()
            
            # Apply modifications
            for line_number, new_content in modification_dict.items():
                if 1 <= line_number <= len(lines):
                    lines[line_number - 1] = new_content + '\n'
                else:
                    print(f"Line {line_number} is out of range. Skipping.")
            
            # Write the modified lines back to the file
            with open(new_file_path, 'w') as file:
                file.writelines(lines)

            print("File modified successfully.")
        
    except FileNotFoundError:
        print(f"Error: The file at '{file_path}' was not found.")
    except Exception as e:
        print(f"An unexpected error occurred: {e}")