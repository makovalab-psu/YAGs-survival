#!/usr/bin/env python3

import json
import sys
import os

def main():
    # Check if enough arguments were provided
    if len(sys.argv) < 3:
        print("Error: Insufficient arguments.")
        print(f"Usage: python {os.path.basename(sys.argv[0])} output_file.txt file1.json file2.json [file3.json ...]")
        return 1
    
    # First argument is the output file
    output_file_path = sys.argv[1]
    
    # Remaining arguments are JSON files to process
    file_paths = sys.argv[2:]
    
    # Dictionary to store loaded JSON data
    data_dict = {}
    
    # Load all JSON files
    for file_path in file_paths:
        try:
            with open(file_path, 'r') as file:
                # Use the filename (without extension) as the key
                filename = os.path.splitext(os.path.basename(file_path))[0]
                data_dict[filename] = json.load(file)
                print(f"Loaded {file_path}")
        except FileNotFoundError:
            print(f"Warning: File {file_path} not found. Skipping.")
        except json.JSONDecodeError:
            print(f"Warning: File {file_path} is not valid JSON. Skipping.")
    
    if not data_dict:
        print("Error: No valid JSON files were loaded.")
        return 1
    
    # Collect all keys from all files
    all_keys = set()
    for data in data_dict.values():
        all_keys.update(data.keys())
    
    try:
        # Open output file for writing
        with open(output_file_path, 'w') as output_file:
            # Write header
            header = "key"
            for filename in data_dict.keys():
                header += f"\t{filename}"
            header += "\tsum"
            output_file.write(header + "\n")
            
            # Write data rows
            for key in sorted(all_keys):
                row = key
                sum_value = 0
                
                for filename, data in data_dict.items():
                    value = data.get(key, 0)
                    row += f"\t{value}"
                    sum_value += value
                
                row += f"\t{sum_value}"
                output_file.write(row + "\n")
        
        print(f"Output successfully written to {output_file_path}")
    
    except IOError as e:
        print(f"Error writing to output file: {e}")
        return 1
    
    return 0

if __name__ == "__main__":
    sys.exit(main())