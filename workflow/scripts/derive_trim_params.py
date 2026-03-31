import json
import zipfile
import re
import os
import sys
from pathlib import Path

def parse_fastqc_data(zip_path):
    """Extracts key modules from FastQC zip."""
    data = {}
    with zipfile.ZipFile(zip_path, 'r') as z:
        # FastQC names the folder inside the zip the same as the zip filename without .zip
        folder_name = Path(zip_path).stem
        data_file = f"{folder_name}/fastqc_data.txt"
        
        with z.open(data_file) as f:
            content = f.read().decode('utf-8')
            
    # Extract modules
    # 1. Basic Statistics (for length)
    m = re.search(r">>Basic Statistics\t(.*?)\n(.*?)\n>>END_MODULE", content, re.DOTALL)
    if m:
        stats = m.group(2)
        match = re.search(r"Sequence length\t(\d+)-?(\d+)?", stats)
        if match:
            # Use the upper bound of length or the single value
            data['max_length'] = int(match.group(2)) if match.group(2) else int(match.group(1))
            data['min_length'] = int(match.group(1))

    # 2. Per base sequence quality (for decay)
    m = re.search(r">>Per base sequence quality\t(.*?)\n(.*?)\n>>END_MODULE", content, re.DOTALL)
    if m:
        qual_lines = m.group(2).strip().split('\n')[1:] # Skip header
        quals = []
        for line in qual_lines:
            parts = line.split('\t')
            # Base, Mean, Median, Lower Quartile, Upper Quartile, 10th Percentile, 90th Percentile
            quals.append({
                'base': parts[0],
                'median': float(parts[2]),
                'q1': float(parts[3]),
                'p10': float(parts[5])
            })
        data['quality_profile'] = quals

    # 3. Adapter Content
    m = re.search(r">>Adapter Content\t(.*?)\n(.*?)\n>>END_MODULE", content, re.DOTALL)
    if m:
        adapter_lines = m.group(2).strip().split('\n')[1:] # Skip header
        adapters = {}
        for line in adapter_lines:
            parts = line.split('\t')
            # Position, Illumina Universal, Nextera, Small RNA, etc.
            # We care about the max value across all positions
            for i, adapter_name in enumerate(["Illumina Universal", "Nextera", "Small RNA"]):
                if i+1 < len(parts):
                    val = float(parts[i+1])
                    adapters[adapter_name] = max(adapters.get(adapter_name, 0), val)
        data['adapters'] = adapters

    return data

def derive_params(fastqc_data_list, is_paired, extra_params=""):
    """Derives Trimmomatic parameters based on FastQC results."""
    
    # 1. Adapter Selection
    # Combine adapter info from all mates
    adapter_scores = {"Illumina Universal": 0, "Nextera": 0, "Small RNA": 0}
    for data in fastqc_data_list:
        if 'adapters' in data:
            for k, v in data['adapters'].items():
                adapter_scores[k] = max(adapter_scores[k], v)
    
    best_adapter = max(adapter_scores, key=adapter_scores.get)
    adapter_file = "TruSeq3" # Default
    if adapter_scores[best_adapter] > 0.1: # Threshold to care about adapters
        if "Nextera" in best_adapter:
            adapter_file = "NexteraPE"
        elif "Small RNA" in best_adapter:
            adapter_file = "TruSeq3" # Trimmomatic TruSeq3 handles small RNA well enough or has specific files
        else:
            adapter_file = "TruSeq3"
    
    suffix = "PE" if is_paired else "SE"
    # Note: Nextera only has -PE.fa in standard installations
    if adapter_file == "NexteraPE":
        adapter_path = f"$CONDA_PREFIX/share/trimmomatic/adapters/{adapter_file}-PE.fa"
    else:
        adapter_path = f"$CONDA_PREFIX/share/trimmomatic/adapters/{adapter_file}-{suffix}.fa"
    
    illuminaclip = f"ILLUMINACLIP:{adapter_path}:2:30:10"
    
    # 2. Quality Trimming
    # We'll use MAXINFO if available, or SLIDINGWINDOW
    # Standard choice: SLIDINGWINDOW:4:20
    # If the library is very high quality (p10 > 25 everywhere), maybe 4:25
    # If it's low quality early on, maybe more aggressive.
    
    overall_p10 = []
    for data in fastqc_data_list:
        if 'quality_profile' in data:
            overall_p10.extend([q['p10'] for q in data['quality_profile']])
    
    avg_p10 = sum(overall_p10) / len(overall_p10) if overall_p10 else 20
    
    # Conservative default
    quality_step = "SLIDINGWINDOW:4:20"
    if avg_p10 < 15:
        quality_step = "MAXINFO:40:0.5" # More aggressive balancing
    
    # 3. Min Length
    lengths = [data.get('max_length', 50) for data in fastqc_data_list]
    max_len = max(lengths)
    min_len_val = max(36, int(max_len * 0.35)) # At least 36 or 35% of max length
    minlen_step = f"MINLEN:{min_len_val}"
    
    # 4. Leading / Trailing
    base_steps = ["LEADING:3", "TRAILING:3"]
    
    # 5. Merge and Assemble
    # If user provided extra_params, we respect them. 
    # We replace our derived steps if the user provided the same "KEY:".
    
    derived_steps = [illuminaclip, quality_step] + base_steps + [minlen_step]
    
    if extra_params:
        user_modules = extra_params.split()
        user_module_types = [m.split(':')[0] for m in user_modules]
        
        final_steps = []
        # Add user modules first
        final_steps.extend(user_modules)
        
        # Add derived modules only if user didn't provide them
        for d in derived_steps:
            d_type = d.split(':')[0]
            if d_type not in user_module_types:
                final_steps.append(d)
        
        command_string = " ".join(final_steps)
    else:
        command_string = " ".join(derived_steps)
        
    return command_string

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("--fastqc", nargs="+", required=True, help="FastQC zip files")
    parser.add_argument("--is-paired", action="store_true", help="Paired-end mode")
    parser.add_argument("--extra-params", default="", help="User provided fixed parameters")
    parser.add_argument("--output", required=True, help="Output JSON file")
    args = parser.parse_args()
    
    fastqc_results = []
    for f in args.fastqc:
        fastqc_results.append(parse_fastqc_data(f))
        
    command_string = derive_params(fastqc_results, args.is_paired, args.extra_params)
    
    with open(args.output, 'w') as f:
        json.dump({"trimmomatic_params": command_string}, f)

if __name__ == "__main__":
    main()
