import sys

selected_chromosomes = set(snakemake.params.selected_chromosomes)

name_map_in = next(f for f in snakemake.input if f.endswith("-tRNAs_name_map.txt"))
name_map_out = next(f for f in snakemake.output if f.endswith("-tRNAs_name_map.txt"))

valid_trnascan_ids = set()
valid_gtrnadb_ids = set()

# 1. Parse and filter the name map
with open(name_map_in, "r") as fin, open(name_map_out, "w") as fout:
    header = fin.readline()
    fout.write(header)
    for line in fin:
        if not line.strip():
            continue
        parts = line.strip().split("\t")
        if len(parts) >= 2:
            trnascan_id = parts[0]
            gtrnadb_id = parts[1]
            chrom = trnascan_id.split(".")[0]
            if chrom in selected_chromosomes:
                valid_trnascan_ids.add(trnascan_id)
                valid_gtrnadb_ids.add(gtrnadb_id)
                fout.write(line)

# 2. Filter other formats
def filter_fasta(fin_path, fout_path):
    with open(fin_path, "r") as fin, open(fout_path, "w") as fout:
        write_record = False
        for line in fin:
            if line.startswith(">"):
                # GtRNAdb headers usually contain the sequence ID in the first word
                header = line[1:].strip()
                tokens = header.split()
                seq_id = tokens[0] if tokens else ""
                
                matched = False
                # Check against GtRNAdb IDs
                for g_id in valid_gtrnadb_ids:
                    if g_id in seq_id:
                        matched = True
                        break
                
                # Check against tRNAscan IDs
                if not matched:
                    for t_id in valid_trnascan_ids:
                        if t_id in header:
                            matched = True
                            break
                            
                write_record = matched
                
            if write_record:
                fout.write(line)

def filter_bed(fin_path, fout_path):
    with open(fin_path, "r") as fin, open(fout_path, "w") as fout:
        for line in fin:
            if not line.strip():
                continue
            parts = line.strip().split("\t")
            # For BED format, we can check chromosome directly from column 1
            if len(parts) > 0 and parts[0] in selected_chromosomes:
                # Also verify the ID matches valid ones to be extremely safe, but chromosome check is enough for BED.
                fout.write(line)

def filter_out(fin_path, fout_path):
    with open(fin_path, "r") as fin, open(fout_path, "w") as fout:
        for line in fin:
            if line.startswith("Sequence") or line.startswith("Name") or line.startswith("--------"):
                fout.write(line)
                continue
            if not line.strip():
                fout.write(line)
                continue
            parts = line.strip().split()
            if len(parts) >= 2:
                chrom = parts[0]
                trna_num = parts[1]
                t_id = f"{chrom}.trna{trna_num}"
                if t_id in valid_trnascan_ids:
                    fout.write(line)

def filter_ss(fin_path, fout_path):
    with open(fin_path, "r") as fin, open(fout_path, "w") as fout:
        block = []
        write_block = False
        for line in fin:
            if line.strip() == "":
                if write_block:
                    fout.write("".join(block))
                    fout.write(line)
                block = []
                write_block = False
            else:
                if not block:
                    # first line of block: "chr1.trna702 (133033962-133034034)     Length: 73 bp"
                    t_id = line.split()[0]
                    if t_id in valid_trnascan_ids:
                        write_block = True
                block.append(line)
        # flush last block
        if write_block and block:
            fout.write("".join(block))

for in_file, out_file in zip(snakemake.input, snakemake.output):
    if in_file.endswith("_name_map.txt"):
        continue  # Already handled
    elif in_file.endswith(".fa"):
        filter_fasta(in_file, out_file)
    elif in_file.endswith(".bed"):
        filter_bed(in_file, out_file)
    elif in_file.endswith(".out"):
        filter_out(in_file, out_file)
    elif in_file.endswith(".ss"):
        filter_ss(in_file, out_file)
