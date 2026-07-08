import sys
from collections import defaultdict

if len(sys.argv) != 2:
    print("Usage: python parse_mito_junctions.py <input.paf>")
    sys.exit(1)

paf_file = sys.argv[1]
reads = defaultdict(list)

# Parse the PAF file
with open(paf_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 11: 
            continue
        
        qname = parts[0]
        qstart = int(parts[2])
        qend = int(parts[3])
        strand = parts[4]
        
        aln_len = qend - qstart
        
        # Filter out tiny fragments
        if aln_len < 8000: 
            continue
            
        if aln_len > 17500:
            unit = "19kb"
        elif aln_len >= 12000 and aln_len <= 17500:
            unit = "15kb"
        else:
            unit = f"frag_{aln_len//1000}k"
            
        reads[qname].append({'start': qstart, 'end': qend, 'unit': unit, 'strand': strand})

# Process and print reads that have MORE THAN ONE block (concatemers)
for qname, blocks in reads.items():
    if len(blocks) < 2:
        continue # Skip single-unit reads since they have no internal junctions
        
    sorted_blocks = sorted(blocks, key=lambda x: x['start'])
    
    structure_parts = []
    for i in range(len(sorted_blocks)):
        b = sorted_blocks[i]
        structure_parts.append(f"{b['unit']}({b['strand']})")
        
        # If there is a next block, calculate the gap between the current block's end and next block's start
        if i < len(sorted_blocks) - 1:
            next_b = sorted_blocks[i+1]
            gap = next_b['start'] - b['end']
            
            # Formatting the gap string
            if gap > 100:
                gap_str = f" --[UNCLEAR GAP: {gap}bp]--> "
            elif gap < -100:
                gap_str = f" --[OVERLAP: {gap}bp]--> "
            else:
                gap_str = f" --[tight seam: {gap}bp]--> "
                
            structure_parts.append(gap_str)
            
    print(f"{qname}\t{''.join(structure_parts)}")
