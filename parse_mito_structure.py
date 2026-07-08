import sys
from collections import defaultdict

if len(sys.argv) != 2:
    print("Usage: python parse_mito_structure.py <input.paf>")
    sys.exit(1)

paf_file = sys.argv[1]
reads = defaultdict(list)

# Parse the PAF file
with open(paf_file, 'r') as f:
    for line in f:
        parts = line.strip().split('\t')
        if len(parts) < 11: 
            continue
        
        qname = parts[0]        # Read name
        qstart = int(parts[2])  # Start coord on read
        qend = int(parts[3])    # End coord on read
        strand = parts[4]       # + or -
        
        aln_len = qend - qstart
        
        # Filter out tiny fragments (noise/edges)
        if aln_len < 8000: 
            continue
            
        # Classify the unit based on how much of the read it consumes
        # The 19kb unit is ~19,039bp. The 15kb unit is ~15,107bp.
        if aln_len > 17500:
            unit = "19kb"
        elif aln_len >= 12000 and aln_len <= 17500:
            unit = "15kb"
        else:
            unit = f"frag_{aln_len//1000}k"
            
        reads[qname].append((qstart, unit, strand))

# Sort alignments by position on the read and print
for qname, blocks in reads.items():
    # Sort by the 5' -> 3' start coordinate on the read
    sorted_blocks = sorted(blocks, key=lambda x: x[0])
    
    # Format the structure string
    structure = " -> ".join([f"{b[1]}({b[2]})" for b in sorted_blocks])
    
    # Print Read ID and its structural formula
    print(f"{qname}\t{structure}")
