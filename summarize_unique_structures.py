import sys
from collections import defaultdict
import re

def reverse_complement(structure):
    """
    Takes a tuple of units like ('19kb(+)', '15kb(+)')
    and returns its reverse complement: ('15kb(-)', '19kb(-)')
    """
    rev = []
    # Reverse the order of the units
    for unit in reversed(structure):
        # Flip the strand orientation
        if '(+)' in unit:
            rev.append(unit.replace('(+)', '(-)'))
        elif '(-)' in unit:
            rev.append(unit.replace('(-)', '(+)'))
        else:
            rev.append(unit)
    return tuple(rev)

def get_canonical_form(structure):
    """
    Returns the lexicographically smaller version of the structure 
    or its reverse complement to act as the standard 'Unique ID'.
    """
    rev_comp = reverse_complement(structure)
    # Alphabetical comparison just to pick one consistent representation
    return structure if structure < rev_comp else rev_comp

def main():
    if len(sys.argv) < 2:
        print("Usage: python summarize_unique_structures.py <file1.txt> <file2.txt> ...")
        sys.exit(1)

    counts = defaultdict(int)

    for filename in sys.argv[1:]:
        with open(filename, 'r') as f:
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) < 2:
                    continue
                
                # parts[1] contains the string like "15kb(+) --[tight seam: -17bp]--> 19kb(-)"
                # This regex extracts just the units: e.g., ['15kb(+)', '19kb(-)']
                units = tuple(re.findall(r'[a-zA-Z0-9_]+\([+-]\)', parts[1]))
                
                if len(units) > 0:
                    # Group symmetrical reads together
                    canonical = get_canonical_form(units)
                    counts[canonical] += 1

    # Print a nicely formatted summary table
    print(f"{'Count':<8} {'Unique Canonical Structure'}")
    print("-" * 60)
    
    # Sort by highest frequency
    for structure, count in sorted(counts.items(), key=lambda x: x[1], reverse=True):
        struct_str = " -> ".join(structure)
        print(f"{count:<8} {struct_str}")

if __name__ == '__main__':
    main()
