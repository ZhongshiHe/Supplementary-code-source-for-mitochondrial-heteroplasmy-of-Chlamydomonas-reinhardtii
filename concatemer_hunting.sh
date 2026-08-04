# Concatemer hunting in Pacbio HiFi and Oxford Nanopore Technology (ONT) reads

# Map HiFi long reads to the 19kb config 1
minimap2 -cx map-hifi cre_cc5325_mitogenome.fasta hifi_mito_long_gt19k.fastq > hifi_structure_ruler.paf

# Map ONT long reads to the 19kb reference
minimap2 -cx map-ont cre_cc5325_mitogenome.fasta ont_mito_long_gt19k.fastq > ont_structure_ruler.paf

# Extract the Junctions (Skipping Single Monomers)
python parse_mito_junctions.py hifi_structure_ruler.paf > hifi_junctions_ruler.txt
python parse_mito_junctions.py ont_structure_ruler.paf > ont_junctions_ruler.txt

# Combine them
cat hifi_junctions_ruler.txt ont_junctions_ruler.txt > all_junctions_ruler.txt

# FULL formulas with monomers
python parse_mito_structure.py hifi_structure_ruler.paf > hifi_formulas.txt
python parse_mito_structure.py ont_structure_ruler.paf > ont_formulas.txt

# Combine them
cat hifi_formulas.txt ont_formulas.txt > all_formulas.txt

# Summarize into proper format
python summarize_unique_structures.py all_formulas.txt > combined_structure_summary_ruler.txt

# Remove duplication and output results (This awk script reads the junctions file, strips out the gap text (like --[tight seam]-->), renames 19kb to config1 and 15kb to config2, safely deduplicates the 3D read direction (by reverse-complementing the chain), and outputs final CSV.)
tail -n +3 combined_structure_summary_ruler.txt | grep -v "frag" | awk '{
    count = $1;
    $1 = "";
    sub(/^[ \t]+/, "");
    n = split($0, arr, " -> ");

    for (i=1; i<=n; i++) {
        if (arr[i] ~ /19kb/) sub(/19kb/, "config1", arr[i]); 
        else if (arr[i] ~ /15kb/) sub(/15kb/, "config2", arr[i]);
    }

    fwd = arr[1];
    for (i=2; i<=n; i++) fwd = fwd " -> " arr[i];

    rev = "";
    for (i=n; i>=1; i--) {
        comp = arr[i];
        if (comp ~ /\(\+\)/) sub(/\(\+\)/, "(-)", comp);
        else if (comp ~ /\(-\)/) sub(/\(-\)/, "(+)", comp);
        
        if (rev == "") rev = comp;
        else rev = rev " -> " comp;
    }

    canonical = (fwd < rev) ? fwd : rev;
    sums[canonical] += count;

} END {
    print "Count,Structure";
    for (s in sums) {
        print sums[s] "," s;
    }
}' | sort -t, -k1,1nr > finalized_structure_counts_FULL.csv
