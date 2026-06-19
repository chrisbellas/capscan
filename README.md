# capscan
Scan for Polinton-like virus, Mriyavirus, Polinton, NCLDV and related MCP genes.
Classifies them and deposits into a file.

Usage: capscan.sh -c <number of cores> <amino_acids.faa>
prerequisites: seqkit, HH-suite

Capscan will first scan each fasta amino acid gene prediction against about 50 condensed MCP profiles (PLV_Viro_Yara_NCLDV_MCPs_V3.1.hmm).
Positive hits will then be scanned against MCP clusters (30% identity clusters) to determine the closest virus cluster (clusscan file).
A consolidated.tsv will be generated showing the closest viral group. 
