# capscan
Scan for Polinton-like virus MCP genes, classify them, deposit in a file.

Usage: capscan.sh <amino_acids.faa>
prerequisites: seqkit, HH-suite

Capscan will first scan each fasta amino acid gene prediction against 50 condensed MCP profiles (PLV_Viro_Yara_NCLDV_MCPs_V3.1.hmm).
Positive hits will then be scanned against MCP clusters (30% identity clusters) to determine the closest virus cluster (clusscan file).
