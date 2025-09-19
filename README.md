# capscan
Scan for Polinton-like virus MCP genes
Scans all gene predictions (.faa files) for PLV or NCLDV-like MCP genes and deposits then in a file.
Usage: Capscan_PLV.sh <genes.faa>
prerequisites: seqkit, HH-suite

Capscan will first scan each fasta amino acid sequence against 50 condensed MCP profiles (PLV_Viro_Yara_NCLDV_MCPs_V3.1.hmm)
Positive hits will then be scanned against MCP clusters (30% identity clusters) to get the closest virus cluster (clusscan file)
