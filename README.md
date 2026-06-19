# capscan

Scan for Polinton-like virus, Mriyavirus, Polinton, NCLDV and related MCP genes.
Classifies them and deposits the results into a file.

## Usage

```bash
capscan.sh -c <number of cores> <amino_acids.faa>
```

## Prerequisites

- [seqkit](https://bioinf.shenwei.me/seqkit/)
- [HH-suite](https://github.com/soedinglab/hh-suite)

## How it works

1. Capscan first scans each FASTA amino acid gene prediction against about 50
   condensed MCP profiles (`PLV_Viro_Yara_NCLDV_MCPs_V3.1.hmm`).
2. Positive hits are then scanned against MCP clusters (30% identity clusters)
   to determine the closest virus cluster (`clusscan` file).
3. A `consolidated.tsv` is generated showing the closest viral group.
