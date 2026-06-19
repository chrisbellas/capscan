#!/bin/bash
# capscan: Scans all gene predictions (.faa files) for PLV or NCLDV-like MCP genes, classifies the PLV and deposits then in a file.
# Usage: capscan.sh <genes.faa>
# prerequisites: seqkit, HH-suite

VERSION="0.2.0"
echo "capscan version $VERSION"

# change location of HMM profile databases
MCPSCAN="PLV_Viro_Yara_NCLDV_MCPs_V3.1.hmm" # 56 profiles to scan for MCP genes (quick)
MCPCLUSTER="VC_MCPs.hmm" # 293 profiles to determine closest viral cluster (Slow and runs on results from above) (Includes 59 Tectivirus profiles - just for fun)
MAJORGROUPS="Major_groups.tsv" # maps each major HMM profile to Virus Phylum/Class and Virus group

# set number of CPUs to use (default: 4, override with -c)
CPU=4

# parse options
while getopts ":c:h" opt; do
    case $opt in
        c) CPU=$OPTARG ;;
        h)
            echo "Usage: capscan.sh [-c cores] <genes.faa>"
            echo "Scans all gene predictions (.faa files) for PLV or NCLDV-like MCP genes, classifies the PLV and deposits them in a file"
            exit 0
            ;;
        :)  echo "Error: option -$OPTARG requires an argument (number of cores)." >&2; exit 1 ;;
        \?) echo "Error: invalid option -$OPTARG." >&2; exit 1 ;;
    esac
done
shift $((OPTIND - 1))

#check for input file
if [ -z "$1" ]; then
    echo "Usage: capscan.sh [-c cores] <genes.faa>"
    echo "Scans all gene predictions (.faa files) for PLV or NCLDV-like MCP genes, classifies the PLV and deposits them in a file"
    exit 1
fi


# Check for seqkit
if ! command -v seqkit &> /dev/null; then
    echo "Error: seqkit is not installed or not in PATH."
    exit 1
fi

# Check for hmmsearch
if ! command -v hmmsearch &> /dev/null; then
    echo "Error: hmmsearch is not installed or not in PATH."
    exit 1
fi

#strip basename
base="${1%.*}"

# scan for capsids
echo 'Scanning for MCP genes'
hmmsearch --tblout $base'_capscan.tblout' --domtblout $base'_capscan.domtblout' -o $base'_capscan_out' -E 1e-9 --cpu $CPU $MCPSCAN $1


# counting MCP hits
count=$(cat $base'_capscan.tblout' | awk '{print $1}' | grep -v '#' | sort -nk 1 | uniq | wc -l)
echo $count ' MCP genes found'

# Pulling MCP genes with seqkit into intermediate file
awk '{print $1}' $base'_capscan.tblout' | sort -nk 1 | uniq | seqkit grep -f - $1 >$base'_MCP.faa'

#parse besthits for MCP gene from hmmsearch results
awk '{print $1, $6, $0}' $base'_capscan.tblout' | sort -k1,1 -k2,2nr |  awk '!seen[$1]++ {$1=""; $2=""; print $0}' | awk '{print $1,$2,$3,$4,$5,$6}' | tr ' ' '\t' >$base'_capscan_besthit.tblout'
echo 'First pass hmmsearch results: '$base'_capscan_besthit.tblout'

# Separating PLVs and GVs
grep -v '#' $base'_capscan_besthit.tblout' |grep GVOG | awk '{print $1}' >$base'_GVMCPS.list'
grep -v '#' $base'_capscan_besthit.tblout' |grep -v GVOG | awk '{print $1}' >$base'_PLV.list'

# filter for hits over 75 bitscore
cat $base'_capscan_besthit.tblout'  | grep -v '#' |grep -v GVOG | awk '$6>=75 {print $1}' >$base'_PLV75.list'

# create list of bad hits - probably junk
if [ -s "$base"_PLV75.list ] && [ -s "$base"_PLV.list ]; then
	grep -v -f "$base"_PLV75.list "$base"_PLV.list >"$base"_PLV_weak_hits.list
else
	# If either file is empty, create an empty weak hits list
	>"$base"_PLV_weak_hits.list
fi

echo "Pulling MCP genes"
seqkit grep --quiet -f $base'_PLV_weak_hits.list' $base'_MCP.faa' >$base'_PLV_junk_hits.faa'
seqkit grep --quiet -f $base'_GVMCPS.list' $base'_MCP.faa' >$base'_GVMCP.faa'
seqkit grep --quiet -f $base'_PLV75.list' $base'_MCP.faa' >$base'_PLVMCP_strong_hits.faa'

count1=$(grep -c ">" "${base}_PLVMCP_strong_hits.faa")
count2=$(grep -c ">" "${base}_GVMCP.faa")

echo "Results:"
echo "$count1 PLV MCPs detected"
echo "$count2 NCLDV MCPs detected"


echo 'Finding MCP cluster from hits'
#Find closest MCP cluster

if [ "$count1" -gt 0 ]; then
		#Perform cluster search
	hmmsearch --tblout $base'_clusscan.tblout' -o $base'_clusscan_out' -E 1e-9 --cpu $CPU $MCPCLUSTER $base'_PLVMCP_strong_hits.faa'
	#parse besthits from  MCP cluster scan
	awk '{print $1, $6, $0}' $base'_clusscan.tblout'  | sort -k1,1 -k2,2nr |  awk '!seen[$1]++ {$1=""; $2=""; print $0}' | awk '{print $1,$2,$3,$4,$5,$6}' | tr ' ' '\t' >$base'_clusscan_besthit.tblout'
	else
		# continue
		echo "PLV Count is 0, Skipping cluster search"
fi

# Cluster results
echo "Virus major group hits"
grep -v '#' $base'_capscan_besthit.tblout' | awk '{print $3}' | sort -nk 1 | uniq -c

#echo "Virus clusters counts"
#grep -v '#' $base'_clusscan_besthit.tblout' | awk '{print $3}' | sort -nk 1 | uniq -c

echo $base'_capscan_besthit.tblout  - capsid supergroup hit file  (evalue, bitscore)'
echo $base'_clusscan_besthit.tblout - closest Virus Cluster  (evalue, bitscore)'
echo $base'_PLVMCP_strong_hits.faa  - PLV MCPs'

# Consolidate results
grep -v '#' $base'_capscan_besthit.tblout' | awk '{print $1,$3,$5,$6}' | tr ' ' '\t' >Round1
grep -v '#' $base'_clusscan_besthit.tblout' | awk '{print $1,$3,$5,$6}' | tr ' ' '\t' >Round2

# Merge Round1 (major HMM profile) and Round2 (virus cluster) on MCP_hit key,
# cross-referencing $MAJORGROUPS to pull Virus Phylum/Class and Virus group for each profile
{
	printf 'MCP_hit\tMajor_HMM_profile\tVirus_Phylum_Class\tVirus_group_hit\tevalue\tbitscore\tVirus_cluster_hit\tCluster_evalue\tCluster_bitscore\n'
	awk -F'\t' -v groups="$MAJORGROUPS" '
		FILENAME==groups{grp[$1]=$2"\t"$3; next}
		FILENAME=="Round2"{r2[$1]=$2"\t"$3"\t"$4; next}
		{print $1"\t"$2"\t"(($2 in grp)?grp[$2]:"NA\tNA")"\t"$3"\t"$4"\t"(($1 in r2)?r2[$1]:"NA\tNA\tNA")}
	' "$MAJORGROUPS" Round2 Round1
} >$base'_consolidated.tsv'
echo $base'_consolidated.tsv         - consolidated results table'

#cleanup
rm -f $base'_PLV_weak_hits.list' $base'_GVMCPS.list' $base'_PLV.list' $base'_PLV75.list' $base'_MCP.faa' Round1 Round2
