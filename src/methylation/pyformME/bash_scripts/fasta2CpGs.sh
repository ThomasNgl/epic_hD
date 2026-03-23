#!/bin/bash

usage() {
    cat <<EOF
Usage: $0 -fa FASTA_PATH -chr CHR_LIST -o SAVE_PATH [-w WINDOW]

Options:
  -fa, --fasta_path   Path to the fasta file
  -chr, --chr_list    Comma-separated list of chromosome names (default: all)
  -o, --save_path     Path to save the results
  -w, --window        Window size (default: 1000)
  -h, --help          Show this help message and exit

Example:
  $0 -fa genome.fa -chr chr1,chr2 -o out_dir -w 1000
EOF
    exit 1
}

# Default values
chr_list=""
window=1000

# Parse arguments
while [[ "$#" -gt 0 ]]; do
    case $1 in
        -fa|--fasta_path)
            fasta_path="$2"
            shift 2
            ;;
        -chr|--chr_list)
            chr_list="$2"
            shift 2
            ;;
        -o|--save_path)
            save_path="$2"
            shift 2
            ;;
        -w|--window)
            window="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown parameter passed: $1"
            usage
            ;;
    esac
done

# Check for required arguments
if [[ -z "$fasta_path" || -z "$save_path" ]]; then
    echo "Missing required arguments!"
    usage
fi

SECONDS=0

# Prepare chr_name_list for Python
if [[ -z "$chr_list" ]]; then
    chr_arg=None
else
    chr_arg="'$chr_list'.split(',')"
fi

cmd="fasta2CpG.fasta2CpGs(fasta_file='$fasta_path', save_path='$save_path', window=$window, chr_name_list=$chr_arg)"

# Run the Python script/function
echo "[$(date)]: Starting command: $cmd"

python -c "
import fasta2CpG
fasta2CpG.fasta2CpG(
    fasta_file='$fasta_path',
    save_path='$save_path',
    window=$window,
    chr_name_list=$chr_arg
)
"

time_elapsed="$SECONDS"
echo "[$(date)]: Command successful. Total time elapsed: $(( $time_elapsed / 3600)) h $(( ($time_elapsed / 60) % 60)) m $(( $time_elapsed % 60 )) s."
