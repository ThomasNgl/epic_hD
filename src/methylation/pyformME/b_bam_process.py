from collections import defaultdict

def deduplicate_reads(reads_list):
    dupe_groups = defaultdict(list)
    for read in reads_list:
        read_chr_name = read.reference_name
        read_start = read.reference_start
        read_end = read.reference_end
        read_length = len(read.query_sequence)
        read_direction = read.is_reverse
        key = (read_chr_name, read_start, read_end, read_length, read_direction)
        dupe_groups[key].append(read)
    
    unique_reads = []
    markdup = 0
    for group in dupe_groups.values():
        best_read = max(group, key=lambda r: r.mapping_quality)
        unique_reads.append(best_read)
        if not best_read.is_duplicate:
            markdup += 1 
    n_reads = len(reads_list)
    n_sg_reads = len(unique_reads)
    dup_qc = {"n_reads":n_reads, "n_sg_reads":n_sg_reads, "markdup":markdup}
    return unique_reads, dup_qc