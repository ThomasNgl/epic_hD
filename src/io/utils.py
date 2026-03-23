def get_chr_names():
    # Mouse chromosomes: 1-19, X, Y, M
    chrs = [f"chr{i}" for i in range(1, 20)]  # chr1 to chr19
    chrs += ["chrX", "chrY", "chrM"]
    return chrs