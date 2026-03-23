import os
from . import io_informME
from . import a_refCpGs_process

def fasta2CpGs_chr(fasta_file,  
                  out_dir,
                  chr_name,
                  window=1000):
    # Load the fasta file in case it is a path and not a fasta obj
    if isinstance(fasta_file, str):
        fasta_file = io_informME.load_fasta(fasta_file)
    # else: it's a valid Fasta object, continue as normal

    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    results = {}
    # Iterate across chromosomes 
    print(f"Process {chr_name}")
    df_CpGs = a_refCpGs_process.get_CpGs_chr(fasta_file, chr_name)
    df_CpGs = a_refCpGs_process.get_CpGs_distance(df_CpGs)
    df_CpGs = a_refCpGs_process.get_CpGs_density(df_CpGs, window=window)
    results[chr_name] = df_CpGs
    print(f"CpGs site positions, distance and density computed for {chr_name}")
    # Save per chromosome as CSV in the specified path
    out_file = os.path.join(out_dir, f"CpGslocation_{chr_name}.csv")
    df_CpGs.to_csv(out_file, index=False)
    print(f"Saved in {os.path.abspath(out_file)}")
    return results


def fasta2CpGs(fasta_file,  
                  out_dir,
                  chr_name_list=None,
                  window=1000):
    """:
    """
    # Load the fasta file in case it is a path and not a fasta obj
    if isinstance(fasta_file, str):
        fasta_file = io_informME.load_fasta(fasta_file)
    # else: it's a valid Fasta object, continue as normal

    # Check the list of chromosomes
    if chr_name_list is None:
        print("Will process all the chromosomes found in the fasta file")
        chr_name_list = list(fasta_file.references)
        print(*chr_name_list, sep='\n')
    
    # Create output directory if it doesn't exist
    os.makedirs(out_dir, exist_ok=True)
    
    results = {}
    # Iterate across chromosomes 
    for chr_name in chr_name_list:
        print(f"Process {chr_name}")
        df_CpGs = a_refCpGs_process.get_CpGs_chr(fasta_file, chr_name)
        if df_CpGs.empty:
            print(f'No CpG sites found in {chr_name}')
            continue
        df_CpGs = a_refCpGs_process.get_CpGs_distance(df_CpGs)
        df_CpGs = a_refCpGs_process.get_CpGs_density(df_CpGs, window=window)
        results[chr_name] = df_CpGs
        print(f"CpGs site positions, distance and density computed for {chr_name}")
        # Save per chromosome as CSV in the specified path
        out_file = os.path.join(out_dir, f"CpGslocation_{chr_name}.csv")
        df_CpGs.to_csv(out_file, index=False)
        print(f"Saved in {os.path.abspath(out_file)}")
    return results

    