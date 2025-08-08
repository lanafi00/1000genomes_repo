#Iterate through simulation analysis data 
#Annotate each mutation as CDS (coding sequence), five prime UTR, three prime UTR, exon, intron, or nongenic

import pandas as pd
import os
import itertools

#Set up list of model ID codes and dfe ID codes which correspond to file name
model_list = ["3I21",
              "3G09",
              "2T12", 
              "9K19", 
              "10J19",
              "4J17"]
dfe_list = ["Z21",
            "K23",
            "K17"]

#Get SLURM task ID
task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))  # defaults to 0 if not running in array.

# Build a list of all combinations (model_id, dfe_id, chromosome)
all_tasks = list(itertools.product(model_list, dfe_list, range(1, 23)))

# Get corresponding parameters
model_id, dfe_id, chrom = all_tasks[task_id]

#Check whether sim has already been annotated, if not then annotate it; this way script can be rerun if some annotations fail
newfile_path = f'annotated_mutations/annotated_{dfe_id}_{model_id}_chr{chrom}.csv'
if os.path.isfile(newfile_path) and os.path.getsize(newfile_path) > 0:
    print(f"Skipping {newfile_path}: already exists and is non-empty.")
else:
    file_path = f"all_chrom_outputs/{dfe_id}_{model_id}_chr{chrom}_freq_differences4.csv"
    
    #Load mutation data, load hapmap data from Ensembl (https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/)
    mutations_df = pd.read_csv(file_path)
    hapmap_df = pd.read_csv('Homo_sapiens.GRCh38.114.chr.gtf.gz', sep='\t', comment='#', header=None,
                            names=['seqname','source','feature','start','end','score','strand','frame','attribute'])
    
    hapmap_df['seqname'] = hapmap_df['seqname'].astype(str)
    
    #Focus on chromosome of interest
    hapmap_chrom = hapmap_df[hapmap_df['seqname'] == str(chrom)]
    
    #Takes input of position, and genome map, and returns the annotation at that position
    def get_annotation(pos, hapmap_chrom):
        hapmap_pos = hapmap_chrom[(hapmap_chrom['start'] <= pos) & (hapmap_chrom['end'] >= pos)]
        hits = hapmap_pos['feature'].unique().tolist()
        if not hits:
            return "nongenic"
        if "CDS" in hits:
            return "CDS"
        if "five_prime_utr" in hits:
            return "five_prime_utr"
        if "three_prime_utr" in hits:
            return "three_prime_utr"
        if "exon" in hits:
            return "exon"
        if "transcript" in hits:
            return "intron"
        else:
            return "other"
         
    #Add annotation for each mutation        
    annotations = []
    for idx, row in mutations_df.iterrows():
        pos = int(float(row['position']))
        annotation = get_annotation(pos,hapmap_chrom)
        annotations.append(annotation)

    #Save to csv
    mutations_df['annotation'] = annotations
    os.makedirs("annotated_mutations", exist_ok=True)
    mutations_df.to_csv(f'annotated_mutations/annotated_{dfe_id}_{model_id}_chr{chrom}.csv', index=False)
                
