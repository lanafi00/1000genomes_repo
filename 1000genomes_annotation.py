#Iterate through 1000 genomes mutation frequency data 
#Annotate each mutation as CDS (coding sequence), five prime UTR, three prime UTR, exon, intron, or nongenic

import pandas as pd
import os
import itertools

#set counter used to get data for correct chromosome
chrom = os.environ.get('SLURM_ARRAY_TASK_ID')

#Load mutation frequency data, load hapmap data from Ensembl (https://ftp.ensembl.org/pub/release-114/gtf/homo_sapiens/)
file_path = f'1000genomes/dedup_chr{chrom}_1000genomes_freqs2.csv'
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
    pos = int(float(row['POS']))
    annotation = get_annotation(pos,hapmap_chrom)
    annotations.append(annotation)

#Save to csv
mutations_df['annotation'] = annotations
os.makedirs("annotated_mutations", exist_ok=True)
mutations_df.to_csv(f'annotated_mutations/annotated_1000genomes_chrom{chrom}.csv', index=False)
            
