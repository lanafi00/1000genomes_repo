#Storing position, ancestral allele, alternative alleles, and frequency of alternative allele, as well as per-population frequency differences
#Then filters out sites fixed for ancestral allele in all 3 pops of interest

import pandas as pd
import gzip
import re
import os 
import numpy as np
import pickle

#set counter used to get data for correct chromosome
i = os.environ.get('SLURM_ARRAY_TASK_ID')

#Load csv of samples if index file exists -- if not then create it
index_file = 'sample_indices.pkl'
if os.path.exists(index_file):
    with open(index_file, 'rb') as f:
        index_map = pickle.load(f)
    ceu_samples = index_map['CEU']
    chb_samples = index_map['CHB']
    yri_samples = index_map['YRI']
else:
    samples_df = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel', sep='\t')
    #make lists of CEU, CHB, YRI sample IDs
    ceu_samples = samples_df[samples_df['pop'] == 'CEU']['sample'].tolist()
    chb_samples = samples_df[samples_df['pop'] == 'CHB']['sample'].tolist()
    yri_samples = samples_df[samples_df['pop'] == 'YRI']['sample'].tolist()
    index_map = {'CEU': ceu_samples, 'CHB': chb_samples, 'YRI': yri_samples}
    with open(index_file, 'wb') as f_out:
        pickle.dump(index_map, f_out)

#Make list to store dicitonary of each variant site and its frequency per each population
results = []

# Debug counters
total_variants = 0
non_snp = 0
too_many_alts = 0
indels = 0
row_count = 0

#Iterate through each variant site present in 1000 Genomes dataset (https://hgdownload.cse.ucsc.edu/gbdb/hg19/1000Genomes/)
with gzip.open(f'ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'rt') as vcf:   
    for line in vcf:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            sample_columns = header[9:]  # first 9 columns are VCF fixed fields
            ceu_indices = [idx for idx, s in enumerate(sample_columns) if s in ceu_samples]
            chb_indices = [idx for idx, s in enumerate(sample_columns) if s in chb_samples]
            yri_indices = [idx for idx, s in enumerate(sample_columns) if s in yri_samples]
            break  # exit header parsing

    for line in vcf:
        fields = line.strip().split('\t')
        ref = fields[3]
        alt = fields[4]
        info = fields[7]
        genotypes = fields[9:]
        #Filter only for biallelic SNPs with a single alternative allele
        if len(ref) != 1:
            non_snp += 1
            continue
        alts = alt.split(',')
        if (len(alts) > 1):
            too_many_alts += 1
            continue
        if any(len(a) > 1 for a in alts):
            indels += 1
            continue
            
        #Counts up the number of samples with the alternative allele among given sample indexes
        def count_alt(indices, alt_code):
            count = 0
            for idx in indices:
                gt = genotypes[idx].split(':')[0]
                alleles = gt.replace('|', '/').split('/')
                count += alleles.count(alt_code)
            return count

        #Gets the ancestral allele for a given line
        def get_ancestral(info):
            match = re.search(r'AA=([ACGT])', info)
            return match.group(1) if match else 'NA'

        #Gets the frequency of alt allele for each population of interest (YRI, CEU, CHB) 
        #Then calculate the difference between YRI and CEU, and YRI and CHB alt allele frequencies
        freqs = {}
        for a in alts:
            #Because we're only looking at biallelic sites, alt_code is always 1
            alt_code = "1"
            ceu_count = count_alt(ceu_indices, alt_code)
            chb_count = count_alt(chb_indices, alt_code)
            yri_count = count_alt(yri_indices, alt_code)

            ceu_freq = ceu_count / (2 * len(ceu_indices)) if ceu_indices else np.nan
            chb_freq = chb_count / (2 * len(chb_indices)) if chb_indices else np.nan
            yri_freq = yri_count / (2 * len(yri_indices)) if yri_indices else np.nan

            freqs['CEU_FREQ_ALT'] = ceu_freq
            freqs['CHB_FREQ_ALT'] = chb_freq
            freqs['YRI_FREQ_ALT'] = yri_freq
            freqs['YRI-CEU'] = yri_freq - ceu_freq 
            freqs['YRI-CHB'] = yri_freq - chb_freq
             
        row = {'CHROM':fields[0], 
              'POS':fields[1],
                'ANCESTRAL': get_ancestral(info),
              'ALT':alts[0]}
        row.update(freqs)
        results.append(row)
        row_count += 1

#Save to df
df = pd.DataFrame(results)

# Create mask: True for rows where all three ancestral freqs == 1.0
mask = (df['CEU_FREQ_ALT'] == 0.0) & (df['CHB_FREQ_ALT'] == 0.0) & (df['YRI_FREQ_ALT'] == 0.0)

# Drop those rows
df_filtered = df[~mask]

os.makedirs("1000genomes", exist_ok=True)
df_filtered.to_csv(f'1000genomes/chr{i}_1000genomes_freqs2.csv',index=False)

#Rows that correspond to the same position in the genome (['POS'] field) are removed on the cluster
print(f"chr{i} - Variants written: {row_count}")
print(f"chr{i} - Non-SNPs filtered: {non_snp}")
print(f"chr{i} - Too many ALT alleles filtered: {too_many_alts}")
print(f"chr{i} - Indels filtered: {indels}")
