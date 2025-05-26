#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#Creating dataframe of mutations from 1000 genomes dataset present in populations of interest
#Later remove rows duplicate for CHROM and POS (repeats) in the terminal
import pandas as pd
import gzip
import re
import os 

#set counter used to get data for correct chromosome
i = os.environ.get('SLURM_ARRAY_TASK_ID')

#Load csv of samples
samples_df = pd.read_csv('integrated_call_samples_v3.20130502.ALL.panel', sep='\t')

#make lists of CEU, CHB, YRI sample IDs
ceu_samples = samples_df[samples_df['pop'] == 'CEU']['sample'].tolist()
chb_samples = samples_df[samples_df['pop'] == 'CHB']['sample'].tolist()
yri_samples = samples_df[samples_df['pop'] == 'YRI']['sample'].tolist()


#Make list to store dicitonary of each variant site and its frequency per each population
results = []

# Debug counters
total_variants = 0
non_snp = 0
too_many_alts = 0
indels = 0
row_count = 0

#Iterate through each variant site
with gzip.open(f'ALL.chr{i}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz', 'rt') as vcf:
    
    for line in vcf:
        if line.startswith('##'):
            continue
        if line.startswith('#CHROM'):
            header = line.strip().split('\t')
            sample_columns = header[9:]  # first 9 columns are VCF fixed fields
            ceu_indices = [i for i, s in enumerate(sample_columns) if s in ceu_samples]
            chb_indices = [i for i, s in enumerate(sample_columns) if s in chb_samples]
            yri_indices = [i for i, s in enumerate(sample_columns) if s in yri_samples]
            break  # exit header parsing

    for line in vcf:
        fields = line.strip().split('\t')
        ref = fields[3]
        alt = fields[4]
        info = fields[7]
        genotypes = fields[9:]
        #Filter only for SNPs with a single altnerative allele
        if len(ref) != 1:
            non_snp += 1
            continue
        alts = alt.split(',')
        if (len(alts) > 2):
            too_many_alts += 1
            continue
        if any(len(a) > 1 for a in alts):
            indels += 1
            continue

        #Create dictionary to store each alternative allele and its index
        alt_dict = {alts[i]:str(i+1) for i in range(len(alts))}
       
            

        #Counts up the number of altnerative alleles among samples of given indices
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

        freqs = {}
        for a in alts:
            alt_code = alt_dict[a]
            ceu_count = count_alt(ceu_indices, alt_code)
            chb_count = count_alt(chb_indices, alt_code)
            yri_count = count_alt(yri_indices, alt_code)

            ceu_freq = ceu_count / (2 * len(ceu_indices)) if ceu_indices else 'NA'
            chb_freq = chb_count / (2 * len(chb_indices)) if chb_indices else 'NA'
            yri_freq = yri_count / (2 * len(yri_indices)) if yri_indices else 'NA'

            freqs[f'CEU_FREQ_{a}'] = ceu_freq
            freqs[f'CHB_FREQ_{a}'] = chb_freq
            freqs[f'YRI_FREQ_{a}'] = yri_freq
            
        row = {'CHROM':fields[0], 
              'POS':fields[1],
               'ANCESTRAL': get_ancestral(info)}
        row.update(freqs)
        results.append(row)
        row_count += 1

df = pd.DataFrame(results)
print(f"chr{i} - Variants written: {row_count}")
print(f"chr{i} - Non-SNPs filtered: {non_snp}")
print(f"chr{i} - Too many ALT alleles filtered: {too_many_alts}")
print(f"chr{i} - Indels filtered: {indels}")

#Save csv for given chromosome. Later combine all csvs in the terminal. 
df.to_csv(f'chr{i}_1000genomes_freqs.csv',index=False)


# In[ ]:


#Adding Ancestral frequency - Alt frequency columns
# Get SLURM array task ID
i = os.environ.get('SLURM_ARRAY_TASK_ID')
df = pd.read_csv(f"chr{i}_1000genomes_freqs.csv")


#list populations and nucleotides
populations = ['CEU', 'CHB', 'YRI']
alleles = ['A', 'C', 'G', 'T']

# Build list of all alt frequency columns that exist in the file
freq_cols = [f"{pop}_FREQ_{allele}" for pop in populations for allele in alleles if f"{pop}_FREQ_{allele}" in df.columns]

# Count number of non-null alt frequencies per site (i.e., per row)
df['ALT_FREQ_COUNT'] = df[freq_cols].notnull().sum(axis=1)

# Classify SNPs
df['SNP_TYPE'] = df['ALT_FREQ_COUNT'].map({3: 'Biallelic', 6: 'Triallelic'})

#Add Ancestral frequency - Alt frequency columns
for pop in populations:
    alt_freq_cols = [f"{pop}_FREQ_{allele}" for allele in alleles if f"{pop}_FREQ_{allele}" in df.columns]

    #Add up alt frequencies (should be 2 if triallelic, 1 if biallelic) 
    df[f"{pop}_ALT_SUM"] = df[alt_freq_cols].sum(axis=1, skipna=True)

    #Subtract all alt frequencies from 1 to get ancestral frequency
    df[f"{pop}_ANCESTRAL_FREQ"] = 1 - df[f"{pop}_ALT_SUM"]

    #Creating column for the Ancestral Frequency - Alternate Frequency for each alt allele present
    for allele in alleles:
        freq_col = f"{pop}_FREQ_{allele}"
        diff_col = f"{pop}_DIFF_{allele}"
        if freq_col in df.columns:
            df[diff_col] = df[f"{pop}_ANCESTRAL_FREQ"] - df[freq_col]
        
# Drop temporary ALT_SUM columns  
df.drop(columns=[f"{pop}_ALT_SUM" for pop in populations] + ['ALT_FREQ_COUNT'], inplace=True)
        
# Save updated CSV
df.to_csv(f"chr{i}_1000genomes_freq_differences.csv", index=False)

# Print summary counts
print("Summary of SNP types:")
print(df['SNP_TYPE'].value_counts(dropna=False))


# In[ ]:


#Filtering out data where the Ancestral allele has a frequency of 1.0 for all populations sampled (fixed mutations)

# Get SLURM array task ID
i = os.environ.get('SLURM_ARRAY_TASK_ID')

# Load CSV
df = pd.read_csv(f'chr{i}_1000genomes_freq_diffs1.csv')

# List of ancestral freq columns
pops_ancestral = ['YRI_ANCESTRAL_FREQ', 'CEU_ANCESTRAL_FREQ','CHB_ANCESTRAL_FREQ']

# Create mask: True for rows where all three ancestral freqs == 1.0
mask = (df[pops_ancestral] == 1.0).all(axis=1)

# Drop those rows
df_filtered = df[~mask]

# Save result
df_filtered.to_csv(f'chr{i}_1000genomes_freq_diffs_filtered.csv', index=False)
print(f"Filtered chromosome {i}: dropped {mask.sum()} rows.")


# In[ ]:


#Adding (YRI ancestral frequency - YRI alt frequency) - (CHB/CEU ancestral frequency - CHB/CEU alt frequency) columns
#In retrospect this was a somewhat roundabout way to accomplish this but what's done is done

import pandas as pd
import os 

#set counter used to get data for correct chromosome
i = os.environ.get('SLURM_ARRAY_TASK_ID')

#Load csv
df = pd.read_csv(f'chr{i}_1000genomes_freq_diffs_filtered.csv')
alleles = ['A', 'C', 'G', 'T']


#For each allele present as an alt allele, create a column (YRI_DIFF) - (CHB_DIFF) and a column (YRI_DIFF) - (CEU_DIFF) 
for allele in alleles:
    diff_yri_col = f'YRI_DIFF_{allele}'
    diff_ceu_col = f'CEU_DIFF_{allele}'
    diff_chb_col = f'CHB_DIFF_{allele}'
    freq_yri_col = f'YRI_FREQ_{allele}'
    freq_ceu_col = f'CEU_FREQ_{allele}'
    freq_chb_col = f'CHB_FREQ_{allele}'
    #Find YRI_FREQ_DIFF - CEU_FREQ_DIFF, and YRI_FREQ_DIFF - CHB_FREQ_DIFF for each population

    if diff_yri_col in df.columns and diff_ceu_col in df.columns:
        df[f'(YRI-CEU)_DIFF_{allele}'] = df[diff_yri_col] - df[diff_ceu_col]

    if diff_yri_col in df.columns and diff_chb_col in df.columns:
        df[f'(YRI-CHB)_DIFF_{allele}'] = df[diff_yri_col] - df[diff_chb_col]

    if freq_yri_col in df.columns and freq_chb_col in df.columns:
        df[f'(YRI-CHB)_{allele}'] = df[freq_yri_col] - df[freq_chb_col]

    if freq_yri_col in df.columns and freq_ceu_col in df.columns:
        df[f'(YRI-CEU)_{allele}'] = df[freq_yri_col] - df[freq_ceu_col]


#Save to new csv file
df.to_csv(f'chr{i}_1000genomes_freq_diffs1.csv',index=False)

