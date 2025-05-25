#!/usr/bin/env python
# coding: utf-8

# In[8]:


#YRI - CEU 
import pandas as pd
dfe_list = ["K23","K17","Z21"]
model_list = ["3I21","3G09","2T12","4J17","10J19"]
for model_id in model_list:
    for dfe_id in dfe_list:
        cutoff = 1.447
        df = pd.read_csv(f"chr20_outputs/{dfe_id}_{model_id}freq_differences3.csv")
        #dataframe only includes frequency differences between YRI (or proxy) population and other population(s)
        plt.figure(figsize=(10, 6))
        plt.hist(df["(a-b)_diff"].abs().dropna(), bins=30, color="skyblue", edgecolor="black", log="true")
        plt.title(f"Histogram of Derived Allele Frequency Differences Between Populations for Model {model} under DFE {dfe}")
        plt.xlabel("Frequency Difference")
        plt.ylabel("Number of Mutations")
        plt.axvline(cutoff, color='red', linestyle='--', linewidth=2, label=f'cutoff = {cutoff}')
        plt.tight_layout()
        plt.show()
        print(f"Model {model_id} with DFE {dfe_id}")

        #Filtering data to only look at differences between YRI and CEU (or proxy) population
        if model_id == "9K19":
            df = df[(df["pop_b"] == 2)]
        else:
            df = df[df["pop_b"] == 1]
        num_extreme = (abs(df["(a-b)_diff"]) > cutoff).sum()
        total_mutations = len(df)
        print(num_extreme)
        print(total_mutations)
        p_value = num_extreme / total_mutations
        print(f"Empirical p-value (YRI (ancestral - alt) - CEU  (ancestral - alt) > Abs[1.447]): {p_value}")
        print()

    
#YRI - CHB (models 3I21 and 2T12 not included because they do not sample from a CHB or proxy population)
import pandas as pd
model_list = ["3G09","4J17","10J19","9K19"]
dfe_list = ["K23","K17","Z21"]
for model_id in model_list:
    for dfe_id in dfe_list:
        cutoff = 1.832
        df = pd.read_csv(f"chr20_outputs/{dfe_id}_{model_id}freq_differences3.csv")
        plt.figure(figsize=(10, 6))
        plt.hist(df["(a-b)_diff"].abs().dropna(), bins=30, color="skyblue", edgecolor="black", log="true")
        plt.title(f"Histogram of Derived Allele Frequency Differences Between Populations for Model {model} under DFE {dfe}")
        plt.xlabel("Frequency Difference")
        plt.ylabel("Number of Mutations")
        plt.tight_layout()
        plt.show()
        #dataframe only includes frequency differences between YRI (or proxy) population and other population(s)
        print(f"Model {model_id} with DFE {dfe_id}")
        
        
        #Filtering data to only look at differences between YRI and CHB (or proxy) population
        if model_id == "9K19":
            df = df[(df["pop_b"] == 5)]
        else:
            df = df[df["pop_b"] == 2]
        num_extreme = (abs(df["(a-b)_diff"]) > 1.832).sum()
        total_mutations = len(df)
        print(num_extreme)
        print(total_mutations)
        p_value = num_extreme / total_mutations
        print(f"Empirical p-value (YRI (ancestral - alt) - CHB (ancestral - alt) > Abs[1.832]): {p_value}")
        print()



# In[27]:


import pandas as pd
model_id = "3G09"
dfe_list = ["K23","K17","Z21"]
for dfe_id in dfe_list:
    #dataframe only includes frequency differences between YRI (or proxy) population and other population(s)
    print(f"Model {model_id} with DFE {dfe_id}")
    df = pd.read_csv(f"{dfe_id}_{model_id}freq_differences1.csv")
    #Filtering data to only look at differences between YRI and CHB (or proxy) population
    df = df[df["pop_b"] == 2]
    num_extreme = (abs(df["freq_diff"]) > 0.916).sum()
    total_mutations = len(df)
    p_value = num_extreme / total_mutations
    print(f"Empirical p-value (YRI - CHB > Abs[0.916]): {p_value}")
    print()


# In[48]:


import tskit
ts = tskit.load("K23_9K19full_chrom0.trees")
for k, pop in enumerate(ts.populations()):
    print(pop.metadata)
    print(
        f"The tree sequence has {len(ts.samples(k))} samples from "
        f"population {k}, which is {pop.metadata['name']}."
    )


# In[47]:


import tskit
ts = tskit.load("K23_9K19full_chrom0.trees")
for k, pop in enumerate(ts.populations()):
    print(pop.metadata)
    print(
        f"The tree sequence has {len(ts.samples(k))} samples from "
        f"population {k}, which is {pop.metadata['name']}."
    )


# In[ ]:




