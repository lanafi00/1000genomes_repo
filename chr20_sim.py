#!/usr/bin/env python
# coding: utf-8

# In[3]:


#Current code used to generate csv file
import stdpopsim
import msprime
import tskit
import os
import pandas as pd
import itertools

model_list = {"3I21":"OutOfAfricaExtendedNeandertalAdmixturePulse_3I21",
              "3G09":"OutOfAfrica_3G09",
              "2T12":"OutOfAfrica_2T12", 
              "9K19":"AncientEurasia_9K19", 
              "10J19":"PapuansOutOfAfrica_10J19",
              "4J17":"OutOfAfrica_4J17"}
dfe_list = {"Z21":"GammaPos_Z21",
            "K23":"Mixed_K23",
            "K17":"Gamma_K17"}

for model_id, model in model_list.items():
    for dfe_id, dfe in dfe_lift.items():
        task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))  # defaults to 0 if not running in array 
        engine = stdpopsim.get_engine("slim")
        species = stdpopsim.get_species("HomSap")
        model = species.get_demographic_model(model)
        model_id = model_id
        dfe_id = dfe_id
        contig = species.get_contig("chr20", mutation_rate=model.mutation_rate, genetic_map="HapMapII_GRCh38")

        #Must change samples based on model used
        if model_id == "3I21":
            samples = {"YRI":108,"CEU":99}
        elif model_id in ["3G09", "10J19", "4J17"]:
            samples = {"YRI":108,"CEU":99,"CHB":103}
        elif model_id == "2T12":
            samples = {"AFR":108,"EUR":99}
        elif model_id == "9K19":
            samples = {"Han":103,"Mbuti":108,"Sardinian":99}
        else:
            raise ValueError(f"Unknown model_id: {model_id}")

        dfe = species.get_dfe(dfe)
        exons = species.get_annotations("ensembl_havana_104_exons")
        exon_intervals = exons.get_chromosome_annotations("chr20")
        contig.add_dfe(intervals=exon_intervals, DFE=dfe)

        ts = engine.simulate(
            model,
            contig,
            samples,
            verbosity=1,
            slim_scaling_factor=1,
        )

        ts.dump(f"{dfe_id}_{model_id}full_chrom{task_id}.trees")


# In[2]:


import pandas as pd
import matplotlib.pyplot as plt
dfe_list = ["K23","K17","Z21"]
model_list = ["3I21","3G09","2T12","9K19","4J17"]
for model in model_list:
    for dfe in dfe_list:

    # Load the CSV
        df = pd.read_csv(f"outputs/{dfe}_{model}freq_differences3.csv")

    # Plot histogram
        plt.figure(figsize=(10, 6))
        plt.hist(df["freq_diff"].abs().dropna(), bins=30, color="skyblue", edgecolor="black", log="true")
        plt.title(f"Histogram of Derived Allele Frequency Differences Between Populations for Model {model} under DFE {dfe}")
        plt.xlabel("Frequency Difference")
        plt.ylabel("Number of Mutations")
        plt.tight_layout()
        plt.show()


# In[ ]:




