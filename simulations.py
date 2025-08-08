#Current code used to generate csv file on compute cluster
import stdpopsim
import msprime
import tskit
import os
import pandas as pd
import itertools
output_dir = "all_simulations"
os.makedirs(output_dir, exist_ok=True)

#Simulate human evolution under all combinations of 6 models, 3 dfes, and 22 chromosomes (1-23) 
#Information on models and dfes found in stdpopsim catalog
model_list = {"3I21":"OutOfAfricaExtendedNeandertalAdmixturePulse_3I21",
              "3G09":"OutOfAfrica_3G09",
              "2T12":"OutOfAfrica_2T12", 
              "9K19":"AncientEurasia_9K19", 
              "10J19":"PapuansOutOfAfrica_10J19",
              "4J17":"OutOfAfrica_4J17"}
dfe_list = {"Z21":"GammaPos_Z21",
            "K23":"Mixed_K23",
            "K17":"Gamma_K17"}

#Get SLURM task ID
task_id = int(os.environ.get("SLURM_ARRAY_TASK_ID", 0))  # defaults to 0 if not running in array.


# Build a list of all combinations (model_id, dfe_id, chromosome)
all_tasks = list(itertools.product(model_list.items(), dfe_list.items(), range(1, 23)))

# Get corresponding parameters
(model_id, model), (dfe_id, dfe), chrom = all_tasks[task_id]
file_path = f"{output_dir}/{dfe_id}_{model_id}_chr{chrom}.trees"

#Check if simulation .trees file already exists. if it does, skip; if not, run the simulation. 
#This way, the script to be resubmitted to the cluster if some jobs time out  
if os.path.isfile(file_path) and os.path.getsize(file_path) > 0:
    print(f"Skipping {file_path}: already exists and is non-empty.")

else:
    engine = stdpopsim.get_engine("slim")
    species = stdpopsim.get_species("HomSap")
    model = species.get_demographic_model(model)
    contig = species.get_contig(f"chr{chrom}", mutation_rate=model.mutation_rate, genetic_map="HapMapII_GRCh38")
    
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

    #Initialize dfe, initialize exon annotation and exon intervals for given chromosome, then add dfe to contig defined above
    dfe = species.get_dfe(dfe)
    exons = species.get_annotations("ensembl_havana_104_exons")
    exon_intervals = exons.get_chromosome_annotations(f"chr{chrom}")
    contig.add_dfe(intervals=exon_intervals, DFE=dfe)

   #Run simulation  
    ts = engine.simulate(
        model,
        contig,
        samples,
        verbosity=1,
        slim_scaling_factor=1,
    )
    print(f"Simulating: Model={model_id}, DFE={dfe_id}, Chromosome=chr{chrom}")
    ts.dump(os.path.join(output_dir, f"{dfe_id}_{model_id}_chr{chrom}.trees"))
