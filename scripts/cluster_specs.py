import pandas as pd
from statistics import mean

typing = pd.read_csv(snakemake.input.typing, sep="\t")
mob = pd.read_csv(snakemake.input.mob, sep="\t")

mob_per_subcomm = {"type":[], "#plasmids":[], "avg_length":[], "rep_types":[], "relaxase_types":[], "mpf_type":[], "predicted_mobility":[]}

subcomms = list(set(typing["type"].to_list()))

for subcomm in subcomms:
    mob_per_subcomm["type"].append(subcomm)
    plasmids = typing[typing["type"]==subcomm]["plasmid"].to_list()

    mob_per_subcomm["#plasmids"].append(len(plasmids))
    mob_per_subcomm["avg_length"].append(round(mean(mob[mob["sample_id"].isin(plasmids)]["size"].to_list())))


    reps = list(set(mob[mob["sample_id"].isin(plasmids)]["rep_type(s)"].to_list()))
    reps = [("|").join(sorted(rep.split(","))) for rep in reps]
    relaxases = list(set(mob[mob["sample_id"].isin(plasmids)]["relaxase_type(s)"].to_list()))
    mpf = list(set(mob[mob["sample_id"].isin(plasmids)]["mpf_type"].to_list()))
    mobility = list(set(mob[mob["sample_id"].isin(plasmids)]["predicted_mobility"].to_list()))
    mob_per_subcomm["rep_types"].append(",".join(reps))
    mob_per_subcomm["relaxase_types"].append(",".join(relaxases))
    mob_per_subcomm["mpf_type"].append(",".join(mpf))
    mob_per_subcomm["predicted_mobility"].append(",".join(mobility))

mob_per_subcomm_df = pd.DataFrame(data=mob_per_subcomm)
mob_per_subcomm_df.sort_values(by=["#plasmids"], inplace=True, ascending=False)
mob_per_subcomm_df.to_csv(snakemake.output.tsv, sep="\t", index=False)

