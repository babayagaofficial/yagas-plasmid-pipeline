import sys
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

import glob
import os
import shutil
import pandas as pd
from skbio import TreeNode


def is_contained(clade, other_clade):
    #for identical plasmid content, keep the larger clade (ties broken by clade name so exactly one is kept)
    plasmids, other_plasmids = clades[clade]["plasmids"], clades[other_clade]["plasmids"]
    if plasmids==other_plasmids:
        return (clades[clade]["size"], clade) < (clades[other_clade]["size"], other_clade)
    return plasmids.issubset(other_plasmids)


typing = pd.read_csv(snakemake.input.typing, sep="\t")
chr_to_plasmid = pd.read_csv(snakemake.input.chr_to_plasmid, sep="\t")
rates = pd.read_csv(snakemake.input.rates)
host_tree = TreeNode.read(snakemake.input.tree, format="newick", convert_underscores=False)

fastafiles_list = [el[0] for el in pd.read_csv(snakemake.input.input_list, header=None).values]
fastafiles = {os.path.splitext(os.path.basename(el))[0]:el for el in fastafiles_list}

subcom = snakemake.params.cluster
min_rate_thresh = snakemake.params.min_rate
avg_rate_thresh = snakemake.params.avg_rate
min_plasmids = snakemake.params.min_plasmids

out_dir = snakemake.output.out_dir
for sub_dir in ["correct_clades", "trees", "lists"]:
    os.makedirs(os.path.join(out_dir, sub_dir), exist_ok=True)

plasmid_to_host = dict(zip(chr_to_plasmid["plasmid"], chr_to_plasmid["chr"]))
plasmids = typing[typing["type"]==subcom]["plasmid"].to_list()

#find clades which pass the rate thresholds
clades = {}
for file in glob.glob(os.path.join(snakemake.params.phylofactor_dir, "clades", "*.txt")):
    clade = os.path.basename(file).replace(".txt", "")
    with open(file, "r") as f:
        genomes = set(line.strip() for line in f if line.strip())

    clade_plasmids = set()
    clade_genomes = set()
    for plasmid in plasmids:
        genome = plasmid_to_host[plasmid]
        if genome in genomes:
            clade_plasmids.add(plasmid)
            clade_genomes.add(genome)

    if len(clade_genomes)==0:
        continue

    clade_rates = rates[rates["Genome_ID"].isin(clade_genomes)]["rate"]
    avg_rate = clade_rates.sum()/len(clade_genomes)
    min_rate = clade_rates.min()

    if min_rate>=min_rate_thresh and avg_rate>=avg_rate_thresh and len(clade_plasmids)>=min_plasmids:
        clades[clade] = {"file": file, "plasmids": clade_plasmids, "genomes": clade_genomes, "size": len(genomes)}

#keep only clades whose plasmids are not contained in another passing clade
for clade, members in clades.items():
    if any(is_contained(clade, other_clade) for other_clade in clades if other_clade!=clade):
        continue

    shutil.copy(members["file"], os.path.join(out_dir, "correct_clades"))

    tree = host_tree.shear(members["genomes"])
    tree.write(os.path.join(out_dir, "trees", f"cl{clade}_{subcom}.nw"), format="newick")

    leaves = [tip.name for tip in tree.tips()]
    ordered_plasmids = sorted(members["plasmids"], key=lambda plasmid: leaves.index(plasmid_to_host[plasmid]))
    with open(os.path.join(out_dir, "lists", f"cl{clade}_{subcom}.txt"), "w") as f:
        for plasmid in ordered_plasmids:
            f.write(fastafiles[plasmid] + "\n")
