import sys
sys.stdout = sys.stderr = open(snakemake.log[0], "w")

from typing import cast
from plasnet.communities import Communities
import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import seaborn as sns
from statistics import mean, median
import itertools

def get_dists(dist_filepath):
    dcj={}
    with open(dist_filepath, "r") as f:
        next(f)
        for line in f:
            plasmid_1, plasmid_2, dist = line.strip().split('\t')
            dcj[frozenset((plasmid_1,plasmid_2))] = int(dist)
    return dcj

communities = cast(Communities, Communities.load(snakemake.input.communities_pickle))

typing = pd.read_csv(snakemake.input.typing, sep="\t")
hubs = set(pd.read_csv(snakemake.input.hubs, sep="\t")["hub_plasmids"].to_list())

dists = get_dists(snakemake.input.dcjs)

big_subcomm_size = snakemake.params.big_subcomm_size

#group big subcommunities by the number of the community they belong to
subcomm_plasmids = typing.groupby("type")["plasmid"].apply(list)
subcomms_per_community = {}
for com, plasmids in subcomm_plasmids.items():
    if len(plasmids)>big_subcomm_size:
        subcomms_per_community.setdefault(com.split("_")[1], []).append(com)

rows = []
for community in communities:
    num=community.label.split("_")[1]
    for com in subcomms_per_community.get(num, []):
        plasmids = subcomm_plasmids[com]

        boundary_dcj = [dists[frozenset(edge)] for edge in nx.edge_boundary(community, plasmids) if edge[1] not in hubs]
        if len(boundary_dcj)>0:
            rows.append({"subcommunity":com, "median":median(boundary_dcj), "mean":mean(boundary_dcj), "location":"boundary"})

        internal_dcj = [dists[frozenset(pair)] for pair in itertools.combinations(plasmids,2)]
        rows.append({"subcommunity":com, "median":median(internal_dcj), "mean":mean(internal_dcj), "location":"internal"})

results = pd.DataFrame(rows, columns=["subcommunity", "median", "mean", "location"])
boundary = results[results["location"]=="boundary"]

for stat in ["median", "mean"]:
    fig, ax = plt.subplots()
    sns.histplot(data=boundary, x=stat, ax=ax, discrete=True)
    plt.savefig(snakemake.output[f"{stat}_hist"])
    plt.close(fig)

    fig, ax = plt.subplots()
    sns.swarmplot(data=results, y=stat, x="location", zorder=0, ax=ax)
    sns.boxplot(data=results, y=stat, x="location", fill=False, ax=ax)
    plt.savefig(snakemake.output[f"{stat}_box"])
    plt.close(fig)

results.to_csv(snakemake.output.tsv, sep="\t", index=False)
