import os
import glob

configfile: "config.yaml"

wildcard_constraints:
    cluster="[^/]+"

#{cluster} in ggCallaroo's output_dir turns each of its rules into a per-cluster rule
GGCALLAROO_DIR = config["output_dir"] + "/ggcallaroo/{cluster}"
GGCALLAROO_OUTPUTS = ["pan_genome_reference.fa", "combined_DNA_CDS.fasta", "combined_protein_CDS.fasta", "gene_data.csv",
                      "gene_presence_absence.csv", "gene_presence_absence_roary.csv", "final_graph.gml", "pre_filt_graph.gml"]

LOG_DIR = config["output_dir"] + "/logs"

PLING_DIR = config["output_dir"] + "/pling_d" + str(config["dcj-indel"]) + "_c" + str(config["containment"]).replace(".", "")

def get_rule_config(rule):
    #rule-specific values override the defaults
    return {**config["resources"]["default"], **config["resources"].get(rule, {})}

def get_resources(rule):
    #memory and time scale with the retry attempt; threads are set separately with get_threads
    res = {key: value for key, value in get_rule_config(rule).items() if key!="threads"}
    return {key: (lambda wildcards, attempt, value=value: value*attempt) for key, value in res.items()}

def get_threads(rule):
    return get_rule_config(rule).get("threads", 1)

#input is given either as a list of fasta files or as one multifasta; the other is created from it
INPUT_LIST = config.get("input_list")
MULTIFASTA = config.get("multifasta")
if bool(INPUT_LIST)==bool(MULTIFASTA):
    raise ValueError("Set exactly one of input_list or multifasta in config.yaml")

def get_multifasta():
    if INPUT_LIST:
        return config["output_dir"] + "/all_plasmids.fna"
    else:
        return MULTIFASTA

def get_input_list():
    if INPUT_LIST:
        return INPUT_LIST
    else:
        return config["output_dir"] + "/plasmid_list.txt"

def get_cluster_list(cluster):
    #calling the checkpoint makes snakemake wait for cluster_lists before evaluating rules which depend on it
    return checkpoints.cluster_lists.get().output.list_dir + f"/{cluster}.txt"

def get_clusters():
    list_dir = checkpoints.cluster_lists.get().output.list_dir
    return [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{list_dir}/*.txt")]

def get_list(cluster):
    files = []
    with open(get_cluster_list(cluster)) as f:
        for line in f:
            if line.strip():
                files.append(line.strip())
    return files

def get_sourmash():
    if config["sourmash"]:
        return "--sourmash"
    else:
        ""

#lightweight rules which run on the submitting node instead of being sent to slurm
localrules: all, cluster_lists, sc_in_chr, dcj_distr, cluster_specs

rule all:
    input:
        ggcallaroo = lambda wildcards: [GGCALLAROO_DIR.format(cluster=cluster) + "/annotated/" + file for cluster in get_clusters() for file in GGCALLAROO_OUTPUTS],
        pangraph = lambda wildcards: [config["output_dir"] + f"/pangraph/{cluster}" for cluster in get_clusters()],
        median_hist = config["output_dir"] + "/boundary/dcj_median.png", mean_hist = config["output_dir"] + "/boundary/dcj_mean.png",   # histograms of boundary values
        median_box = config["output_dir"] + "/boundary/internal_vs_boundary_median.png",  mean_box = config["output_dir"] + "/boundary/internal_vs_boundary_mean.png",    # swarm + box plots, internal vs boundary
        dcj_avg_tsv = config["output_dir"] + "/boundary/dcj_averages.tsv",
        cluster_specs_tsv = config["output_dir"] + "/cluster_specs.tsv",
        plot = config["output_dir"] + "/dcj_distr/hist_plot.png",
        stats = config["output_dir"] + "/dcj_distr/stats.txt",
        parsnp_dir = lambda wildcards: [config["output_dir"] + f"/parsnp/{cluster}" for cluster in get_clusters()],
        dcj_trees = PLING_DIR + "/submatrices",
        rel_core_plot = config["output_dir"] + "/rel_core/rel_core_plot.png",
        rel_core_tsv = config["output_dir"] + "/rel_core/rel_core.tsv",
        post_phylofactor = lambda wildcards: [config["output_dir"] + f"/post_phylofactor/{cluster}" for cluster in get_clusters()]

if MULTIFASTA:
    rule separate_fastas:
        input:
            multi = MULTIFASTA
        output:
            fasta_dir = directory(config["output_dir"] + "/fastas"),
            fasta_list = config["output_dir"] + "/plasmid_list.txt"
        resources: **get_resources("separate_fastas")
        threads: get_threads("separate_fastas")
        run:
            from Bio import SeqIO
            from Bio.SeqRecord import SeqRecord

            os.makedirs(output.fasta_dir, exist_ok=True)
            with open(output.fasta_list, "w") as fasta_list:
                for record in SeqIO.parse(input.multi, "fasta"):
                    sep_record = SeqRecord(record.seq, record.id, "")
                    with open(f"{output.fasta_dir}/{record.id}.fna", "w") as output_handle:
                        SeqIO.write(sep_record, output_handle, "fasta")
                    fasta_list.write(f"{output.fasta_dir}/{record.id}.fna\n")

if INPUT_LIST:
    rule cat_fastas:
        input:
            fastas = INPUT_LIST
        output:
            multifasta = config["output_dir"] + "/all_plasmids.fna"
        resources: **get_resources("cat_fastas")
        threads: get_threads("cat_fastas")
        run:
            from Bio import SeqIO

            with open(output.multifasta, "w") as multi:
                with open(input.fastas, "r") as f:
                    for line in f:
                        path = line.strip()
                        record = SeqIO.read(path, "fasta")
                        SeqIO.write(record, multi, "fasta")

rule pling:
    input:
        fastas = get_input_list()
    output:
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        hubs = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/hub_plasmids.csv",
        communities_pickle = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/communities.pkl",
        dcjs = PLING_DIR + "/all_plasmids_distances.tsv"
    params:
        dcj = int(config["dcj-indel"]),
        containment = float(config["containment"]),
        sourmash = get_sourmash(),
        batch_size = config[batch_size],
        pling_out = PLING_DIR
    conda:
        "envs/pling.yaml"
    resources: **get_resources("pling")
    threads: get_threads("pling")
    log:
        LOG_DIR + "/pling.log"
    shell:
        "pling cluster align {input.fastas} {params.pling_out} --cores {threads} --dcj {params.dcj} --containment_distance {params.containment} {params.sourmash} {params.batch_size}> {log} 2>&1"

rule mobtyper:
    input:
        fastas = get_multifasta()
    output:
        mob = config["output_dir"] + "/mobtyper_results.txt"
    conda:
        "envs/mobsuite.yaml"
    resources: **get_resources("mobtyper")
    threads: get_threads("mobtyper")
    log:
        LOG_DIR + "/mobtyper.log"
    shell:
        "mob_typer --multi --infile {input.fastas} --out_file {output.mob} > {log} 2>&1"

checkpoint cluster_lists:
    input:
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        input_list = get_input_list()
    output:
        list_dir = directory(config["output_dir"] + "/cluster_lists")
    params:
        min_cluster_size = config["big_subcomm_size"]
    run:
        import pandas as pd
        import os

        clusters_df = pd.read_csv(input.typing, sep="\t")
        clusters = list(set(clusters_df["type"].values))

        fastafiles_list = [el[0] for el in pd.read_csv(input.input_list, header=None).values]
        fastafiles = {os.path.splitext(os.path.basename(el))[0]:el for el in fastafiles_list}

        os.makedirs(output.list_dir, exist_ok=True)
        for cluster in clusters:
            if len(clusters_df[clusters_df["type"]==cluster])>params.min_cluster_size:
                with open(f"{output.list_dir}/{cluster}.txt", "w") as f:
                    for name in clusters_df[clusters_df["type"]==cluster]["plasmid"].values:
                        f.write(fastafiles[name] + "\n")

GGCALLAROO_CONFIG = {
    "output_dir": GGCALLAROO_DIR,
    "refs": lambda wildcards: get_cluster_list(wildcards.cluster),
    "reads": None,
    "ggcaller_cli_args": config.get("ggcaller_cli_args", "--save"),
    "panaroo_cli_args": config.get("panaroo_cli_args", "--clean-mode moderate"),
    "bakta_db": config["bakta_db"]
}

module ggcallaroo:
    #pinned to a commit, since the rules used below are referred to by name
    snakefile: github("samhorsfield96/ggCallaroo", path="Snakefile", commit="4ffa1068bc6fd06f90abcab4a7e1d6277e45ba9b")
    config: GGCALLAROO_CONFIG

#ggCallaroo's own rule all is excluded, since its inputs contain the {cluster} wildcard
use rule translate_representatives, annotate_pan_ref, annotate_dna_CDS, annotate_dna_prot, annotate_gene_data, annotate_gpa, annotate_gpa_roary, annotate_gml, annotate_gml_pref_filt from ggcallaroo as ggcallaroo_* with:
    resources: **get_resources("default")
    threads: get_threads("default")

use rule ggcaller from ggcallaroo as ggcallaroo_ggcaller with:
    resources: **get_resources("ggcallaroo_ggcaller")
    threads: get_threads("ggcallaroo_ggcaller")

use rule panaroo from ggcallaroo as ggcallaroo_panaroo with:
    resources: **get_resources("ggcallaroo_panaroo")
    threads: get_threads("ggcallaroo_panaroo")

use rule bakta_proteins from ggcallaroo as ggcallaroo_bakta_proteins with:
    resources: **get_resources("ggcallaroo_bakta_proteins")
    threads: get_threads("ggcallaroo_bakta_proteins")

rule pangraph:
    input:
        fastas = lambda wildcards: get_list(wildcards.cluster)
    output:
        ann_dir = directory(config["output_dir"] + "/pangraph/{cluster}")
    resources: **get_resources("pangraph")
    threads: get_threads("pangraph")
    log:
        LOG_DIR + "/pangraph/{cluster}.log"
    shell:
        """
        mkdir -p {output.ann_dir}
        pangraph build --circular -k minimap2 -s 20 -b 5 --len 200 -j {threads} {input.fastas} > {output.ann_dir}/pangraph.json 2> {log}
        pangraph export gfa --output {output.ann_dir}/pangraph.gfa --minimum-length 200 -j {threads} {output.ann_dir}/pangraph.json >> {log} 2>&1
        """

rule rel_core_sizes:
    input:
        panaroo = lambda wildcards: [GGCALLAROO_DIR.format(cluster=cluster) + "/panaroo/pan_genome_reference.fa" for cluster in get_clusters()],
        list_dir = config["output_dir"] + "/cluster_lists",
        mob = config["output_dir"] + "/mobtyper_results.txt"
    output:
        plot = config["output_dir"] + "/rel_core/rel_core_plot.png",
        tsv = config["output_dir"] + "/rel_core/rel_core.tsv"
    params:
        ggcaller_dir = config["output_dir"] + "/ggcallaroo"
    log:
        LOG_DIR + "/rel_core_sizes.log"
    resources: **get_resources("rel_core_sizes")
    threads: get_threads("rel_core_sizes")
    conda:
        "envs/python.yaml"
    script:
        "scripts/get_core_sizes.py"

rule sc_in_chr:
    input:
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        chr_to_plasmid = config["chr_to_plasmid"]
    output:
        plasmid_presence_absence = config["output_dir"] + "/host_presence/presence_per_host.tsv"
    params:
        big = config["big_subcomm_size"]
    run:
        import pandas as pd
        typing = pd.read_csv(input.typing, sep="\t")
        plasmid_presence = pd.read_csv(input.chr_to_plasmid, sep="\t")
        subcomm_sizes = typing["type"].value_counts()
        big_subcomms = subcomm_sizes[subcomm_sizes>params.big].index
        hosts = plasmid_presence["chr"].unique()

        #1 if a host carries at least one plasmid of the subcommunity, 0 otherwise
        merged = plasmid_presence.merge(typing[typing["type"].isin(big_subcomms)], on="plasmid")
        presence_df = pd.crosstab(merged["chr"], merged["type"]).clip(upper=1)
        presence_df = presence_df.reindex(index=hosts, columns=big_subcomms, fill_value=0).rename_axis(index=None, columns=None)
        presence_df.sort_index(inplace=True)
        presence_df.sort_index(axis=1, inplace=True)
        presence_df.to_csv(output.plasmid_presence_absence, sep="\t")

rule phylofactor:
    input:
        tree = config["host_tree"],
        traits = config["output_dir"] + "/host_presence/presence_per_host.tsv"
    output:
        tree_vis = config["output_dir"] + "/phylofactor/{cluster}/tree.pdf",
        rates = config["output_dir"] + "/phylofactor/{cluster}/rates.csv"
    params:
        cluster = lambda wildcards: wildcards.cluster,
        out_dir = config["output_dir"] + "/phylofactor/{cluster}"
    conda:
        "envs/phylofactor.yaml"
    resources: **get_resources("phylofactor")
    threads: get_threads("phylofactor")
    log:
        LOG_DIR + "/phylofactor/{cluster}.log"
    shell:
        "R < scripts/phylofactor.R {input.tree} {input.traits} {params.cluster} {params.out_dir} --no-save > {log} 2>&1"

rule post_phylofactor:
    input:
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        chr_to_plasmid = config["chr_to_plasmid"],
        tree = config["host_tree"],
        input_list = get_input_list(),
        rates = config["output_dir"] + "/phylofactor/{cluster}/rates.csv"
    output:
        out_dir = directory(config["output_dir"] + "/post_phylofactor/{cluster}")
    params:
        cluster = lambda wildcards: wildcards.cluster,
        phylofactor_dir = config["output_dir"] + "/phylofactor/{cluster}",
        min_rate = 0.4,
        avg_rate = 0.5,
        min_plasmids = 4
    log:
        LOG_DIR + "/post_phylofactor/{cluster}.log"
    resources: **get_resources("post_phylofactor")
    threads: get_threads("post_phylofactor")
    conda:
        "envs/python.yaml"
    script:
        "scripts/filter_phylofactor.py"

rule dcj_trees:
    input:
        #depend on pling's declared outputs, since no rule outputs PLING_DIR itself
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        dcjs = PLING_DIR + "/all_plasmids_distances.tsv"
    output:
        submatrices_dir = directory(PLING_DIR + "/submatrices")
    params:
        pling_out = PLING_DIR
    conda:
        "envs/pling.yaml"
    resources: **get_resources("dcj_trees")
    threads: get_threads("dcj_trees")
    log:
        LOG_DIR + "/dcj_trees.log"
    shell:
        "pling submatrix {params.pling_out} --vis_trees > {log} 2>&1"

rule parsnp:
    input:
        fastas = lambda wildcards: get_list(wildcards.cluster)
    output:
        parsnp_dir = directory(config["output_dir"] + "/parsnp/{cluster}")
    params:
        cluster = lambda wildcards: wildcards.cluster
    resources: **get_resources("parsnp")
    threads: get_threads("parsnp")
    shadow: "shallow"
    log:
        LOG_DIR + "/parsnp/{cluster}.log"
    conda:
        "envs/parsnp.yaml"
    shell:
        """
        mkdir {params.cluster}
        cp {input.fastas} {params.cluster}/
        parsnp -d {params.cluster} -c -p {threads} -o {output.parsnp_dir} > {log} 2>&1
        """


rule dcj_distr:
    input:
        dcj_dists = PLING_DIR + "/all_plasmids_distances.tsv"
    output:
        plot = config["output_dir"] + "/dcj_distr/hist_plot.png",
        stats = config["output_dir"] + "/dcj_distr/stats.txt"
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        dists = pd.read_csv(input.dcj_dists, sep="\t")
        fig, ax = plt.subplots()
        sns.histplot(data=dists,x="distance", ax=ax, discrete=True)
        plt.savefig(output.plot)

        with open(output.stats, "w") as f:
            f.write("mode:"+str(dists["distance"].mode().tolist())+"\n")
            f.write("median:"+str(dists["distance"].median())+"\n")
            f.write("mean:"+str(dists["distance"].mean())+"\n")

rule cluster_specs:
    input:
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        mob = config["output_dir"] + "/mobtyper_results.txt"
    output:
        tsv = config["output_dir"] + "/cluster_specs.tsv"
    log:
        LOG_DIR + "/cluster_specs.log"
    conda:
        "envs/python.yaml"
    script:
        "scripts/cluster_specs.py"

rule boundary:
    input:
        typing = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/typing.tsv",
        hubs = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/hub_plasmids.csv",
        communities_pickle = PLING_DIR + "/dcj_thresh_" + str(config["dcj-indel"]) + "_graph/objects/communities.pkl",
        dcjs = PLING_DIR + "/all_plasmids_distances.tsv"
    params:
        big_subcomm_size = config["big_subcomm_size"]
    output:
        median_hist = config["output_dir"] + "/boundary/dcj_median.png", mean_hist = config["output_dir"] + "/boundary/dcj_mean.png",   # histograms of boundary values
        median_box = config["output_dir"] + "/boundary/internal_vs_boundary_median.png",  mean_box = config["output_dir"] + "/boundary/internal_vs_boundary_mean.png",    # swarm + box plots, internal vs boundary
        tsv = config["output_dir"] + "/boundary/dcj_averages.tsv"
    log:
        LOG_DIR + "/boundary.log"
    resources: **get_resources("boundary")
    threads: get_threads("boundary")
    conda:
        "envs/python.yaml"
    script:
        "scripts/boundary.py"

