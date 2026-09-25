import os
import glob
import shutil
import subprocess

configfile: "config.yaml"

wildcard_constraints:
    cluster="[^/]+"

PLING_DIR = config["output_dir"] + "/pling_d" + str(config["dcj-indel"]) + "_c" + str(config["containment"]).replace(".", "")

def get_multifasta():
    if config["input_list"]:
        return config["output_dir"] + "/all_plasmids.fna"
    else:
        return config["multifasta"]

def get_input_list():
    if config["input_list"]:
        return config["input_list"]
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


rule all:
    input:
        ggcallaroo = lambda wildcards: [config["output_dir"] + f"/ggcallaroo/{cluster}" for cluster in get_clusters()],
        pangraph = lambda wildcards: [config["output_dir"] + f"/pangraph/{cluster}" for cluster in get_clusters()],
        median_hist = config["output_dir"] + "/boundary/dcj_median.png", mean_hist = config["output_dir"] + "/boundary/dcj_mean.png",   # histograms of boundary values
        median_box = config["output_dir"] + "/boundary/internal_vs_boundary_median.png",  mean_box = config["output_dir"] + "/boundary/internal_vs_boundary_mean.png",    # swarm + box plots, internal vs boundary
        dcj_avg_tsv = config["output_dir"] + "/boundary/dcj_averages.tsv",
        cluster_specs_tsv = config["output_dir"] + "/cluster_specs.tsv",
        plot = config["output_dir"] + "/dcj_distr/hist_plot.png",
        stats = config["output_dir"] + "/dcj_distr/stats.txt",
        parsnp_dir = lambda wildcards: [config["output_dir"] + f"/parsnp/{cluster}" for cluster in get_clusters()],
        dcj_trees = PLING_DIR + "/submatrices"

rule separate_fastas:
    input:
        multi = config["multifasta"]
    output:
        fasta_dir = directory(config["output_dir"] + "/fastas"),
        fasta_list = config["output_dir"] + "/plasmid_list.txt"
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

rule cat_fastas:
    input:
        fastas = config["input_list"]
    output:
        multifasta = config["output_dir"] + "/all_plasmids.fna"
    run:
        from Bio import SeqIO

        with open(output.multifasta, "w") as multi:
            with open(input.fastas, "r") as f:
                for line in f:
                    path = line.strip()
                    record = SeqIO.read(path, "fasta")
                    SeqIO.write(multi, record)

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
        pling_out = PLING_DIR
    conda:
        "pling"
    resources:
        pass
    threads: config["pling_threads"]
    shell:
        "pling cluster align {input.fastas} {params.pling_out} --cores {threads} --dcj {params.dcj} --containment_distance {params.containment}"

rule mobtyper:
    input:
        fastas = get_multifasta()
    output:
        mob = config["output_dir"] + "/mobtyper_results.txt"
    conda:
        "mobsuite"
    resources:
        pass
    shell:
        "mob_typer --multi --infile {input.fastas} --out_file {output.mob}"

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

rule ggcallaroo:
    input:
        fasta_list = lambda wildcards: get_cluster_list(wildcards.cluster)
    output:
        ann_dir = directory(config["output_dir"] + "/ggcallaroo/{cluster}")
    conda:
        "ggcallaroo"
    resources:
        pass
    threads: 8
    params:
        ggcallaroo_path = config["ggcallaroo"],
        bakta_db = config["bakta_db"],
        ggcaller_cli_args = "",
        panaroo_cli_args = ""
    shell:
        """
        snakemake --cores {threads} --use-conda --snakefile {params.ggcallaroo_path}/Snakefile --directory {params.ggcallaroo_path} --config refs={input.fasta_list} output_dir={output.ann_dir} bakta_db={params.bakta_db}
        """

rule pangraph:
    input:
        fastas = lambda wildcards: get_list(wildcards.cluster)
    output:
        ann_dir = directory(config["output_dir"] + "/pangraph/{cluster}")
    resources:
        mem_mb=lambda wildcards, attempt: 40000*attempt
    threads: 8
    shell:
        """
        mkdir -p {output.ann_dir}
        pangraph build --circular -k minimap2 -s 20 -b 5 --len 200 {input.fastas} > {output.ann_dir}/pangraph.json
        pangraph export gfa --output {output.ann_dir}/pangraph.gfa --minimum-length 200 {output.ann_dir}/pangraph.json
        """

rule rel_core_sizes:
    input:
        ggcallaroo_dirs = lambda wildcards: [config["output_dir"] + f"/ggcallaroo/{cluster}" for cluster in get_clusters()],
        list_dir = config["output_dir"] + "/cluster_lists",
        mob = config["output_dir"] + "/mobtyper_results.txt"
    output:
        plot = config["output_dir"] + "/rel_core/rel_core_plot.png",
        tsv = config["output_dir"] + "/rel_core/rel_core.tsv"
    params:
        ggcaller_dir = config["output_dir"] + "/ggcallaroo"
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
    conda: "phylofactor"
    resources:
        mem_mb=lambda wildcards, attempt: 20000*attempt
    shell:
        "R < scripts/phylofactor.R {input.tree} {input.traits} {params.cluster} {params.out_dir} --no-save"

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
        "pling"
    resources:
        pass
    shell:
        "pling submatrix {params.pling_out} --vis_trees"

rule parsnp:
    input:
        fastas = lambda wildcards: get_list(wildcards.cluster)
    output:
        parsnp_dir = directory(config["output_dir"] + "/parsnp/{cluster}")
    params:
        cluster = lambda wildcards: wildcards.cluster
    resources:
        mem_mb=lambda wildcards, attempt: 40000*attempt
    threads: config["parsnp_threads"]
    shadow: "shallow"
    run:
        os.mkdir(params.cluster)
        for file in input.fastas:
            shutil.copy(file, params.cluster)
        try:
            subprocess.run(f"parsnp -c {params.cluster} -p {threads} -o {output.parsnp_dir}", shell=True, check=True, capture_output=True)
        except subprocess.CalledProcessError as e:
            print(e.stderr.decode())
            print(e)
            raise e


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
    script:
        "scripts/boundary.py"

