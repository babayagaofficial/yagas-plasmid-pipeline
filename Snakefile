import os
import glob

configfile: "config.yaml"

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

def get_clusters():
    clusterpath = "blub"
    return [os.path.basename(el).replace('.txt','') for el in glob.glob(f"{cluster_path}/*.txt")]

def get_list(cluster):
    files = []
    cluster_path = "blub"
    with open(f"{cluster_path}/{cluster}.txt") as f:
        for line in f:
            files.append(line.strip())
    return files


rule all:
    input:
        ggcaller = [config["output_dir"] + f"/ggcaller/{cluster}" for cluster in get_clusters()]

rule separate_fastas:
    input:
        multi = config["multifasta"]
    output:
        fasta_dir = directory(config["output_dir"] + "fastas")
        fasta_list = config["output_dir"] + "/plasmid_list.txt"
    run:
        from Bio import SeqIO
        from Bio.SeqRecord import SeqRecord

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

checkpoint pling:
    input:
        fastas = get_input_list()
    output:
        pling_out = config["output_dir"] + "/pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '')
        typing = config["output_dir"] + "/pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '') + "/dcj_thresh_" + config["dcj-indel"] + "_graph/objects/typing.tsv"
    params:
        dcj = int(config["dcj-indel"])
        containment = float(config["containment"])
    conda:
        "pling"
    resources:
        pass
    threads: config["pling_threads"]
    shell:
        "pling cluster align {input.fastas} {output.pling_out} --cores {threads} --dcj {params.dcj} --containment_distance {params.containment}"

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

rule cluster_lists:
    input:
        typing = config["output_dir"] + "/pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '') + "/dcj_thresh_" + config["dcj-indel"] + "_graph/objects/typing.tsv"
        input_list = get_input_list
    output:
        list_dir = config["output_dir"] + "/cluster_lists"
    params:
        min_cluster_size = config["min_cluster_size"]
    run:
        import pandas as pd
        import os

        clusters_df = pd.read_csv(input.typing, sep="\t")
        clusters = list(set(clusters_df["type"].values))

        fastafiles_list = [el[0] for el in pd.read_csv(input.input_list, header=None).values]
        fastafiles = {os.path.splitext(os.path.basename(el))[0]:el for el in fastafiles_list}

        for cluster in clusters:
            if len(clusters_df[clusters_df["type"]==cluster])>params.min_cluster_size:
                with open(f"{output.list_dir}/{cluster}.txt", "w") as f:
                    for name in clusters_df[clusters_df["type"]==cluster]["plasmid"].values:
                        f.write(fastafiles[name] + "\n")

rule ggcallaroo:
    input:
        fasta_list = lambda wildcards: config["output_dir"] + "/cluster_lists/" + wildcards.cluster + ".txt"
    output:
        ann_dir = directory(config["output_dir"] + "/ggcallaroo/{cluster}")
    conda:
        "ggcallaroo"
    resources:
        pass
    threads: 8
    params:
        gcallaroo_path = config["ggcallaroo"],
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
        ann_dir = directory(f"{out}/{{cluster}}")
    conda:
        "pangraph"
    resources:
        mem_mb=lambda wildcards, attempt: 40000*attempt
    threads: 8
    shell:
        """
        export JULIA_NUM_THREADS={threads}
        julia --project=. src/PanGraph.jl build --circular -k minimap2 -s 20 -b 5 {input.fastas} > {output}/pangraph.json
        julia --project=. src/PanGraph.jl export --no-duplications --output-directory {output} pangraph.json
        """

rule rel_core_sizes:
    input:
        ggcaller_dir = config["output_dir"] + "/ggcaller/{cluster}" #this changes with ggcallaroo
        list_dir = config["output_dir"] + "/cluster_lists"
        mob = config["output_dir"] + "/mobtyper_results.txt"
    output:
        plot = config["output_dir"] + "rel_core/rel_core_plot.png"
        tsv = config["output_dir"] + "rel_core/rel_core.tsv"
    script:
        "scripts/get_core_sizes.py"

rule sc_in_chr:
    input:
        typing = config["output_dir"] + "pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '') + "/dcj_thresh_" + config["dcj-indel"] + "_graph/objects/typing.tsv"
        chr_to_plasmid = config["chr_to_plasmid"]
    output:
        plasmid_presence_absence = config["output_dir"] + "/host_presence/presence_per_host.tsv"
    params:
        big = config["big_subcomm_sizes"]
    run:
        import pandas as pd
        typing = pd.read_csv(input.typing, sep="\t")
        plasmid_presence = pd.read_csv(input.chr_to_plasmid, sep="\t")
        subcomms = list(set(typing["type"].to_list()))
        big_subcomms = [subcomm for subcomm in subcomms if len(typing[typing["type"]==subcomm])>params.big]
        subcomm_presence = {subcomm:[] for subcomm in big_subcomms}
        names=[]
        for host in plasmids_presence["chr"].to_list():
            names.append(host)
            plasmids = plasmid_presence[plasmid_presence["chr"]==host]["plasmid"].to_list()

            for subcomm in big_subcomms:
                abs_presence = len(typing[typing["plasmid"].isin(plasmids) & (typing["type"]==subcomm)]["plasmid"].to_list())
                if abs_presence == 0:
                    subcomm_presence[subcomm].append(0)
                else:
                    subcomm_presence[subcomm].append(1)

        presence_df = pd.DataFrame(data=presence_absence, index=names)
        presence_df.sort_index(inplace=True)
        presence_df.sort_index(axis=1, inplace=True)
        presence_df.to_csv(output.plasmid_presence_absence, sep="\t")

rule phylofactor:
    input:
        tree = config["host_tree"],
        traits = config["output_dir"] + "/host_presence/presence_per_host.tsv"
    output:
        tree_vis = config["output"] + "phylofactor/{cluster}/tree.pdf"
        out_dir = config["output"] + "phylofactor/{cluster}"
    params:
        cluster = lambda wildcards: wildcards.cluster
    conda: "phylofactor"
    resources:
        mem_mb=lambda wildcards, attempt: 20000*attempt
    shell:
        "R < scripts/phylofactor.R {input.tree} {input.traits} {params.cluster} {output.out_dir} --no-save"

rule post_phylofactor:
    pass

rule dcj_trees:
    input:
        pling_out = config["output_dir"] + "/pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '')
    output:
        submatrices_dir = directory(config["output_dir"]"/pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '') + "/submatrices")
    conda:
        "pling"
    resources:
        pass
    shell:
        "pling submatrix {input.pling_out} --vis_trees"

rule parsnp:
    input:
        fastas = lambda wildcards: get_list(wildcards.cluster)
    output:
        parsnp_dir = directory(config["output_dir"] + "/{cluster}")
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
        dcj_dists = config["output_dir"] + "/pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '') + "/all_plasmids_distances.tsv"
    output:
        plot = "dcj_distr/hist_plot.png"
        stats = "dcj_distr/stats.txt"
    run:
        import pandas as pd
        import seaborn as sns
        import matplotlib.pyplot as plt

        dists = pd.read_csv(input.dcj_dists, sep="\t")
        fig, ax = plt.subplots()
        sns.histplot(data=dists,x="distance", ax=ax, discrete=True)
        plt.savefig(output.plot)

        with open(output.stats) as f:
            f.write("mode:"+str(dists["distance"].mode())+"\n")
            f.write("median:"+str(dists["distance"].median())+"\n")
            f.write("mean:"+str(dists["distance"].mean())+"\n")

rule cluster_specs:
    input:
        typing = config["output_dir"] + "pling_d" + config["dcj-indel"] + "_c" + config["containment"].replace(".", '') + "/dcj_thresh_" + config["dcj-indel"] + "_graph/objects/typing.tsv"
        mob = config["output_dir"] + "/mobtyper_results.txt"
    output:
        tsv = config["output_dir"] + "/cluster_specs.tsv"
    script:
        "scripts/cluster_specs.py"

rule rep_types:
    pass

rule boundary:
    pass
