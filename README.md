I keep rerunning the same stuff on data, so why not bodge it all together ¯\\\_(ツ)\_/¯

## What it does

Starting from a set of plasmid assemblies, the pipeline:

- clusters plasmids into communities and subcommunities with [pling](https://github.com/iqbal-lab-org/pling), and summarises DCJ-Indel distances, NJ trees and boundary vs internal distances per subcommunity
- types plasmids with [MOB-suite](https://github.com/phac-nml/mob-suite) and summarises each subcommunity
- for every subcommunity larger than `big_subcomm_size`:
  - builds a pangenome with [ggCallaroo](https://github.com/samhorsfield96/ggCallaroo) (ggCaller + Panaroo + Bakta) and computes relative core genome sizes
  - builds a pangenome graph with [PanGraph](https://github.com/neherlab/pangraph)
  - builds a core genome alignment with [Parsnp](https://github.com/marbl/parsnp)
  - finds host clades enriched for the subcommunity with [phylofactor](https://github.com/reptalex/phylofactor), and filters them

## Installation

### 1. Clone the repository

```bash
git clone https://github.com/babayagaofficial/yagas-plasmid-pipeline.git
cd yagas-plasmid-pipeline
```

### 2. Create the Snakemake environment

```bash
mamba env create -f envs/snakemake.yaml
conda activate yagas-plasmid-pipeline
```

This installs Snakemake 9, the SLURM executor plugin, and the Python packages used by the rules which run inside Snakemake's own environment. All other tools get their own conda environments (in `envs/`), which Snakemake creates automatically on the first run.

### 3. Install PanGraph

PanGraph is not available on conda, so **you need to download the binary yourself and add it to your `PATH`**. Binaries for each platform are on the [PanGraph releases page](https://github.com/neherlab/pangraph/releases); the pipeline was set up with version 1.4.0.

### 4. Download a Bakta database

ggCallaroo annotates the pangenome with Bakta, which needs a local database. See the [Bakta documentation](https://github.com/oschwengers/bakta#database) for how to download one, and set `bakta_db` in `config.yaml` to its path.

## Input

Set the following in `config.yaml`:

| Key | Description |
|---|---|
| `input_list` | Text file with one plasmid FASTA path per line. The file name without its extension is used as the plasmid ID. |
| `multifasta` | Alternatively, a single FASTA file containing all plasmids; record IDs are used as plasmid IDs. Set exactly one of `input_list` and `multifasta`, and leave the other empty. |
| `output_dir` | Output directory; use an absolute path here. |
| `chr_to_plasmid` | Tab-separated file with columns `chr` and `plasmid`, relating each plasmid ID to the ID of its host chromosome. |
| `host_tree` | Newick tree of the host chromosomes; tip labels must match the `chr` column of `chr_to_plasmid`. |
| `dcj-indel`, `containment` | pling's DCJ-Indel and containment distance thresholds. |
| `big_subcomm_size` | Subcommunities with more plasmids than this are analysed individually. |
| `ggcaller_cli_args`, `panaroo_cli_args` | Extra command-line arguments passed to ggCaller and Panaroo. |
| `bakta_db` | Path to the Bakta database. |
| `resources` | Memory (`mem_mb`), time (`runtime`, in minutes) and `threads` per rule; see [Resources](#resources). |

## Usage

### Locally

```bash
snakemake --cores 16 --sdm conda
```

### On a SLURM cluster

```bash
snakemake --executor slurm --jobs 100 --retries 2 --sdm conda --latency-wait 60 \
  --default-resources slurm_account=<account> slurm_partition=<partition>
```

Each rule is submitted as its own job, with the memory, time and threads set in `config.yaml`. A few lightweight rules run on the node you launch Snakemake from instead of being submitted.

Building the conda environments and loading the ggCallaroo workflow requires internet access (ggCallaroo is imported from GitHub, and phylofactor is installed from GitHub when its environment is built). If your compute nodes are offline, build the environments on a login node first:

```bash
snakemake --sdm conda --conda-create-envs-only
```

## Resources

Resources are set per rule in the `resources` section of `config.yaml`, and rules which are not listed use `default`:

```yaml
resources:
  default:          {mem_mb: 4000,  runtime: 60,  threads: 1}
  pling:            {mem_mb: 32000, runtime: 720, threads: 16}
```

When a job fails and is retried (`--retries`), its `mem_mb` and `runtime` are multiplied by the attempt number, so a job which runs out of memory or time is resubmitted with twice, then three times, as much. The values in the repository are placeholders; adjust them to your data.

## Output

All results are written to `output_dir`. Subcommunities analysed individually are referred to as `{cluster}`.

| Path | Contents |
|---|---|
| `pling_d{dcj}_c{containment}/` | pling output: communities, subcommunities (`dcj_thresh_{dcj}_graph/objects/typing.tsv`), distances, and DCJ-Indel trees in `submatrices/` |
| `mobtyper_results.txt` | MOB-typer results for all plasmids |
| `cluster_specs.tsv` | Size, average length, replicon, relaxase, MPF types and predicted mobility per subcommunity |
| `dcj_distr/` | Histogram and summary statistics of all DCJ-Indel distances |
| `boundary/` | Median and mean DCJ-Indel distances on subcommunity boundaries vs within subcommunities (`dcj_averages.tsv` and plots) |
| `cluster_lists/` | FASTA lists of the subcommunities analysed individually |
| `ggcallaroo/{cluster}/` | ggCaller, Panaroo and Bakta results; annotated pangenome files are in `annotated/` |
| `rel_core/` | Relative core genome sizes per subcommunity |
| `pangraph/{cluster}/` | PanGraph pangenome graph (`pangraph.json`, `pangraph.gfa`) |
| `parsnp/{cluster}/` | Parsnp core genome alignment |
| `host_presence/presence_per_host.tsv` | Presence/absence of each subcommunity in each host |
| `phylofactor/{cluster}/` | Host clades found by phylofactor, per-host rates and a tree plot |
| `post_phylofactor/{cluster}/` | Filtered clades, pruned host trees and plasmid FASTA lists for each clade |
| `logs/` | Log files for each rule; ggCallaroo's logs are in `ggcallaroo/{cluster}/logs/` |
