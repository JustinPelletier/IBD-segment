# IBD-segment

A Nextflow pipeline for detecting identity-by-descent (IBD) segments with
[Hap-IBD](https://github.com/browning-lab/hap-ibd) from chromosome-split,
phased VCF files. The pipeline performs variant-level quality control, prepares
genetic maps, detects and summarizes IBD segments, and optionally identifies
IBD-sharing communities with Louvain and Leiden clustering.

The pipeline is configured for SLURM execution and was developed for the
Digital Research Alliance of Canada Narval cluster.

## Workflow

For each configured chromosome, the pipeline:

1. filters variants according to minor allele frequency and variant missingness;
2. confirms that all retained genotypes are phased and non-missing;
3. confirms that sample IDs and their order are identical across chromosomes;
4. optionally generates a marker-matched genetic map with PLINK;
5. detects IBD segments with Hap-IBD;
6. optionally removes complete segments overlapping genome-gap intervals;
7. generates chromosome-specific and genome-wide per-pair summaries;
8. optionally performs recursive Louvain and Leiden clustering.

PhaseIBD and the former `MergeIBD` process are not included in version 2.

## Repository structure

```text
IBD-segment/
├── IBD_Pipeline.nf
├── nextflow.config
├── Run_IBD.sh
├── assets/
│   └── empty.gap.bed
└── bin/
    ├── hap-ibd.jar
    ├── IBD_per_pair.py
    ├── cluster_ibd_graph.py
    └── HapMap/
        ├── GRCh37/
        └── GRCh38/
```

## Input requirements

Input data must be split by chromosome and supplied as bgzip-compressed VCF
files. Each input VCF must:

- contain phased diploid genotypes using `|` as the allele separator;
- contain the same samples in the same order across chromosomes;
- contain no missing alleles after quality control;
- use chromosome names matching the selected genetic maps, such as `1` or
  `chr1`.

For example:

```text
cohort.chr1.vcf.gz
cohort.chr2.vcf.gz
...
cohort.chr22.vcf.gz
```

The input pattern is specified with a `{chr}` placeholder:

```groovy
input_pattern = '/path/to/data/cohort.chr{chr}.vcf.gz'
```

The placeholder can occur anywhere in the path or filename.

## Software requirements on Narval

The Nextflow processes load the required Narval modules or activate the Python
environment through `beforeScript` directives.

| Process | Required software |
|---|---|
| `QC_VCF` | `module load bcftools` |
| `VALIDATE_SAMPLES` | Standard Linux utilities |
| `GENETIC_MAP` | `module load StdEnv/2020` and `module load plink/1.9b_6.21-x86_64` |
| `HAP_IBD` | `module load java/25.36`, `module load bcftools`, and `module load bedtools` |
| Summary processes | Python virtual environment |
| Clustering | Python virtual environment with `igraph` |

Nextflow itself must be loaded before launching the workflow:

```bash
module load nextflow
```

### Python virtual environment

Create the environment once before running the pipeline. The clustering script
requires `python-igraph`; the per-pair summarization script otherwise uses only
the Python standard library.

```bash
module load python

virtualenv --no-download ~/virtualenvs/ibd_pipeline
source ~/virtualenvs/ibd_pipeline/bin/activate

pip install --no-index igraph
```

Confirm that igraph is available:

```bash
python3 -c "import igraph; print(igraph.__version__)"
```

The corresponding path in `nextflow.config` is:

```groovy
python_venv = "${System.getenv('HOME')}/virtualenvs/ibd_pipeline"
```

## Configuration

All user-controlled parameters are defined in `nextflow.config`.

### Input and output

```groovy
input_pattern = '/path/to/data/cohort.chr{chr}.vcf.gz'
chromosomes = 1..22
outdir = 'results'
```

### Variant-level quality control

```groovy
qc {
    min_maf = 0.01
    max_variant_missingness = 0.0
}
```

Hap-IBD requires phased genotypes without missing alleles. A missingness value
greater than zero is therefore valid only when every variant retained after QC
still has complete genotypes. The pipeline explicitly checks this requirement
before running Hap-IBD.

### Genetic maps

```groovy
assembly = 'GRCh38'
genetic_map_pattern = "${projectDir}/bin/HapMap/GRCh38/no_chr.chr{chr}.GRCh38.gmap"
generate_custom_maps = true
```

Both `chr1`-style and `1`-style HapMap maps for GRCh37 and GRCh38 are included.
The chromosome naming convention must match the VCF.

When `generate_custom_maps = true`, PLINK interpolates the genetic position of
each variant retained after QC and produces a map containing exactly those
markers. When it is `false`, the selected HapMap file is passed directly to
Hap-IBD.

### Hap-IBD parameters

```groovy
hapibd {
    min_seed_cm = 2.0
    min_extend_cm = 1.0
    min_output_cm = 2.0
    min_markers = 200
    min_mac = 2
    max_gap_bp = 1000
}
```

These correspond to the Hap-IBD arguments `min-seed`, `min-extend`,
`min-output`, `min-markers`, `min-mac`, and `max-gap`.

### Genome-gap filtering

```groovy
remove_gaps = true
gap_pattern = '/path/to/gaps/GRCh38.chr{chr}.gap.bed'
```

When enabled, any complete IBD segment overlapping at least one interval in the
chromosome-specific BED file is excluded. BED files must use zero-based,
half-open coordinates and chromosome names matching the Hap-IBD output.

### Per-pair summary threshold

```groovy
summary {
    segment_threshold_cm = 5.0
}
```

In addition to statistics across all segments, the per-pair summary reports the
number and total length of segments greater than or equal to this threshold.

### Clustering

```groovy
Louvain = true
Leiden = true

n_louvain = 3
n_leiden = 3
```

`Louvain` and `Leiden` determine whether each clustering process is launched.
The `n_louvain` and `n_leiden` values specify maximum recursive refinement
depths, not repeated stochastic runs:

- level 0 places the complete cohort in one initial community;
- level 1 clusters the complete IBD-sharing graph;
- level 2 reclusters each eligible level-1 community independently;
- subsequent levels repeat this refinement.

Fine-tuning parameters are defined as:

```groovy
clustering {
    weight_column = 'total_IBD_length_cM'
    min_cluster_size = 20
    min_modularity_gain = 0.0001
    auto_select = true
    leiden_resolution = 1.0
    seed = 2026
}
```

`min_cluster_size` specifies the minimum size required for a community to be
submitted to another refinement level. Smaller communities remain in the
output but are not refined further.

When `auto_select = true`, recursive refinement stops after a level whose
global weighted modularity gain is smaller than `min_modularity_gain`. Among
the computed levels, the level with the highest modularity is selected. When
`auto_select = false`, the deepest successfully computed level is selected.

#### Edge weights

The `weight_column` parameter selects a column from the genome-wide per-pair
summary. Larger values indicate stronger IBD sharing and increase the tendency
of two samples to be placed in the same community.

Available and potentially useful values include:

| Value | Interpretation |
|---|---|
| `total_IBD_length_cM` | Total length of all shared segments; recommended default |
| `number_of_IBD_segments` | Number of shared segments, regardless of length |
| `mean_IBD_length_cM` | Mean shared-segment length |
| `max_IBD_length_cM` | Longest segment shared by the pair |
| `number_of_segments_ge_5_cM` | Count of segments at or above the configured threshold |
| `total_length_of_segments_ge_5_cM` | Total length of segments at or above the threshold |

The names of the last two columns change when `segment_threshold_cm` changes.
For example, a threshold of `3.0` produces columns ending in `_ge_3_cM`.

`min_IBD_length_cM` and `sd_IBD_length_cM` are technically valid numeric
columns but are generally not meaningful measures of total pairwise IBD
sharing and are not recommended as graph weights.

### Computational resources

Resources can be adjusted independently for each process:

```groovy
resources {
    qc         { cpus = 2; memory = '8 GB';   time = '2h' }
    map        { cpus = 1; memory = '8 GB';   time = '2h' }
    hapibd     { cpus = 8; memory = '175 GB'; time = '4h' }
    summary    { cpus = 1; memory = '8 GB';   time = '12h' }
    clustering { cpus = 4; memory = '32 GB';  time = '12h' }
}
```

## Running the pipeline

Clone the repository and enter it:

```bash
git clone https://github.com/JustinPelletier/IBD-segment.git
cd IBD-segment
```

Edit `nextflow.config`, then launch the workflow:

```bash
module load nextflow

nextflow run IBD_Pipeline.nf \
    -c nextflow.config \
    -w ~/scratch/IBD-segment-work
```

Resume an interrupted or partially completed run with:

```bash
nextflow run IBD_Pipeline.nf \
    -c nextflow.config \
    -resume \
    -w ~/scratch/IBD-segment-work
```

The Nextflow work directory should be placed under `~/scratch` on Narval to
support high-throughput temporary I/O. Do not delete it until the run has
completed successfully and all required outputs have been published.

## Outputs

```text
results/
├── qc/
│   ├── chr*.qc.stats.txt
│   └── cohort.samples.txt
├── maps/
│   └── chr*.custom.map
├── segments/
│   └── chr*.hapibd.ibd.gz
├── logs/
│   └── chr*.hapibd.log
├── per_pair/
│   ├── genomewide.per_pair.tsv.gz
│   └── by_chromosome/
│       └── chr*.per_pair.tsv.gz
└── clustering/
    ├── louvain/
    └── leiden/
```

### Hap-IBD segment files

The `chr*.hapibd.ibd.gz` files are headerless and preserve the standard
eight-column Hap-IBD format:

1. first sample ID;
2. first haplotype index;
3. second sample ID;
4. second haplotype index;
5. chromosome;
6. segment start position in base pairs;
7. segment end position in base pairs;
8. segment length in centimorgans.

### Per-pair summaries

Each per-pair summary contains:

- `ID1` and `ID2`;
- `number_of_IBD_segments`;
- `total_IBD_length_cM`;
- `mean_IBD_length_cM`;
- `min_IBD_length_cM`;
- `max_IBD_length_cM`;
- `sd_IBD_length_cM`;
- number of segments at or above the configured threshold;
- total length of segments at or above the threshold.

Only pairs sharing at least one detected segment are written. Pairs absent from
the summary implicitly have zero IBD sharing. The complete `cohort.samples.txt`
file is supplied to the clustering script so participants without detected
edges remain present as isolated graph vertices.

### Clustering outputs

Each enabled algorithm produces:

| File | Description |
|---|---|
| `{method}.membership.tsv.gz` | Community membership at every computed level |
| `{method}.selected_membership.tsv.gz` | Membership at the selected level |
| `{method}.diagnostics.tsv` | Modularity, gain, cluster count, splits and stopping reason by level |
| `{method}.cluster_sizes.tsv` | Community sizes at every level |

Automatic modularity selection is a practical model-selection heuristic, not a
formal estimate of the true number of populations. Final clusters should also
be evaluated for stability, size, genetic composition, geography, ancestry and
biological interpretability.

## Standalone script usage

Create a summary for one chromosome:

```bash
source ~/virtualenvs/ibd_pipeline/bin/activate

python3 bin/IBD_per_pair.py \
    --input chr1.hapibd.ibd.gz \
    --output chr1.per_pair.tsv.gz \
    --threshold-cm 5
```

Run Louvain clustering directly:

```bash
python3 bin/cluster_ibd_graph.py \
    --edges genomewide.per_pair.tsv.gz \
    --samples cohort.samples.txt \
    --method louvain \
    --max-levels 3 \
    --weight-column total_IBD_length_cM \
    --min-cluster-size 20 \
    --min-modularity-gain 0.0001 \
    --resolution 1.0 \
    --seed 2026 \
    --auto-select true \
    --output-prefix louvain
```

## References

If Hap-IBD is used in a publication, cite:

> Zhou Y, Browning SR, Browning BL. A fast and simple method for detecting
> identity-by-descent segments in large-scale data. *American Journal of Human
> Genetics*. 2020;106(4):426–437. https://doi.org/10.1016/j.ajhg.2020.02.010

Relevant clustering references:

> Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E. Fast unfolding of
> communities in large networks. *Journal of Statistical Mechanics: Theory and
> Experiment*. 2008;2008(10):P10008.

> Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing
> well-connected communities. *Scientific Reports*. 2019;9:5233.

## Author

- Justin Pelletier
- McGill University
- justin.pelletier2@mcgill.ca
