# IBD-segment

A Nextflow pipeline for detecting identity-by-descent (IBD) segments with
[Hap-IBD](https://github.com/browning-lab/hap-ibd), summarizing genome-wide IBD
sharing, identifying nested IBD-sharing communities with Louvain and Leiden,
and consolidating genetically similar terminal communities by iterative
genotype-based \(Fst\) clumping.

The pipeline is configured for SLURM and was developed for the Digital Research
Alliance of Canada Narval cluster.

## Workflow

For each configured autosome, the pipeline:

1. filters variants by minor allele frequency and variant missingness;
2. verifies that retained genotypes are phased and non-missing;
3. verifies that sample IDs and their order match across chromosomes;
4. detects IBD segments with Hap-IBD using a chromosome-specific genetic map;
5. optionally removes complete segments overlapping excluded genomic regions;
6. calculates chromosome-specific per-pair IBD summary statistics;
7. combines chromosome summaries into an exact genome-wide pair summary;
8. performs fixed-depth recursive Louvain and/or Leiden clustering;
9. prepares an independent common-SNP, QC-filtered, LD-pruned genotype dataset;
10. iteratively merges eligible terminal clusters when Hudson
    \(F_{ST}<0.0005\), recalculating \(F_{ST}\) after every merge.

The genotype-preparation branch can run concurrently with Hap-IBD and the IBD
summary branch. Both clustering methods reuse the same pruned genotype dataset.

PhaseIBD and the former `MergeIBD` process are not included in version 2.

## Repository structure

```text
IBD-segment/
├── IBD_Pipeline.nf
├── nextflow.config
├── Run_IBD.sh
├── assets/
│   ├── chr_in_chrom_field/
│   │   ├── plink.chrchr1.GRCh38.map
│   │   └── ...
│   ├── no_chr_in_chrom_field/
│   │   ├── plink.chr1.GRCh38.map
│   │   └── ...
│   └── genome_gap_hg38_and_MHC.noheader.bed
└── bin/
    ├── hap-ibd.jar
    ├── IBD_per_pair.py
    ├── cluster_ibd_graph.py
    └── fst_clump.py
```

## Input requirements

Input data must consist of chromosome-split, bgzip-compressed VCF files and
their tabix indexes. Each VCF must:

- contain phased diploid genotypes using `|` as the allele separator;
- contain the same samples in the same order across chromosomes;
- contain no missing alleles after Hap-IBD QC;
- use chromosome identifiers compatible with the selected maps;
- be accompanied by a `.tbi` index.

Example:

```text
cohort.chr1.vcf.gz
cohort.chr1.vcf.gz.tbi
...
cohort.chr22.vcf.gz
cohort.chr22.vcf.gz.tbi
```

The paths use a `{chr}` placeholder:

```groovy
input_pattern = '/path/to/cohort.chr{chr}.vcf.gz'
input_index_pattern = '/path/to/cohort.chr{chr}.vcf.gz.tbi'
chromosomes = 1..22
```

## Software requirements

| Process | Required software |
|---|---|
| `QC_VCF` | bcftools and tabix |
| `VALIDATE_SAMPLES` | Standard Linux utilities |
| `HAP_IBD` | Java, bcftools, bedtools and Hap-IBD |
| `PER_PAIR_CHROMOSOME` | Python standard library |
| `PER_PAIR_GENOMEWIDE` | Python standard library |
| `CLUSTER_GRAPH` | Python with `igraph` |
| `PREPARE_FST_DATA` | bcftools and a recent PLINK 2 |
| `FST_CLUMP` | Python standard library and a recent PLINK 2 |

Load Nextflow before launching the workflow:

```bash
module load nextflow
```

### Python environment

The clustering script requires `python-igraph`. The summary and FST-clumping
scripts otherwise use only the Python standard library.

```bash
module load python/3.14.2

python3 -m venv ~/virtualenvs/ibd_pipeline
source ~/virtualenvs/ibd_pipeline/bin/activate

python3 -m pip install --upgrade pip
python3 -m pip install igraph

python3 -c "import igraph; print(igraph.__version__)"
```

Set the environment path in `nextflow.config`:

```groovy
python_venv = "${System.getenv('HOME')}/virtualenvs/ibd_pipeline"
```

### PLINK 2 module

`PREPARE_FST_DATA` and `FST_CLUMP` require a module that supplies the `plink2`
executable with PGEN and `--fst` support. Identify an available version on the
cluster and verify it before running:

```bash
module spider plink
module load AVAILABLE_PLINK2_MODULE
command -v plink2
plink2 --version
plink2 --help fst
```

The configuration value must contain only the module name, not the words
`module load`:

```groovy
plink2_module = 'AVAILABLE_PLINK2_MODULE'
```

The process itself runs `module load ${params.plink2_module}`. A value such as
`'module load plink/...'` would incorrectly generate `module load module load
plink/...`.

## Configuration

All user-controlled parameters are defined under `params` in
`nextflow.config`.

### Programs

```groovy
hapibd_jar = "${projectDir}/bin/hap-ibd.jar"
per_pair_script = "${projectDir}/bin/IBD_per_pair.py"
clustering_script = "${projectDir}/bin/cluster_ibd_graph.py"
fst_clump_script = "${projectDir}/bin/fst_clump.py"
```

### Genetic maps and chromosome names

```groovy
chr_in_chrom_field = true
```

Set this according to the VCF `CHROM` field:

- `true` for `chr1`, ..., `chr22`;
- `false` for `1`, ..., `22`.

The pipeline selects the corresponding bundled GRCh38 maps from `assets/`.

### Hap-IBD variant QC

```groovy
qc {
    min_maf = 0.001
    max_variant_missingness = 0.0
}
```

Hap-IBD requires phased genotypes without missing alleles. The pipeline checks
this explicitly after variant filtering.

### Hap-IBD parameters

```groovy
hapibd {
    min_seed_cm = 2.0
    min_extend_cm = 1.0
    min_output_cm = 3.0
    min_markers = 200
    min_mac = 2
    max_gap_bp = 1000
}
```

These correspond to Hap-IBD's `min-seed`, `min-extend`, `min-output`,
`min-markers`, `min-mac`, and `max-gap` arguments.

### Excluded genomic regions

```groovy
remove_gaps = true
gap_file = "${projectDir}/assets/genome_gap_hg38_and_MHC.noheader.bed"
```

When enabled, a complete IBD segment is removed if it overlaps at least one BED
interval. The BED file uses zero-based, half-open coordinates. Chromosome names
are normalized internally to match the VCF and Hap-IBD output.

### Per-pair summaries

```groovy
summary {
    segment_threshold_cm = 5.0
}
```

The summary reports statistics across all segments and separately counts and
sums segments at or above the configured threshold.

Chromosome summaries retain sufficient statistics, including summed squared
segment lengths. `IBD_per_pair.py --aggregate-summaries` combines them without
rereading the much larger raw Hap-IBD files and calculates genome-wide means
and standard deviations correctly.

Temporary SQLite databases are placed in node-local `$SLURM_TMPDIR`. This is
important because SQLite random I/O can be extremely slow on Lustre.

### Fixed-depth recursive clustering

```groovy
Louvain = true
Leiden = true

n_louvain = 3
n_leiden = 3

clustering {
    weight_column = 'total_IBD_length_x_mean_IBD_length_cM2'
    min_cluster_size = 20

    // Retained temporarily for CLI compatibility; ignored by the revised script.
    min_modularity_gain = 0.0
    auto_select = false

    leiden_resolution = 1.0
    seed = 2026
}
```

With a depth of three:

- `level_0` assigns the full cohort to `C1`;
- `level_1` clusters the complete weighted graph;
- `level_2` independently reclusters eligible level-1 communities;
- `level_3` independently reclusters eligible level-2 communities.

A community is eligible when its size is at least `min_cluster_size`.
Ineligible, edgeless, or unsplittable communities are carried forward unchanged
so every sample always has labels through `level_3`.

The script always materializes the requested depth. Modularity and modularity
gain remain diagnostics only: they do not stop refinement or choose a terminal
level. The deepest requested level is used for FST clumping.

Hierarchical labels are deterministic within a run structure. For example:

```text
ID        level_0  level_1  level_2  level_3
11100001  C1       C1.1     C1.1.1   C1.1.1.1
11100003  C1       C1.1     C1.1.2   C1.1.2.1
```

#### Edge weights

`weight_column` selects a numeric column from the genome-wide pair summary.
Larger values represent stronger IBD sharing.

| Value | Interpretation |
|---|---|
| `total_IBD_length_cM` | Total length of all shared segments |
| `total_IBD_length_x_mean_IBD_length_cM2` | Total length multiplied by mean segment length |
| `number_of_IBD_segments` | Number of shared segments |
| `mean_IBD_length_cM` | Mean segment length |
| `max_IBD_length_cM` | Longest shared segment |
| `number_of_segments_ge_5_cM` | Number of segments at least 5 cM long |
| `total_length_of_segments_ge_5_cM` | Total length of segments at least 5 cM long |
| `total_length_x_mean_length_of_segments_ge_5_cM_cM2` | Threshold-specific total length multiplied by mean length |

Threshold-specific column names change with `segment_threshold_cm`.
`sum_squared_IBD_length_cM`, `min_IBD_length_cM`, and `sd_IBD_length_cM` are
valid numeric columns but are not recommended as measures of total graph-edge
strength.

### FST genotype preparation and clumping

```groovy
fst {
    enabled = true

    min_maf = 0.05
    max_variant_missingness = 0.01
    hwe_pvalue = 1e-12

    prune_window_kb = 1000
    prune_step_variants = 1
    prune_r2 = 0.1

    cluster_column = 'level_3'
    clump_threshold = 0.0005
    min_cluster_size = 20
}
```

`PREPARE_FST_DATA` concatenates the chromosome VCFs in numeric order, retains
autosomal biallelic A/C/G/T SNPs, applies FST-specific MAF, missingness, and HWE
filters, performs LD pruning, and writes one reusable PGEN dataset.

The HWE threshold is deliberately lenient because a pooled, structured cohort
can deviate from HWE for biological reasons such as the Wahlund effect. Its
purpose here is to remove only extreme failures.

For each enabled graph method, `FST_CLUMP` then:

1. reads the configured terminal membership column;
2. computes all pairwise Hudson FST estimates with PLINK 2;
3. identifies the eligible pair with the lowest finite FST;
4. merges it if its FST is below `clump_threshold`;
5. recalculates all pairwise FST estimates using the updated memberships;
6. repeats until no eligible pair remains below the threshold.

Both clusters must contain at least `fst.min_cluster_size` samples to be
eligible. Smaller terminal clusters remain unchanged. Recalculation after each
merge avoids transitive connected-component chaining.

Final labels (`Cluster_1`, `Cluster_2`, ...) are assigned deterministically by
decreasing cluster size, with the retained internal label used to break ties.

### Resources

```groovy
resources {
    qc          { cpus = 2; memory = '8 GB';   time = '2h' }
    hapibd      { cpus = 8; memory = '175 GB'; time = '2h' }
    summary     { cpus = 1; memory = '16 GB';  time = '6h' }
    clustering  { cpus = 4; memory = '86 GB';  time = '4h' }
    fst_prepare { cpus = 8; memory = '64 GB';  time = '24h' }
    fst_clump   { cpus = 8; memory = '32 GB';  time = '48h' }
}
```

The summary implementation is single-threaded. Node-local SQLite storage is
more important than additional CPUs for this step. FST preparation may require
substantial node-local disk space for the temporary concatenated BCF.

## Running the pipeline

```bash
git clone https://github.com/JustinPelletier/IBD-segment.git
cd IBD-segment

module load nextflow

nextflow run IBD_Pipeline.nf \
    -c nextflow.config \
    -w ~/scratch/IBD-segment-work
```

Resume an interrupted run with:

```bash
nextflow run IBD_Pipeline.nf \
    -c nextflow.config \
    -resume \
    -w ~/scratch/IBD-segment-work
```

Keep the work directory under scratch and do not delete it until the workflow
and output validation are complete.

If an offline node cannot reach `www.nextflow.io`, Nextflow may print a version
check `curl` timeout before launching the installed version. This is separate
from pipeline process failures when Nextflow continues to launch normally.

## Outputs

```text
results/
├── qc/
│   ├── chr*.qc.stats.txt
│   ├── chr*.qc.summary.tsv
│   └── cohort.samples.txt
├── segments/
│   └── chr*.hapibd.ibd.gz
├── logs/
│   └── chr*.hapibd.log
├── per_pair/
│   ├── genomewide.per_pair.tsv.gz
│   └── by_chromosome/
│       └── chr*.per_pair.tsv.gz
├── fst/
│   └── genotypes/
│       ├── fst.pruned.pgen
│       ├── fst.pruned.pvar.zst
│       ├── fst.pruned.psam
│       ├── fst.prune.in
│       ├── fst.prepare.summary.tsv
│       └── fst.prepare.log
└── clustering/
    ├── louvain.membership.tsv.gz
    ├── louvain.selected_membership.tsv.gz
    ├── louvain.diagnostics.tsv
    ├── louvain.cluster_sizes.tsv
    ├── leiden.membership.tsv.gz
    ├── leiden.selected_membership.tsv.gz
    ├── leiden.diagnostics.tsv
    ├── leiden.cluster_sizes.tsv
    ├── louvain/
    │   ├── louvain.pairwise_fst.tsv.gz
    │   ├── louvain.fst_merge_history.tsv
    │   ├── louvain.final_membership.tsv.gz
    │   └── louvain.fst_clumping_summary.tsv
    └── leiden/
        ├── leiden.pairwise_fst.tsv.gz
        ├── leiden.fst_merge_history.tsv
        ├── leiden.final_membership.tsv.gz
        └── leiden.fst_clumping_summary.tsv
```

### Hap-IBD segments

The `chr*.hapibd.ibd.gz` files are headerless and retain Hap-IBD's eight
columns: two sample IDs and haplotype indices, chromosome, start position, end
position, and segment length in centimorgans.

### Per-pair summaries

Each pair summary contains:

- `ID1` and `ID2`;
- `number_of_IBD_segments`;
- `total_IBD_length_cM`;
- `sum_squared_IBD_length_cM`;
- `mean_IBD_length_cM`;
- `min_IBD_length_cM`;
- `max_IBD_length_cM`;
- `sd_IBD_length_cM`;
- threshold-specific segment count, total length, and derived weight columns.

Only pairs sharing at least one detected segment are written. Absent pairs
implicitly have zero IBD sharing. The complete sample list is supplied to the
clustering script, so samples without detected edges remain graph vertices.

### Clustering outputs

| File | Description |
|---|---|
| `{method}.membership.tsv.gz` | Membership at `level_0` through the requested terminal level |
| `{method}.selected_membership.tsv.gz` | Two-column membership for the deepest requested level |
| `{method}.diagnostics.tsv` | Cluster count, modularity, gain, eligible/split/carried communities, and terminal status by level |
| `{method}.cluster_sizes.tsv` | Cluster size and further-refinement eligibility at each level |

### FST clumping outputs

| File | Description |
|---|---|
| `{method}.pairwise_fst.tsv.gz` | Pairwise Hudson FST estimates at every recalculation iteration |
| `{method}.fst_merge_history.tsv` | Ordered merge history with pre-merge FST and cluster sizes |
| `{method}.final_membership.tsv.gz` | Original hierarchical membership plus `pre_fst_cluster` and `final_cluster` |
| `{method}.fst_clumping_summary.tsv` | Initial/final cluster counts, merges, threshold, and sample count |

## Validation

### Genome-wide IBD summary

Verify gzip integrity:

```bash
gzip -t results/per_pair/genomewide.per_pair.tsv.gz && \
echo 'PASS: genome-wide gzip file is valid'
```

The segment total must match the sum of chromosome summaries:

```bash
chromosome_segments=$(
    for file in results/per_pair/by_chromosome/chr*.per_pair.tsv.gz
    do
        zcat "$file" |
        awk 'NR > 1 {sum += $3} END {printf "%.0f\n", sum}'
    done |
    awk '{sum += $1} END {printf "%.0f\n", sum}'
)

genomewide_segments=$(
    zcat results/per_pair/genomewide.per_pair.tsv.gz |
    awk 'NR > 1 {sum += $3} END {printf "%.0f\n", sum}'
)

printf 'Chromosome total: %s\n' "$chromosome_segments"
printf 'Genome-wide total: %s\n' "$genomewide_segments"

[[ "$chromosome_segments" == "$genomewide_segments" ]] && \
echo 'PASS: all IBD segments were aggregated'
```

Total IBD length should also agree, allowing a negligible floating-point
difference caused by decimal formatting when summaries are written and reread.

### Clustering and FST results

Confirm that both membership files contain `level_3`:

```bash
for file in results/clustering/{louvain,leiden}.membership.tsv.gz
do
    zcat "$file" | head -n 1
done
```

Confirm that every genotype sample received a final label:

```bash
for method in louvain leiden
do
    printf '%s\t' "$method"
    zcat "results/clustering/${method}/${method}.final_membership.tsv.gz" |
    awk -F '\t' 'END {print NR - 1}'
done
```

Inspect the clumping summaries and merge histories:

```bash
column -t -s $'\t' results/clustering/*/*.fst_clumping_summary.tsv
column -t -s $'\t' results/clustering/louvain/louvain.fst_merge_history.tsv | head
```

## Standalone usage

### Chromosome and genome-wide summaries

```bash
python3 bin/IBD_per_pair.py \
    --input chr1.hapibd.ibd.gz \
    --output chr1.per_pair.tsv.gz \
    --threshold-cm 5 \
    --temporary-directory "${SLURM_TMPDIR:-.}"

python3 bin/IBD_per_pair.py \
    --input chr*.per_pair.tsv.gz \
    --output genomewide.per_pair.tsv.gz \
    --threshold-cm 5 \
    --aggregate-summaries \
    --temporary-directory "${SLURM_TMPDIR:-.}"
```

### Graph clustering

The deprecated modularity arguments are accepted for compatibility but ignored.

```bash
python3 bin/cluster_ibd_graph.py \
    --edges genomewide.per_pair.tsv.gz \
    --samples cohort.samples.txt \
    --method louvain \
    --max-levels 3 \
    --weight-column total_IBD_length_x_mean_IBD_length_cM2 \
    --min-cluster-size 20 \
    --min-modularity-gain 0 \
    --resolution 1.0 \
    --seed 2026 \
    --auto-select false \
    --output-prefix louvain
```

### Iterative FST clumping

```bash
python3 bin/fst_clump.py \
    --pgen fst.pruned.pgen \
    --pvar fst.pruned.pvar.zst \
    --psam fst.pruned.psam \
    --membership louvain.membership.tsv.gz \
    --cluster-column level_3 \
    --threshold 0.0005 \
    --min-cluster-size 20 \
    --plink2 plink2 \
    --threads 8 \
    --memory-mb 32000 \
    --output-prefix louvain
```

## Troubleshooting

### `No such variable: method`

Dynamic method-specific publication directories must use a closure:

```groovy
publishDir {
    "${params.outdir}/clustering/${method}"
},
    mode: 'copy',
    overwrite: true
```

### Lmod reports an unknown module named `module`

Set `plink2_module` to the module name only:

```groovy
plink2_module = 'AVAILABLE_PLINK2_MODULE'
```

Do not include `module load` in the parameter value.

### Summary task is unexpectedly slow

Verify that `.command.sh` contains:

```bash
--temporary-directory "${SLURM_TMPDIR}"
```

The temporary `.sqlite.tmp` file should be on node-local storage rather than in
the Nextflow work directory on Lustre.

## References

If Hap-IBD is used in a publication, cite:

> Zhou Y, Browning SR, Browning BL. A fast and simple method for detecting
> identity-by-descent segments in large-scale data. *American Journal of Human
> Genetics*. 2020;106(4):426–437.
> https://doi.org/10.1016/j.ajhg.2020.02.010

Relevant clustering and FST references:

> Blondel VD, Guillaume J-L, Lambiotte R, Lefebvre E. Fast unfolding of
> communities in large networks. *Journal of Statistical Mechanics: Theory and
> Experiment*. 2008;2008(10):P10008.

> Traag VA, Waltman L, van Eck NJ. From Louvain to Leiden: guaranteeing
> well-connected communities. *Scientific Reports*. 2019;9:5233.

> Bhatia G, Patterson N, Sankararaman S, Price AL. Estimating and interpreting
> FST: The impact of rare variants. *Genome Research*. 2013;23(9):1514–1521.

## Author

- Justin Pelletier
- McGill University
- justin.pelletier2@mcgill.ca
