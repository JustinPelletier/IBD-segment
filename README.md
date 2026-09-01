# IBD-segment

A Nextflow pipeline for detecting identity-by-descent (IBD) segments with
[Hap-IBD](https://github.com/browning-lab/hap-ibd), summarizing genome-wide IBD
sharing, identifying nested IBD-sharing communities with Louvain and Leiden,
and consolidating genetically similar terminal communities by iterative
genotype-based Hudson \(F_{ST}\) clumping. The pipeline also summarizes
within- and between-cluster IBD sharing and produces annotated heatmaps of
final-cluster IBD sharing and Hudson \(F_{ST}\).

The pipeline is configured for SLURM and was developed for the Digital Research
Alliance of Canada Narval cluster.

## Workflow

The pipeline:

1. filters variants by minor allele frequency and variant missingness;
2. verifies that retained genotypes are phased and non-missing;
3. verifies that sample IDs and their order match across chromosomes;
4. optionally identifies close relatives using KING and a preliminary
   full-cohort Hap-IBD pass;
5. selects and applies one consistent unrelated sample set to every chromosome;
6. detects final IBD segments with Hap-IBD using chromosome-specific maps;
7. optionally removes complete segments overlapping excluded genomic regions;
8. calculates chromosome-specific and exact genome-wide pair summaries;
9. performs fixed-depth recursive Louvain and/or Leiden clustering;
10. prepares an independent common-SNP, QC-filtered, LD-pruned genotype dataset;
11. iteratively merges eligible terminal clusters below the configured Hudson
    \(F_{ST}\) threshold, recalculating \(F_{ST}\) after every merge;
12. summarizes within- and between-cluster IBD sharing after merging; and
13. recomputes final pairwise Hudson \(F_{ST}\) and generates annotated heatmaps.

The genotype-preparation branch can run concurrently with Hap-IBD and the IBD
summary branch. Both clustering methods reuse the same pruned genotype dataset
and final analysis cohort.

PhaseIBD and the former `MergeIBD` process are not included in version 2.

## Repository structure

```text
IBD-segment/
├── IBD_Pipeline.nf
├── nextflow.example.config
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
    ├── fst_clump.py
    ├── relatedness_selector.py
    ├── final_ibd_sharing.py
    ├── final_fst_heatmap.py
    ├── run_ibdne_clusters.py
    └── ibdne.23Apr20.ae9.jar
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
| `HAP_IBD_RELATEDNESS` | Java, bcftools, bedtools and Hap-IBD |
| `KING_RELATEDNESS` | bcftools and a recent PLINK 2 |
| `SELECT_UNRELATED_SAMPLES` | Python standard library |
| `APPLY_UNRELATED_SAMPLES` | bcftools and tabix |
| `HAP_IBD` | Java, bcftools, bedtools and Hap-IBD |
| `PER_PAIR_CHROMOSOME` | Python standard library |
| `PER_PAIR_GENOMEWIDE` | Python standard library |
| `CLUSTER_GRAPH` | Python with `igraph` |
| `PREPARE_FST_DATA` | bcftools and a recent PLINK 2 |
| `FST_CLUMP` | Python standard library and a recent PLINK 2 |
| `FINAL_IBD_SHARING` | Python with NumPy and Matplotlib |
| `FINAL_FST_HEATMAP` | Python with NumPy and Matplotlib, and a recent PLINK 2 |
| `RUN_IBDNE_FINAL_CLUSTERS` | Java, Python standard library and IBDNe |

Load Nextflow before launching the workflow:

```bash
module load nextflow
```

### Python environment

The clustering script requires `python-igraph`. Final plots require NumPy and
Matplotlib. The summary, relatedness-selection, and FST-clumping scripts
otherwise use only the Python standard library.

```bash
module load python/3.14.2

python3 -m venv ~/virtualenvs/ibd_pipeline
source ~/virtualenvs/ibd_pipeline/bin/activate

python3 -m pip install --upgrade pip
python3 -m pip install igraph numpy matplotlib

python3 -c "import igraph; print(igraph.__version__)"
python3 -c "import numpy, matplotlib; print(numpy.__version__, matplotlib.__version__)"
```

Set the environment path in your run configuration:

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

All user-controlled parameters are defined under `params`. Copy the bundled
example before adapting paths and resources:

```bash
cp nextflow.example.config nextflow.config
```

### Programs

```groovy
hapibd_jar = "${projectDir}/bin/hap-ibd.jar"
per_pair_script = "${projectDir}/bin/IBD_per_pair.py"
clustering_script = "${projectDir}/bin/cluster_ibd_graph.py"
fst_clump_script = "${projectDir}/bin/fst_clump.py"
relatedness_selector_script = "${projectDir}/bin/relatedness_selector.py"
final_ibd_script = "${projectDir}/bin/final_ibd_sharing.py"
final_fst_heatmap_script = "${projectDir}/bin/final_fst_heatmap.py"
ibdne_runner_script = "${projectDir}/bin/run_ibdne_clusters.py"
ibdne_jar = "${projectDir}/bin/ibdne.23Apr20.ae9.jar"
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

### Relatedness filtering

```groovy
relatedness {
    mode = 'discover_apply'
    sample_list = null
    removal_mode = 'optimal'
    king_cutoff = 0.0884
    ibd_total_cm_cutoff = 1000.0

    min_maf = 0.05
    max_variant_missingness = 0.01
    prune_window_variants = 50
    prune_step_variants = 10
    prune_r2 = 0.1
    indep_order = 1
}
```

Three modes are available:

- `discover_apply` runs KING and a preliminary full-cohort Hap-IBD pass,
  selects an unrelated subset, and applies it to all final analyses;
- `reuse` skips discovery and uses `relatedness.sample_list`;
- `disabled` retains the complete cohort.

`removal_mode = 'optimal'` retains a large mutually unrelated subset.
`remove_all` removes every participant belonging to at least one flagged pair.
The Hap-IBD cutoff is based on total genome-wide sharing from the preliminary
pass and is independent of the final clustering threshold.

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
    segment_threshold_cm = 3.0
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
    weight_column = 'total_length_of_segments_ge_3_cM'
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

A parent community is eligible for reclustering when its size is at least
`min_cluster_size`. After a proposed split, child communities smaller than the
minimum are pooled into one residual child. The split is accepted only when the
pooled residual and every retained child contain at least the minimum number of
samples. If the pooled residual remains undersized but at least two valid large
children exist, it is absorbed into the large child with which it shares the
greatest total IBD edge weight. Ties are resolved deterministically by child
size and sample ID. The proposed split is rejected only when fewer than two
valid children can remain, in which case the parent label is carried forward.

For example, proposed children of sizes 45, 30, 10, and 15 become children of
sizes 45, 30, and 25. The two undersized children are combined into the final
residual child. Edgeless and otherwise unsplittable parents are also carried
forward, so every sample always has labels through `level_3` and no accepted
split creates a cluster smaller than `min_cluster_size`.

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
| `number_of_segments_ge_{threshold}_cM` | Number of segments at or above the configured threshold |
| `total_length_of_segments_ge_{threshold}_cM` | Total length of segments at or above the configured threshold |
| `total_length_x_mean_length_of_segments_ge_{threshold}_cM_cM2` | Threshold-specific total length multiplied by mean length |

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
    clump_threshold = 0.001
    min_cluster_size = 15
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

### Final-cluster IBDNe

```groovy
ibdne {
    enabled = true
    min_cluster_size = 20

    mincm = 2.0
    nits = 1000
    nboots = 80
    filtersamples = true
    seed = 2026

    parallel_clusters = 2
    java_heap_gb = 24
}
```

`RUN_IBDNE_FINAL_CLUSTERS` runs after FST clumping. For each Louvain and/or
Leiden result, it streams the chromosome Hap-IBD files once and retains only
segments whose two samples belong to the same `final_cluster`. Final clusters
with fewer than `ibdne.min_cluster_size` participants are skipped and recorded
in the status table.

The number of unordered haplotype pairs is supplied explicitly as

\[
\mathrm{npairs}=\frac{2N(2N-2)}{2},
\]

where \(N\) is the complete final-cluster sample count. This prevents cluster
members with no retained within-cluster segment from being omitted from the
IBDNe sampling denominator. Cluster-specific segment files are created in
node-local `$SLURM_TMPDIR` and are not published.

The process runs up to `parallel_clusters` IBDNe analyses concurrently within
each clustering-method task. `java_heap_gb × parallel_clusters`, Python/Java
overhead, and filesystem buffers must fit within `resources.ibdne.memory`.

Download the official IBDNe JAR and place it at the configured path before
running. The current implementation targets `ibdne.23Apr20.ae9.jar`.

### Resources

```groovy
resources {
    qc                  { cpus = 2; memory = '4 GB';   time = '1h' }
    relatedness_hapibd  { cpus = 8; memory = '175 GB'; time = '2h' }
    relatedness_king    { cpus = 8; memory = '64 GB';  time = '12h' }
    relatedness_select  { cpus = 1; memory = '16 GB';  time = '2h' }
    hapibd              { cpus = 8; memory = '175 GB'; time = '2h' }
    summary             { cpus = 1; memory = '12 GB';  time = '6h' }
    clustering          { cpus = 4; memory = '64 GB';  time = '4h' }
    fst_prepare         { cpus = 8; memory = '8 GB';   time = '2h' }
    fst_clump           { cpus = 8; memory = '4 GB';   time = '2h' }
    final_ibd           { cpus = 1; memory = '16 GB';  time = '6h' }
    final_fst           { cpus = 8; memory = '32 GB';  time = '12h' }
    ibdne               { cpus = 16; memory = '64 GB'; time = '48h' }
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

cp nextflow.example.config nextflow.config

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
│   ├── cohort.samples.txt
│   └── unrelated/
│       └── chr*.unrelated.summary.tsv
├── relatedness/
│   ├── unrelated.samples.txt
│   ├── related_samples_removed.txt
│   ├── related_pairs.tsv.gz
│   ├── related_components.tsv.gz
│   ├── relatedness_selection_summary.tsv
│   ├── king.relatedness.summary.tsv
│   ├── king.kin0
│   ├── king.prune.in
│   ├── king.relatedness.log
│   ├── relatedness.genomewide.per_pair.tsv.gz
│   ├── discovery_segments/
│   │   └── chr*.relatedness.hapibd.ibd.gz
│   └── logs/
│       └── chr*.relatedness.hapibd.log
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
    │   ├── louvain.fst_clumping_summary.tsv
    │   ├── ibd_sharing/
    │   │   ├── louvain.ibd_sharing.tsv
    │   │   ├── louvain.within_cluster_distribution.tsv
    │   │   ├── louvain.ibd_mean_matrix.tsv
    │   │   ├── louvain.ibd_heatmap.png
    │   │   └── louvain.within_ibd_boxplot.png
    │   ├── final_fst/
    │   │   ├── louvain.final_fst.tsv
    │   │   ├── louvain.final_fst_matrix.tsv
    │   │   └── louvain.final_fst_heatmap.png
    │   └── ibdne/
    │       ├── louvain.ibdne_status.tsv
    │       └── Cluster_*/
    │           ├── Cluster_*.ne
    │           ├── Cluster_*.boot
    │           ├── Cluster_*.log
    │           ├── Cluster_*.pair.excl
    │           ├── Cluster_*.region.excl
    │           ├── Cluster_*.command.txt
    │           └── Cluster_*.driver.log
    └── leiden/
        ├── leiden.pairwise_fst.tsv.gz
        ├── leiden.fst_merge_history.tsv
        ├── leiden.final_membership.tsv.gz
        ├── leiden.fst_clumping_summary.tsv
        ├── ibd_sharing/
        │   ├── leiden.ibd_sharing.tsv
        │   ├── leiden.within_cluster_distribution.tsv
        │   ├── leiden.ibd_mean_matrix.tsv
        │   ├── leiden.ibd_heatmap.png
        │   └── leiden.within_ibd_boxplot.png
        ├── final_fst/
        │   ├── leiden.final_fst.tsv
        │   ├── leiden.final_fst_matrix.tsv
        │   └── leiden.final_fst_heatmap.png
        └── ibdne/
            ├── leiden.ibdne_status.tsv
            └── Cluster_*/
                ├── Cluster_*.ne
                ├── Cluster_*.boot
                ├── Cluster_*.log
                ├── Cluster_*.pair.excl
                ├── Cluster_*.region.excl
                ├── Cluster_*.command.txt
                └── Cluster_*.driver.log
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
| `{method}.diagnostics.tsv` | Cluster count, modularity, gain, eligible/split/carried communities, pooled small children, absorbed undersized residuals, rejected splits, and terminal status by level |
| `{method}.cluster_sizes.tsv` | Cluster size and further-refinement eligibility at each level |

### FST clumping outputs

| File | Description |
|---|---|
| `{method}.pairwise_fst.tsv.gz` | Pairwise Hudson FST estimates at every recalculation iteration |
| `{method}.fst_merge_history.tsv` | Ordered merge history with pre-merge FST and cluster sizes |
| `{method}.final_membership.tsv.gz` | Original hierarchical membership plus `pre_fst_cluster` and `final_cluster` |
| `{method}.fst_clumping_summary.tsv` | Initial/final cluster counts, merges, threshold, and sample count |

### Final IBD-sharing outputs

`FINAL_IBD_SHARING` uses each method's post-clumping `final_cluster`
assignments and streams through the genome-wide pair table once.

For a pair of different clusters \(A\) and \(B\), the heatmap value is:

\[
\frac{\text{total observed IBD between }A\text{ and }B}{n_A n_B}.
\]

For the diagonal of cluster \(A\), it is:

\[
\frac{\text{total observed IBD within }A}{\binom{n_A}{2}}.
\]

Pairs absent from the genome-wide edge table are therefore included as zero
without being materialized. This makes heatmap cells comparable across
clusters of different sizes.

| File | Description |
|---|---|
| `{method}.ibd_sharing.tsv` | Long-format within- and between-cluster totals, possible/detected pair counts, detected-pair fraction, all-pair mean, and detected-pair mean |
| `{method}.ibd_mean_matrix.tsv` | Symmetric matrix of mean IBD length per possible pair, including implicit zeros |
| `{method}.ibd_heatmap.png` | Annotated heatmap of the matrix |
| `{method}.within_cluster_distribution.tsv` | Detected within-cluster pair count and mean, minimum, quartiles, median, and maximum IBD length |
| `{method}.within_ibd_boxplot.png` | Distribution of total IBD length among detected within-cluster sample pairs; black dots show means and `n` labels show cluster sample sizes |

The heatmap and boxplot deliberately use different denominators. Heatmap means
include all possible pairs, including zero-sharing pairs. Boxplots describe the
distribution conditional on a within-cluster pair having detected IBD sharing.

### Final Hudson FST outputs

`FINAL_FST_HEATMAP` reruns PLINK 2 once per enabled clustering method using the
post-clumping `final_cluster` assignments. This explicitly estimates Hudson
\(F_{ST}\) for the final merged clusters rather than extracting values from an
earlier clumping iteration.

| File | Description |
|---|---|
| `{method}.final_fst.tsv` | Long-format exact pairwise Hudson FST estimates |
| `{method}.final_fst_matrix.tsv` | Symmetric final-cluster FST matrix |
| `{method}.final_fst_heatmap.png` | Annotated `YlOrRd` heatmap with adaptive text color |

### Final-cluster IBDNe outputs

| File | Description |
|---|---|
| `{method}.ibdne_status.tsv` | Cluster sample size, retained segment count, explicit `npairs`, completion/skip status, and failure message |
| `Cluster_*.ne` | Effective population-size estimate and 95% bootstrap confidence interval by generation |
| `Cluster_*.boot` | Original and bootstrap effective population-size histories |
| `Cluster_*.pair.excl` | Close sample pairs excluded by IBDNe |
| `Cluster_*.region.excl` | Genomic regions excluded by IBDNe |
| `Cluster_*.log` | Native IBDNe run log |
| `Cluster_*.command.txt` | Exact IBDNe command and parameter values used for the cluster |
| `Cluster_*.driver.log` | Standard output/error captured by the pipeline driver |



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
