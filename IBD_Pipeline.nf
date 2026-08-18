/*
 * Hap-IBD segment-detection pipeline
 *
 * Author: Justin Pelletier
 * Version: 2.2
 */

nextflow.enable.dsl = 2


/*
 * Filter variants according to MAF and variant missingness.
 *
 * The process also verifies that:
 * - variants remain after QC;
 * - the number of samples remains unchanged;
 * - VCF and genetic-map chromosome identifiers match;
 * - every retained genotype is phased and non-missing.
 */
process QC_VCF {
    tag "chr${chromosome}"

    cpus params.resources.qc.cpus
    memory params.resources.qc.memory
    time params.resources.qc.time

    errorStrategy 'terminate'

    beforeScript """
    module load bcftools
    """

    publishDir "${params.outdir}/qc",
        pattern: '*.qc.stats.txt',
        mode: 'copy',
        overwrite: true

    publishDir "${params.outdir}/qc",
        pattern: '*.qc.summary.tsv',
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chromosome),
        path(vcf),
        path(vcf_index),
        path(genetic_map),
        path(gap_file)

    output:
    tuple val(chromosome),
          path("chr${chromosome}.qc.vcf.gz"),
          path("chr${chromosome}.qc.vcf.gz.tbi"),
          path(genetic_map),
          path(gap_file),
          emit: vcfs

    tuple val(chromosome),
          path("chr${chromosome}.samples.txt"),
          emit: samples

    path "chr${chromosome}.qc.stats.txt",
          emit: stats

    path "chr${chromosome}.qc.summary.tsv",
          emit: summary

    script:
    """
    set -euo pipefail


    # ------------------------------------------------------------------
    # Count variants and samples before QC.
    # ------------------------------------------------------------------

    variants_before_qc=\$(
        bcftools index \
            -n \
            ${vcf}
    )

    samples_before_qc=\$(
        bcftools query \
            -l \
            ${vcf} |
        wc -l
    )

    if ! [[ "\${variants_before_qc}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to determine the number of variants before QC." >&2
        echo "Observed value: \${variants_before_qc}" >&2
        exit 1
    fi

    if ! [[ "\${samples_before_qc}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to determine the number of samples before QC." >&2
        echo "Observed value: \${samples_before_qc}" >&2
        exit 1
    fi

    if [[ \${variants_before_qc} -eq 0 ]]
    then
        echo "ERROR: the chromosome ${chromosome} input VCF contains no variants." >&2
        exit 1
    fi

    if [[ \${samples_before_qc} -eq 0 ]]
    then
        echo "ERROR: the chromosome ${chromosome} input VCF contains no samples." >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Apply variant-level QC.
    # ------------------------------------------------------------------

    bcftools +fill-tags \
        ${vcf} \
        -Ou \
        -- \
        -t MAF,F_MISSING |
    bcftools view \
        -i 'MAF>=${params.qc.min_maf} && F_MISSING<=${params.qc.max_variant_missingness}' \
        -Oz \
        -o chr${chromosome}.qc.vcf.gz

    tabix \
        -f \
        -p vcf \
        chr${chromosome}.qc.vcf.gz


    # ------------------------------------------------------------------
    # Count variants and samples after QC.
    # ------------------------------------------------------------------

    variants_after_qc=\$(
        bcftools index \
            -n \
            chr${chromosome}.qc.vcf.gz
    )

    samples_after_qc=\$(
        bcftools query \
            -l \
            chr${chromosome}.qc.vcf.gz |
        wc -l
    )

    if ! [[ "\${variants_after_qc}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to determine the number of variants after QC." >&2
        echo "Observed value: \${variants_after_qc}" >&2
        exit 1
    fi

    if ! [[ "\${samples_after_qc}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to determine the number of samples after QC." >&2
        echo "Observed value: \${samples_after_qc}" >&2
        exit 1
    fi

    if [[ \${variants_after_qc} -eq 0 ]]
    then
        echo "ERROR: no variants remained after QC on chromosome ${chromosome}." >&2
        exit 1
    fi

    if [[ \${samples_after_qc} -eq 0 ]]
    then
        echo "ERROR: no samples remained after QC on chromosome ${chromosome}." >&2
        exit 1
    fi

    if [[ \${samples_before_qc} -ne \${samples_after_qc} ]]
    then
        echo "ERROR: the number of samples changed during variant-level QC." >&2
        echo "Before QC: \${samples_before_qc}" >&2
        echo "After QC:  \${samples_after_qc}" >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Calculate variant-retention statistics.
    # ------------------------------------------------------------------

    variants_removed=\$(
        awk \
            -v before="\${variants_before_qc}" \
            -v after="\${variants_after_qc}" \
            'BEGIN { print before - after }'
    )

    retention_percent=\$(
        awk \
            -v before="\${variants_before_qc}" \
            -v after="\${variants_after_qc}" \
            '
            BEGIN {
                printf "%.2f", 100 * after / before
            }
            '
    )


    # ------------------------------------------------------------------
    # Write the complete sample list.
    # ------------------------------------------------------------------

    bcftools query \
        -l \
        chr${chromosome}.qc.vcf.gz \
        > chr${chromosome}.samples.txt

    if [[ ! -s chr${chromosome}.samples.txt ]]
    then
        echo "ERROR: the chromosome ${chromosome} sample list is empty." >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Verify VCF and genetic-map chromosome identifiers.
    # ------------------------------------------------------------------

    vcf_chromosome=\$(
        bcftools index \
            -s \
            chr${chromosome}.qc.vcf.gz |
        awk 'NR == 1 { print \$1 }'
    )

    map_chromosome=\$(
        awk '
            NF > 0 && \$1 !~ /^#/ {
                print \$1
                exit
            }
        ' ${genetic_map}
    )

    if [[ -z "\${vcf_chromosome}" ]]
    then
        echo "ERROR: unable to read the chromosome identifier from the QC-filtered VCF." >&2
        exit 1
    fi

    if [[ -z "\${map_chromosome}" ]]
    then
        echo "ERROR: unable to read the chromosome identifier from the genetic map." >&2
        exit 1
    fi

    if [[ "\${vcf_chromosome}" != "\${map_chromosome}" ]]
    then
        echo "ERROR: chromosome identifiers do not match for chromosome ${chromosome}." >&2
        echo "VCF chromosome: \${vcf_chromosome}" >&2
        echo "Map chromosome: \${map_chromosome}" >&2
        echo "Set params.chr_in_chrom_field to match the VCF convention." >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Verify that every retained genotype is phased and non-missing.
    #
    # AWK scans the complete stream rather than exiting early. This
    # prevents bcftools from receiving SIGPIPE under `set -o pipefail`.
    # ------------------------------------------------------------------

    if bcftools query \
        -f '[%GT\\n]' \
        chr${chromosome}.qc.vcf.gz |
    awk '
        index(\$0, "/") || index(\$0, ".") {
            invalid = 1
        }

        END {
            exit !invalid
        }
    '
    then
        echo "ERROR: chr${chromosome}.qc.vcf.gz contains an unphased or missing genotype." >&2
        echo "Hap-IBD requires every genotype to be phased and non-missing." >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Write the concise QC summary.
    # ------------------------------------------------------------------

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        chromosome \
        samples \
        variants_before_QC \
        variants_after_QC \
        variants_removed \
        retention_percent \
        min_maf \
        max_variant_missingness \
        > chr${chromosome}.qc.summary.tsv

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        "\${vcf_chromosome}" \
        "\${samples_after_qc}" \
        "\${variants_before_qc}" \
        "\${variants_after_qc}" \
        "\${variants_removed}" \
        "\${retention_percent}" \
        "${params.qc.min_maf}" \
        "${params.qc.max_variant_missingness}" \
        >> chr${chromosome}.qc.summary.tsv


    # ------------------------------------------------------------------
    # Write the detailed bcftools statistics report.
    # ------------------------------------------------------------------

    bcftools stats \
        chr${chromosome}.qc.vcf.gz \
        > chr${chromosome}.qc.stats.txt
    """
}

/*
 * Confirm that sample IDs and their ordering are identical across
 * all chromosome-specific VCF files.
 *
 * cohort.samples.txt is also used to include participants without
 * detected IBD edges as isolated clustering vertices.
 */
process VALIDATE_SAMPLES {
    tag 'validate sample lists'

    cpus 1
    memory '1 GB'
    time '30m'

    errorStrategy 'terminate'

    publishDir "${params.outdir}/qc",
        pattern: 'cohort.samples.txt',
        mode: 'copy',
        overwrite: true

    input:
    path sample_files

    output:
    path 'cohort.samples.txt'

    script:
    """
    set -euo pipefail

    cp \
        ${sample_files[0]} \
        cohort.samples.txt

    for sample_file in ${sample_files.join(' ')}
    do
        if ! cmp -s cohort.samples.txt \${sample_file}
        then
            echo "ERROR: sample IDs or sample ordering differ between chromosomes." >&2
            echo "Reference sample list: ${sample_files[0]}" >&2
            echo "Non-matching sample list: \${sample_file}" >&2
            exit 1
        fi
    done
    """
}


/*
 * Detect IBD segments using Hap-IBD.
 *
 * Hap-IBD directly interpolates genetic positions from the supplied
 * chromosome-specific genetic map. A marker-matched PLINK map is
 * therefore not required.
 *
 * If excluded-region filtering is enabled, complete IBD segments
 * overlapping at least one listed interval are removed.
 */
process HAP_IBD {
    tag "chr${chromosome}"

    cpus params.resources.hapibd.cpus
    memory params.resources.hapibd.memory
    time params.resources.hapibd.time

    errorStrategy 'retry'
    maxRetries 1

    beforeScript """
    module load java/25.36
    module load bcftools
    module load bedtools
    """

    publishDir "${params.outdir}/segments",
        pattern: '*.hapibd.ibd.gz',
        mode: 'copy',
        overwrite: true

    publishDir "${params.outdir}/logs",
        pattern: '*.hapibd.log',
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chromosome),
          path(vcf),
          path(vcf_index),
          path(genetic_map),
          path(gap_file)

    path hapibd_jar

    output:
    tuple val(chromosome),
          path("chr${chromosome}.hapibd.ibd.gz"),
          emit: segments

    path "chr${chromosome}.hapibd.log",
          emit: logs

    script:
    def javaGb = Math.max(
        1,
        (task.memory.giga * 0.90) as int
    )

    def filterGaps = params.remove_gaps \
        ? """
          vcf_chromosome=\$(
              bcftools index -s ${vcf} |
              awk 'NR == 1 { print \$1 }'
          )

          if [[ -z "\${vcf_chromosome}" ]]
          then
              echo "ERROR: unable to determine the chromosome name from ${vcf}." >&2
              exit 1
          fi

          normalized_chromosome=\$(
              echo "\${vcf_chromosome}" |
              sed 's/^chr//'
          )

          # Select the excluded regions for the current chromosome.
          #
          # The chromosome label is rewritten to match the VCF and
          # Hap-IBD output, allowing either "1" or "chr1" conventions.
          awk \
              -v target="\${vcf_chromosome}" \
              -v normalized="\${normalized_chromosome}" \
              '
              BEGIN {
                  OFS = "\\t"
              }

              {
                  bed_chromosome = \$1
                  sub(/^chr/, "", bed_chromosome)

                  if (bed_chromosome == normalized) {
                      print target, \$2, \$3, \$4
                  }
              }
              ' ${gap_file} \
              > chr${chromosome}.normalized.gaps.bed

          if [[ ! -s chr${chromosome}.normalized.gaps.bed ]]
          then
              echo "WARNING: no excluded regions were found for chromosome \${vcf_chromosome}." >&2

              cp \
                  chr${chromosome}.raw.ibd \
                  chr${chromosome}.filtered.ibd
          else
              # Hap-IBD reports one-based inclusive positions.
              # BED uses zero-based, half-open coordinates.
              awk '
                  BEGIN {
                      OFS = "\\t"
                  }

                  {
                      bed_start = \$6 - 1

                      if (bed_start < 0) {
                          bed_start = 0
                      }

                      print \$5, bed_start, \$7, \$0
                  }
              ' chr${chromosome}.raw.ibd |
                  bedtools intersect \
                      -v \
                      -a - \
                      -b chr${chromosome}.normalized.gaps.bed |
                  cut -f4- \
                  > chr${chromosome}.filtered.ibd
          fi
          """ \
        : """
          cp \
              chr${chromosome}.raw.ibd \
              chr${chromosome}.filtered.ibd
          """

    """
    set -euo pipefail

    java \
        -Xmx${javaGb}g \
        -jar ${hapibd_jar} \
        gt=${vcf} \
        map=${genetic_map} \
        out=chr${chromosome}.raw \
        min-seed=${params.hapibd.min_seed_cm} \
        min-extend=${params.hapibd.min_extend_cm} \
        min-output=${params.hapibd.min_output_cm} \
        min-markers=${params.hapibd.min_markers} \
        min-mac=${params.hapibd.min_mac} \
        max-gap=${params.hapibd.max_gap_bp} \
        nthreads=${task.cpus}

    mv \
        chr${chromosome}.raw.log \
        chr${chromosome}.hapibd.log

    gunzip -c \
        chr${chromosome}.raw.ibd.gz \
        > chr${chromosome}.raw.ibd

    ${filterGaps}

    bgzip -c \
        chr${chromosome}.filtered.ibd \
        > chr${chromosome}.hapibd.ibd.gz
    """
}



/*
 * Generate pair-level IBD summary statistics for each chromosome.
 */
process PER_PAIR_CHROMOSOME {
    tag "chr${chromosome}"

    cpus params.resources.summary.cpus
    memory params.resources.summary.memory
    time params.resources.summary.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/per_pair/by_chromosome",
        pattern: '*.per_pair.tsv.gz',
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chromosome),
          path(ibd_file)

    path summary_script

    output:
    tuple val(chromosome),
          path("chr${chromosome}.per_pair.tsv.gz")

    script:
    """
    set -euo pipefail

    if [[ -z "\${SLURM_TMPDIR:-}" ]]
    then
        echo "ERROR: SLURM_TMPDIR is not defined." >&2
        exit 1
    fi

    python3 ${summary_script} \
        --input ${ibd_file} \
        --output chr${chromosome}.per_pair.tsv.gz \
        --threshold-cm ${params.summary.segment_threshold_cm} \
        --temporary-directory "\${SLURM_TMPDIR}"
    """
}


/*
 * Combine chromosome-specific pair summaries into a genome-wide
 * pair-level summary.
 *
 * This avoids rereading all raw Hap-IBD segment files.
 */
process PER_PAIR_GENOMEWIDE {
    tag 'genome-wide summary'

    cpus params.resources.summary.cpus
    memory params.resources.summary.memory
    time params.resources.summary.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/per_pair",
        pattern: 'genomewide.per_pair.tsv.gz',
        mode: 'copy',
        overwrite: true

    input:
    path chromosome_summaries
    path summary_script

    output:
    path 'genomewide.per_pair.tsv.gz'

    script:
    """
    set -euo pipefail

    if [[ -z "\${SLURM_TMPDIR:-}" ]]
    then
        echo "ERROR: SLURM_TMPDIR is not defined." >&2
        exit 1
    fi

    python3 ${summary_script} \
        --input ${chromosome_summaries.join(' ')} \
        --output genomewide.per_pair.tsv.gz \
        --threshold-cm ${params.summary.segment_threshold_cm} \
        --aggregate-summaries \
        --temporary-directory "\${SLURM_TMPDIR}"
    """
}


/*
 * Apply recursive Louvain or Leiden clustering to the genome-wide
 * weighted IBD-sharing graph.
 */
process CLUSTER_GRAPH {
    tag "${method} clustering"

    cpus params.resources.clustering.cpus
    memory params.resources.clustering.memory
    time params.resources.clustering.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/clustering",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(method),
          val(max_levels)

    path pair_summary
    path samples
    path clustering_script

    output:
    tuple val(method),
          path("${method}.membership.tsv.gz"),
          emit: membership

    path "${method}.selected_membership.tsv.gz"
    path "${method}.diagnostics.tsv"
    path "${method}.cluster_sizes.tsv"

    script:
    """
    set -euo pipefail

    python3 ${clustering_script} \
        --edges ${pair_summary} \
        --samples ${samples} \
        --method ${method} \
        --max-levels ${max_levels} \
        --weight-column ${params.clustering.weight_column} \
        --min-cluster-size ${params.clustering.min_cluster_size} \
        --min-modularity-gain ${params.clustering.min_modularity_gain} \
        --resolution ${params.clustering.leiden_resolution} \
        --seed ${params.clustering.seed} \
        --auto-select ${params.clustering.auto_select} \
        --output-prefix ${method}
    """
}


/*
 * Prepare an autosomal, common, high-quality, LD-pruned PLINK 2
 * dataset for genotype-based FST estimation.
 *
 * The chromosome VCFs have already passed the Hap-IBD QC step. This
 * process applies a separate, configurable FST-specific filter because
 * FST estimation and graph construction have different requirements.
 * The concatenated BCF is written to node-local storage and removed at
 * the end of the task.
 */
process PREPARE_FST_DATA {
    tag 'prepare FST genotype data'

    cpus params.resources.fst_prepare.cpus
    memory params.resources.fst_prepare.memory
    time params.resources.fst_prepare.time

    errorStrategy 'terminate'

    beforeScript """
    module load bcftools
    module load ${params.plink2_module}
    """

    publishDir "${params.outdir}/fst/genotypes",
        mode: 'copy',
        overwrite: true

    input:
    path qc_vcfs

    output:
    tuple path('fst.pruned.pgen'),
          path('fst.pruned.pvar.zst'),
          path('fst.pruned.psam'),
          emit: dataset

    path 'fst.prune.in',
         emit: variants

    path 'fst.prepare.summary.tsv',
         emit: summary

    path 'fst.prepare.log',
         emit: log

    script:
    """
    set -euo pipefail

    if [[ -z "\${SLURM_TMPDIR:-}" ]]
    then
        echo "ERROR: SLURM_TMPDIR is not defined." >&2
        exit 1
    fi

    merged_bcf="\${SLURM_TMPDIR}/fst.autosomes.bcf"

    bcftools concat \
        --output-type b \
        --output "\${merged_bcf}" \
        ${qc_vcfs.join(' ')}

    bcftools index \
        --force \
        "\${merged_bcf}"

    plink2 \
        --bcf "\${merged_bcf}" \
        --autosome \
        --snps-only just-acgt \
        --max-alleles 2 \
        --set-all-var-ids '@:#:\$r:\$a' \
        --new-id-max-allele-len 100 truncate \
        --maf ${params.fst.min_maf} \
        --geno ${params.fst.max_variant_missingness} \
        --hwe ${params.fst.hwe_pvalue} midp \
        --make-pgen vzs \
        --sort-vars \
        --threads ${task.cpus} \
        --memory ${task.memory.mega as int} \
        --out fst.qc \
        > fst.prepare.log 2>&1

    plink2 \
        --pfile fst.qc vzs \
        --indep-pairwise \
            ${params.fst.prune_window_kb}kb \
            ${params.fst.prune_step_variants} \
            ${params.fst.prune_r2} \
        --threads ${task.cpus} \
        --memory ${task.memory.mega as int} \
        --out fst \
        >> fst.prepare.log 2>&1

    if [[ ! -s fst.prune.in ]]
    then
        echo "ERROR: LD pruning retained no variants." >&2
        exit 1
    fi

    plink2 \
        --pfile fst.qc vzs \
        --extract fst.prune.in \
        --make-pgen vzs \
        --threads ${task.cpus} \
        --memory ${task.memory.mega as int} \
        --out fst.pruned \
        >> fst.prepare.log 2>&1

    sample_count=\$(awk 'NR > 1 { count++ } END { print count + 0 }' fst.pruned.psam)
    variant_count=\$(wc -l < fst.prune.in)

    printf 'samples\tvariants\tmin_maf\tmax_variant_missingness\thwe_pvalue\tprune_window_kb\tprune_step_variants\tprune_r2\n' \
        > fst.prepare.summary.tsv

    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
        "\${sample_count}" \
        "\${variant_count}" \
        '${params.fst.min_maf}' \
        '${params.fst.max_variant_missingness}' \
        '${params.fst.hwe_pvalue}' \
        '${params.fst.prune_window_kb}' \
        '${params.fst.prune_step_variants}' \
        '${params.fst.prune_r2}' \
        >> fst.prepare.summary.tsv

    rm -f \
        "\${merged_bcf}" \
        "\${merged_bcf}.csi"
    """
}


/*
 * Iteratively clump terminal graph communities using Hudson FST.
 *
 * The Python driver writes the current categorical cluster phenotype,
 * calls PLINK 2, merges the eligible pair with the smallest FST below
 * the configured threshold, and repeats until no eligible pair remains.
 */
process FST_CLUMP {
    tag "${method} FST clumping"

    cpus params.resources.fst_clump.cpus
    memory params.resources.fst_clump.memory
    time params.resources.fst_clump.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    module load ${params.plink2_module}
    source ${params.python_venv}/bin/activate
    """

    publishDir {
        "${params.outdir}/clustering/${method}"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(method),
          path(membership)

    tuple path(pgen),
          path(pvar),
          path(psam)

    path fst_clump_script

    output:
    path "${method}.pairwise_fst.tsv.gz"
    path "${method}.fst_merge_history.tsv"
    path "${method}.final_membership.tsv.gz"
    path "${method}.fst_clumping_summary.tsv"

    script:
    """
    set -euo pipefail

    python3 ${fst_clump_script} \
        --pgen ${pgen} \
        --pvar ${pvar} \
        --psam ${psam} \
        --membership ${membership} \
        --cluster-column ${params.fst.cluster_column} \
        --threshold ${params.fst.clump_threshold} \
        --min-cluster-size ${params.fst.min_cluster_size} \
        --plink2 plink2 \
        --threads ${task.cpus} \
        --memory-mb ${task.memory.mega as int} \
        --output-prefix ${method}
    """
}


workflow {
    /*
     * Validate configuration parameters.
     */
    if (!params.input_pattern) {
        error 'Set params.input_pattern.'
    }

    if (!params.input_index_pattern) {
        error 'Set params.input_index_pattern.'
    }

    if (!params.hapibd_jar) {
        error 'Set params.hapibd_jar.'
    }

    if (!params.per_pair_script) {
        error 'Set params.per_pair_script.'
    }

    if (!params.python_venv) {
        error 'Set params.python_venv.'
    }

    if (params.fst.enabled && !params.plink2_module) {
        error 'Set params.plink2_module when FST clumping is enabled.'
    }

    if (params.fst.enabled && !params.fst_clump_script) {
        error 'Set params.fst_clump_script when FST clumping is enabled.'
    }

    if (
        params.fst.enabled &&
        !(params.Louvain || params.Leiden)
    ) {
        error 'Enable Louvain and/or Leiden when FST clumping is enabled.'
    }

    if (
        (params.Louvain || params.Leiden) &&
        !params.clustering_script
    ) {
        error 'Set params.clustering_script when clustering is enabled.'
    }

    if (!params.gap_file) {
        error 'Set params.gap_file.'
    }

    if (!(params.chromosomes as List)) {
        error 'params.chromosomes must contain at least one chromosome.'
    }

    if (!(params.chr_in_chrom_field instanceof Boolean)) {
        error 'params.chr_in_chrom_field must be true or false.'
    }

    if (
        params.qc.min_maf < 0 ||
        params.qc.min_maf > 0.5
    ) {
        error 'params.qc.min_maf must be between 0 and 0.5.'
    }

    if (
        params.qc.max_variant_missingness < 0 ||
        params.qc.max_variant_missingness > 1
    ) {
        error 'params.qc.max_variant_missingness must be between 0 and 1.'
    }

    if (params.hapibd.min_seed_cm <= 0) {
        error 'params.hapibd.min_seed_cm must be greater than zero.'
    }

    if (
        params.hapibd.min_extend_cm <= 0 ||
        params.hapibd.min_extend_cm > params.hapibd.min_seed_cm
    ) {
        error 'params.hapibd.min_extend_cm must be greater than zero and no greater than min_seed_cm.'
    }

    if (params.hapibd.min_output_cm <= 0) {
        error 'params.hapibd.min_output_cm must be greater than zero.'
    }

    if (params.hapibd.min_markers < 1) {
        error 'params.hapibd.min_markers must be at least one.'
    }

    if (params.hapibd.min_mac < 1) {
        error 'params.hapibd.min_mac must be at least one.'
    }

    if (params.hapibd.max_gap_bp < -1) {
        error 'params.hapibd.max_gap_bp must be at least -1.'
    }

    if (params.summary.segment_threshold_cm < 0) {
        error 'params.summary.segment_threshold_cm must be greater than or equal to zero.'
    }

    if (
        params.n_louvain < 0 ||
        params.n_leiden < 0
    ) {
        error 'Clustering refinement depths must be greater than or equal to zero.'
    }

    if (
        (params.Louvain || params.Leiden) &&
        params.clustering.min_cluster_size < 2
    ) {
        error 'params.clustering.min_cluster_size must be at least two.'
    }

    if (
        (params.Louvain || params.Leiden) &&
        params.clustering.min_modularity_gain < 0
    ) {
        error 'params.clustering.min_modularity_gain must be greater than or equal to zero.'
    }

    if (
        params.Leiden &&
        params.clustering.leiden_resolution <= 0
    ) {
        error 'params.clustering.leiden_resolution must be greater than zero.'
    }

    if (
        params.fst.enabled &&
        (params.fst.min_maf < 0 || params.fst.min_maf > 0.5)
    ) {
        error 'params.fst.min_maf must be between 0 and 0.5.'
    }

    if (
        params.fst.enabled &&
        (params.fst.max_variant_missingness < 0 || params.fst.max_variant_missingness > 1)
    ) {
        error 'params.fst.max_variant_missingness must be between 0 and 1.'
    }

    if (
        params.fst.enabled &&
        (params.fst.hwe_pvalue <= 0 || params.fst.hwe_pvalue > 1)
    ) {
        error 'params.fst.hwe_pvalue must be greater than zero and no greater than one.'
    }

    if (params.fst.enabled && params.fst.prune_window_kb < 1) {
        error 'params.fst.prune_window_kb must be at least one.'
    }

    if (params.fst.enabled && params.fst.prune_step_variants < 1) {
        error 'params.fst.prune_step_variants must be at least one.'
    }

    if (
        params.fst.enabled &&
        (params.fst.prune_r2 <= 0 || params.fst.prune_r2 >= 1)
    ) {
        error 'params.fst.prune_r2 must be greater than zero and smaller than one.'
    }

    if (params.fst.enabled && params.fst.clump_threshold < 0) {
        error 'params.fst.clump_threshold must be greater than or equal to zero.'
    }

    if (params.fst.enabled && params.fst.min_cluster_size < 2) {
        error 'params.fst.min_cluster_size must be at least two.'
    }

    if (params.fst.enabled && !params.fst.cluster_column) {
        error 'Set params.fst.cluster_column.'
    }


    /*
     * Resolve fixed program and asset files.
     */
    hapibdJar = file(
        params.hapibd_jar,
        checkIfExists: true
    )

    summaryScript = file(
        params.per_pair_script,
        checkIfExists: true
    )

    excludedRegions = file(
        params.gap_file,
        checkIfExists: true
    )

    if (params.fst.enabled) {
        fstClumpScript = file(
            params.fst_clump_script,
            checkIfExists: true
        )
    }


    /*
     * Select the bundled GRCh38 PLINK maps without modifying or
     * renaming them.
     *
     * false:
     *   VCF chromosomes are 1, 2, ..., 22.
     *   Maps are:
     *   assets/no_chr_in_chrom_field/plink.chrN.GRCh38.map
     *
     * true:
     *   VCF chromosomes are chr1, chr2, ..., chr22.
     *   Maps are:
     *   assets/chr_in_chrom_field/plink.chrchrN.GRCh38.map
     */
    def chrInChromField = params.chr_in_chrom_field

    def geneticMapDirectory = chrInChromField \
        ? file("${projectDir}/assets/chr_in_chrom_field")
        : file("${projectDir}/assets/no_chr_in_chrom_field")

    def geneticMapPrefix = chrInChromField \
        ? 'plink.chrchr'
        : 'plink.chr'


    /*
     * Construct one input tuple for each chromosome.
     */
    inputs = Channel
        .fromList(params.chromosomes as List)
        .map { chromosome ->
            def chromosomeString = chromosome as String

            def vcf = file(
                params.input_pattern.replace(
                    '{chr}',
                    chromosomeString
                ),
                checkIfExists: true
            )

            def vcfIndex = file(
                params.input_index_pattern.replace(
                    '{chr}',
                    chromosomeString
                ),
                checkIfExists: true
            )

            def geneticMap = file(
                "${geneticMapDirectory}/${geneticMapPrefix}${chromosomeString}.GRCh38.map",
                checkIfExists: true
            )

            tuple(
                chromosome,
                vcf,
                vcfIndex,
                geneticMap,
                excludedRegions
            )
        }


    /*
     * Perform variant QC and verify chromosome sample lists.
     */
    qc = QC_VCF(
        inputs
    )

    validatedSamples = VALIDATE_SAMPLES(
        qc.samples
            .map { chromosome, sampleFile ->
                sampleFile
            }
            .collect()
    )


    /*
     * Prepare one reusable, LD-pruned genotype dataset for all enabled
     * clustering methods. This branch can run concurrently with Hap-IBD.
     */
    if (params.fst.enabled) {
        fstDataset = PREPARE_FST_DATA(
            qc.vcfs
                .toSortedList { first, second ->
                    (first[0] as Integer) <=> (second[0] as Integer)
                }
                .map { chromosomeTuples ->
                    chromosomeTuples.collect { chromosomeTuple ->
                        chromosomeTuple[1]
                    }
                }
        )
    }


    /*
     * Detect IBD segments using the official GRCh38 PLINK maps.
     */
    hapIBD = HAP_IBD(
        qc.vcfs,
        hapibdJar
    )


    /*
     * Generate chromosome-specific and genome-wide pair summaries.
     */
    chromosomePairSummaries = PER_PAIR_CHROMOSOME(
        hapIBD.segments,
        summaryScript
    )

    genomewideSummary = PER_PAIR_GENOMEWIDE(
        chromosomePairSummaries
            .map { chromosome, summaryFile ->
                summaryFile
            }
            .collect(),
        summaryScript
    )


    /*
     * Build the requested clustering-method channel.
     */
    clusteringMethods = []

    if (params.Louvain) {
        clusteringMethods << tuple(
            'louvain',
            params.n_louvain as Integer
        )
    }

    if (params.Leiden) {
        clusteringMethods << tuple(
            'leiden',
            params.n_leiden as Integer
        )
    }


    /*
     * Run each requested clustering algorithm.
     *
     * genomewideSummary and validatedSamples are reusable value
     * channels because their upstream inputs were created with collect().
     */
    if (clusteringMethods) {
        clusteringScript = file(
            params.clustering_script,
            checkIfExists: true
        )

        clusteringResults = CLUSTER_GRAPH(
            Channel.fromList(clusteringMethods),
            genomewideSummary,
            validatedSamples,
            clusteringScript
        )

        if (params.fst.enabled) {
            FST_CLUMP(
                clusteringResults.membership,
                fstDataset.dataset,
                fstClumpScript
            )
        }
    }
}
