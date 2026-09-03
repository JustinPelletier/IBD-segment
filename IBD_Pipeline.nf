/*
 * Hap-IBD segment-detection, relatedness-filtering and clustering pipeline
 *
 * Author: Justin Pelletier
 * Version: 2.4
 */

nextflow.enable.dsl = 2


/*
 * Apply variant-level QC independently to each chromosome.
 *
 * The process:
 * - filters variants by MAF and missingness;
 * - confirms variants and samples remain;
 * - confirms that sample count is unchanged;
 * - confirms agreement between VCF and genetic-map chromosome labels;
 * - confirms that every retained genotype is phased and non-missing;
 * - writes the chromosome sample list.
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
        exit 1
    fi

    if ! [[ "\${samples_before_qc}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to determine the number of samples before QC." >&2
        exit 1
    fi

    if [[ \${variants_before_qc} -eq 0 ]]
    then
        echo "ERROR: chromosome ${chromosome} contains no variants." >&2
        exit 1
    fi

    if [[ \${samples_before_qc} -eq 0 ]]
    then
        echo "ERROR: chromosome ${chromosome} contains no samples." >&2
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
        exit 1
    fi

    if ! [[ "\${samples_after_qc}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to determine the number of samples after QC." >&2
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
        echo "ERROR: sample count changed during variant-level QC." >&2
        echo "Before QC: \${samples_before_qc}" >&2
        echo "After QC:  \${samples_after_qc}" >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Write the complete chromosome sample list.
    # ------------------------------------------------------------------
    bcftools query \
        -l \
        chr${chromosome}.qc.vcf.gz \
        > chr${chromosome}.samples.txt

    if [[ ! -s chr${chromosome}.samples.txt ]]
    then
        echo "ERROR: chromosome ${chromosome} sample list is empty." >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Confirm that VCF and genetic-map chromosome labels agree.
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
        echo "ERROR: unable to determine the VCF chromosome label." >&2
        exit 1
    fi

    if [[ -z "\${map_chromosome}" ]]
    then
        echo "ERROR: unable to determine the genetic-map chromosome label." >&2
        exit 1
    fi

    if [[ "\${vcf_chromosome}" != "\${map_chromosome}" ]]
    then
        echo "ERROR: VCF and genetic-map chromosome labels differ." >&2
        echo "VCF chromosome: \${vcf_chromosome}" >&2
        echo "Map chromosome: \${map_chromosome}" >&2
        echo "Set params.chr_in_chrom_field to match the VCF convention." >&2
        exit 1
    fi


    # ------------------------------------------------------------------
    # Confirm that all retained genotypes are phased and non-missing.
    # ------------------------------------------------------------------

    invalid_genotypes=\$(
        bcftools query \
            -f '[%GT\\n]' \
            chr${chromosome}.qc.vcf.gz |
        awk '
            {
                if (index(\$0, "/") > 0 || index(\$0, ".") > 0) {
                    invalid++
                }
            }

            END {
                print invalid + 0
            }
        '
    )

    if ! [[ "\${invalid_genotypes}" =~ ^[0-9]+\$ ]]
    then
        echo "ERROR: unable to validate retained genotypes." >&2
        exit 1
    fi

    if [[ \${invalid_genotypes} -ne 0 ]]
    then
        echo "ERROR: chromosome ${chromosome} contains \${invalid_genotypes} unphased or missing genotypes after QC." >&2
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
    # Write the QC summary.
    # ------------------------------------------------------------------

    printf 'chromosome\\tvariants_before_qc\\tvariants_after_qc\\tvariants_removed\\tretention_percent\\tsamples\\tmin_maf\\tmax_variant_missingness\\n' \
        > chr${chromosome}.qc.summary.tsv

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        '${chromosome}' \
        "\${variants_before_qc}" \
        "\${variants_after_qc}" \
        "\${variants_removed}" \
        "\${retention_percent}" \
        "\${samples_after_qc}" \
        '${params.qc.min_maf}' \
        '${params.qc.max_variant_missingness}' \
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
 * Validate the complete-cohort chromosome sample lists.
 *
 * This process is called once, before relatedness discovery.
 */
process VALIDATE_SAMPLES {
    tag 'validate full-cohort sample lists'

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
        if ! cmp -s cohort.samples.txt "\${sample_file}"
        then
            echo "ERROR: full-cohort sample IDs or ordering differ between chromosomes." >&2
            echo "Reference: ${sample_files[0]}" >&2
            echo "Different: \${sample_file}" >&2
            exit 1
        fi
    done

    if [[ ! -s cohort.samples.txt ]]
    then
        echo "ERROR: validated full-cohort sample list is empty." >&2
        exit 1
    fi
    """
}


/*
 * Validate chromosome sample lists after unrelated-sample subsetting.
 *
 * A separate process name is required because a DSL2 process cannot be
 * invoked twice from the same workflow.
 */
process VALIDATE_ANALYSIS_SAMPLES {
    tag 'validate unrelated-cohort sample lists'

    cpus 1
    memory '1 GB'
    time '30m'

    errorStrategy 'terminate'

    publishDir "${params.outdir}/qc/unrelated",
        pattern: 'analysis.samples.txt',
        mode: 'copy',
        overwrite: true

    input:
    path sample_files

    output:
    path 'analysis.samples.txt'

    script:
    """
    set -euo pipefail


    cp \
        ${sample_files[0]} \
        analysis.samples.txt

    for sample_file in ${sample_files.join(' ')}
    do
        if ! cmp -s analysis.samples.txt "\${sample_file}"
        then
            echo "ERROR: unrelated-cohort sample IDs or ordering differ between chromosomes." >&2
            echo "Reference: ${sample_files[0]}" >&2
            echo "Different: \${sample_file}" >&2
            exit 1
        fi
    done

    if [[ ! -s analysis.samples.txt ]]
    then
        echo "ERROR: validated unrelated-cohort sample list is empty." >&2
        exit 1
    fi
    """
}


/*
 * Preliminary full-cohort Hap-IBD analysis for relatedness discovery.
 *
 * This process uses params.relatedness_hapibd. The final unrelated-cohort
 * Hap-IBD process will separately use params.hapibd.
 */
process HAP_IBD_RELATEDNESS {
    tag "relatedness chr${chromosome}"

    cpus params.resources.relatedness_hapibd.cpus
    memory params.resources.relatedness_hapibd.memory
    time params.resources.relatedness_hapibd.time

    errorStrategy 'retry'
    maxRetries 1

    beforeScript """
    module load java/25.36
    module load bcftools
    module load bedtools
    """

    publishDir "${params.outdir}/relatedness/discovery_segments",
        pattern: '*.relatedness.hapibd.ibd.gz',
        mode: 'copy',
        overwrite: true

    publishDir "${params.outdir}/relatedness/logs",
        pattern: '*.relatedness.hapibd.log',
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
          path("chr${chromosome}.relatedness.hapibd.ibd.gz"),
          emit: segments

    path "chr${chromosome}.relatedness.hapibd.log",
         emit: logs

    script:
    def javaGb = Math.max(
        1,
        (task.memory.giga * 0.90) as int
    )

    def filterGaps = params.remove_gaps \
        ? """
          vcf_chromosome=\$(
              bcftools index \
                  -s \
                  ${vcf} |
              awk 'NR == 1 { print \$1 }'
          )

          if [[ -z "\${vcf_chromosome}" ]]
          then
              echo "ERROR: unable to determine chromosome from ${vcf}." >&2
              exit 1
          fi

          normalized_chromosome=\$(
              echo "\${vcf_chromosome}" |
              sed 's/^chr//'
          )

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
              > chr${chromosome}.relatedness.gaps.bed

          if [[ -s chr${chromosome}.relatedness.gaps.bed ]]
          then
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
              ' chr${chromosome}.relatedness.raw.ibd |
              bedtools intersect \
                  -v \
                  -a - \
                  -b chr${chromosome}.relatedness.gaps.bed |
              cut -f4- \
                  > chr${chromosome}.relatedness.filtered.ibd
          else
              cp \
                  chr${chromosome}.relatedness.raw.ibd \
                  chr${chromosome}.relatedness.filtered.ibd
          fi
          """ \
        : """
          cp \
              chr${chromosome}.relatedness.raw.ibd \
              chr${chromosome}.relatedness.filtered.ibd
          """

    """
    set -euo pipefail

    java \
        -Xmx${javaGb}g \
        -jar ${hapibd_jar} \
        gt=${vcf} \
        map=${genetic_map} \
        out=chr${chromosome}.relatedness.raw \
        min-seed=${params.relatedness_hapibd.min_seed_cm} \
        min-extend=${params.relatedness_hapibd.min_extend_cm} \
        min-output=${params.relatedness_hapibd.min_output_cm} \
        min-markers=${params.relatedness_hapibd.min_markers} \
        min-mac=${params.relatedness_hapibd.min_mac} \
        max-gap=${params.relatedness_hapibd.max_gap_bp} \
        nthreads=${task.cpus}

    mv \
        chr${chromosome}.relatedness.raw.log \
        chr${chromosome}.relatedness.hapibd.log

    gunzip -c \
        chr${chromosome}.relatedness.raw.ibd.gz \
        > chr${chromosome}.relatedness.raw.ibd

    ${filterGaps}

    bgzip -c \
        chr${chromosome}.relatedness.filtered.ibd \
        > chr${chromosome}.relatedness.hapibd.ibd.gz
    """
}


/*
 * Calculate chromosome-level pair summaries for the preliminary run.
 */
process PER_PAIR_CHROMOSOME_RELATEDNESS {
    tag "relatedness chr${chromosome}"

    cpus params.resources.summary.cpus
    memory params.resources.summary.memory
    time params.resources.summary.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    input:
    tuple val(chromosome),
          path(ibd_file)

    path summary_script

    output:
    tuple val(chromosome),
          path("chr${chromosome}.relatedness.per_pair.tsv.gz")

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
        --output chr${chromosome}.relatedness.per_pair.tsv.gz \
        --threshold-cm ${params.relatedness_summary.segment_threshold_cm} \
        --temporary-directory "\${SLURM_TMPDIR}"
    """
}


/*
 * Aggregate preliminary chromosome summaries genome-wide.
 *
 * The selector uses total_IBD_length_cM from this table.
 */
process PER_PAIR_GENOMEWIDE_RELATEDNESS {
    tag 'relatedness genome-wide summary'

    cpus params.resources.summary.cpus
    memory params.resources.summary.memory
    time params.resources.summary.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/relatedness",
        pattern: 'relatedness.genomewide.per_pair.tsv.gz',
        mode: 'copy',
        overwrite: true

    input:
    path chromosome_summaries
    path summary_script

    output:
    path 'relatedness.genomewide.per_pair.tsv.gz'

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
        --output relatedness.genomewide.per_pair.tsv.gz \
        --threshold-cm ${params.relatedness_summary.segment_threshold_cm} \
        --aggregate-summaries \
        --temporary-directory "\${SLURM_TMPDIR}"
    """
}


/*
 * Prepare common LD-pruned SNPs and calculate KING kinship coefficients.
 */
process KING_RELATEDNESS {
    tag 'KING relatedness discovery'

    cpus params.resources.relatedness_king.cpus
    memory params.resources.relatedness_king.memory
    time params.resources.relatedness_king.time

    errorStrategy 'terminate'

    beforeScript """
    module load bcftools
    module load ${params.plink2_module}
    """

    publishDir "${params.outdir}/relatedness",
        mode: 'copy',
        overwrite: true

    input:
    path qc_vcfs

    output:
    path 'king.kin0',
         emit: pairs

    path 'king.prune.in',
         emit: variants

    path 'king.relatedness.summary.tsv',
         emit: summary

    path 'king.relatedness.log',
         emit: log

    script:
    """
    set -euo pipefail

    if [[ -z "\${SLURM_TMPDIR:-}" ]]
    then
        echo "ERROR: SLURM_TMPDIR is not defined." >&2
        exit 1
    fi

    merged_bcf="\${SLURM_TMPDIR}/king.autosomes.bcf"

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
        --maf ${params.relatedness.min_maf} \
        --geno ${params.relatedness.max_variant_missingness} \
        --make-pgen vzs \
        --sort-vars \
        --threads ${task.cpus} \
        --memory ${task.memory.mega as int} \
        --out king.qc \
        > king.relatedness.log 2>&1

    plink2 \
        --pfile king.qc vzs \
        --indep-pairwise \
            ${params.relatedness.prune_window_variants} \
            ${params.relatedness.prune_step_variants} \
            ${params.relatedness.prune_r2} \
        --indep-order ${params.relatedness.indep_order} \
        --threads ${task.cpus} \
        --memory ${task.memory.mega as int} \
        --out king \
        >> king.relatedness.log 2>&1

    if [[ ! -s king.prune.in ]]
    then
        echo "ERROR: KING LD pruning retained no variants." >&2
        exit 1
    fi

    plink2 \
        --pfile king.qc vzs \
        --extract king.prune.in \
        --make-king-table \
        --king-table-filter ${params.relatedness.king_cutoff} \
        --threads ${task.cpus} \
        --memory ${task.memory.mega as int} \
        --out king \
        >> king.relatedness.log 2>&1

    if [[ ! -f king.kin0 ]]
    then
        echo "ERROR: PLINK 2 did not create king.kin0." >&2
        exit 1
    fi

    sample_count=\$(
        awk '
            NR > 1 {
                count++
            }

            END {
                print count + 0
            }
        ' king.qc.psam
    )

    variant_count=\$(
        wc -l \
            < king.prune.in
    )

    printf 'samples\\tvariants\\tmin_maf\\tmax_variant_missingness\\tking_cutoff\\tprune_window_variants\\tprune_step_variants\\tprune_r2\\tindep_order\\n' \
        > king.relatedness.summary.tsv

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
        "\${sample_count}" \
        "\${variant_count}" \
        '${params.relatedness.min_maf}' \
        '${params.relatedness.max_variant_missingness}' \
        '${params.relatedness.king_cutoff}' \
        '${params.relatedness.prune_window_variants}' \
        '${params.relatedness.prune_step_variants}' \
        '${params.relatedness.prune_r2}' \
        '${params.relatedness.indep_order}' \
        >> king.relatedness.summary.tsv

    rm -f \
        "\${merged_bcf}" \
        "\${merged_bcf}.csi"
    """
}


/*
 * Combine KING and Hap-IBD related pairs and select the final unrelated set.
 *
 * removal_mode:
 * - optimal: retain a maximum/practical maximal mutually unrelated subset;
 * - remove_all: exclude every participant belonging to a flagged pair.
 */
process SELECT_UNRELATED_SAMPLES {
    tag 'select unrelated samples'

    cpus params.resources.relatedness_select.cpus
    memory params.resources.relatedness_select.memory
    time params.resources.relatedness_select.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/relatedness",
        mode: 'copy',
        overwrite: true

    input:
    path cohort_samples
    path king_pairs
    path genomewide_summary
    path selector_script

    output:
    path 'unrelated.samples.txt',
         emit: unrelated

    path 'related_samples_removed.txt',
         emit: removed

    path 'related_pairs.tsv.gz',
         emit: pairs

    path 'related_components.tsv.gz',
         emit: components

    path 'relatedness_selection_summary.tsv',
         emit: summary

    script:
    """
    set -euo pipefail

    python3 ${selector_script} \
        --samples ${cohort_samples} \
        --king ${king_pairs} \
        --ibd ${genomewide_summary} \
        --king-cutoff ${params.relatedness.king_cutoff} \
        --ibd-total-cm-cutoff ${params.relatedness.ibd_total_cm_cutoff} \
        --removal-mode ${params.relatedness.removal_mode} \
        --output-prefix relatedness

    mv \
        relatedness.unrelated.samples.txt \
        unrelated.samples.txt

    mv \
        relatedness.related_samples_removed.txt \
        related_samples_removed.txt

    mv \
        relatedness.related_pairs.tsv.gz \
        related_pairs.tsv.gz

    mv \
        relatedness.related_components.tsv.gz \
        related_components.tsv.gz

    mv \
        relatedness.selection_summary.tsv \
        relatedness_selection_summary.tsv
    """
}


/*
 * Subset each QC VCF to the selected unrelated cohort.
 *
 * This occurs before the final Hap-IBD and FST branches.
 */
process APPLY_UNRELATED_SAMPLES {
    tag "unrelated chr${chromosome}"

    cpus params.resources.qc.cpus
    memory params.resources.qc.memory
    time params.resources.qc.time

    errorStrategy 'terminate'

    beforeScript """
    module load bcftools
    """

    publishDir "${params.outdir}/qc/unrelated",
        pattern: '*.unrelated.summary.tsv',
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chromosome),
          path(vcf),
          path(vcf_index),
          path(genetic_map),
          path(gap_file)

    path unrelated_samples

    output:
    tuple val(chromosome),
          path("chr${chromosome}.unrelated.vcf.gz"),
          path("chr${chromosome}.unrelated.vcf.gz.tbi"),
          path(genetic_map),
          path(gap_file),
          emit: vcfs

    tuple val(chromosome),
          path("chr${chromosome}.samples.txt"),
          emit: samples

    path "chr${chromosome}.unrelated.summary.tsv",
         emit: summary

    script:
    """
    set -euo pipefail

    if [[ ! -s ${unrelated_samples} ]]
    then
        echo "ERROR: the unrelated sample list is empty." >&2
        exit 1
    fi

    duplicate_count=\$(
        sort ${unrelated_samples} |
        uniq -d |
        wc -l
    )

    if [[ \${duplicate_count} -ne 0 ]]
    then
        echo "ERROR: the unrelated sample list contains duplicate IDs." >&2
        exit 1
    fi

    bcftools query \
        -l \
        ${vcf} \
        > input.samples.txt

    missing_count=\$(
        comm \
            -23 \
            <(sort ${unrelated_samples}) \
            <(sort input.samples.txt) |
        wc -l
    )

    if [[ \${missing_count} -ne 0 ]]
    then
        echo "ERROR: \${missing_count} unrelated-list IDs are absent from chromosome ${chromosome}." >&2

        comm \
            -23 \
            <(sort ${unrelated_samples}) \
            <(sort input.samples.txt) >&2

        exit 1
    fi

    bcftools view \
        --samples-file ${unrelated_samples} \
        --output-type z \
        --output chr${chromosome}.unrelated.vcf.gz \
        ${vcf}

    tabix \
        -f \
        -p vcf \
        chr${chromosome}.unrelated.vcf.gz

    bcftools query \
        -l \
        chr${chromosome}.unrelated.vcf.gz \
        > chr${chromosome}.samples.txt

    input_count=\$(wc -l < input.samples.txt)
    requested_count=\$(wc -l < ${unrelated_samples})
    retained_count=\$(wc -l < chr${chromosome}.samples.txt)

    if [[ \${retained_count} -ne \${requested_count} ]]
    then
        echo "ERROR: retained sample count does not match the unrelated list." >&2
        exit 1
    fi

    printf 'chromosome\\tinput_samples\\tretained_samples\\tremoved_samples\\n' \
        > chr${chromosome}.unrelated.summary.tsv

    printf '%s\\t%s\\t%s\\t%s\\n' \
        '${chromosome}' \
        "\${input_count}" \
        "\${retained_count}" \
        "\$((input_count - retained_count))" \
        >> chr${chromosome}.unrelated.summary.tsv
    """
}



/*
 * Final Hap-IBD analysis on the unrelated cohort.
 *
 * This process uses params.hapibd, independently of the preliminary
 * relatedness_hapibd parameters.
 */
process HAP_IBD {
    tag "final chr${chromosome}"

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
              bcftools index \
                  -s \
                  ${vcf} |
              awk 'NR == 1 { print \$1 }'
          )

          if [[ -z "\${vcf_chromosome}" ]]
          then
              echo "ERROR: unable to determine chromosome from ${vcf}." >&2
              exit 1
          fi

          normalized_chromosome=\$(
              echo "\${vcf_chromosome}" |
              sed 's/^chr//'
          )

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
              > chr${chromosome}.gaps.bed

          if [[ -s chr${chromosome}.gaps.bed ]]
          then
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
                  -b chr${chromosome}.gaps.bed |
              cut -f4- \
                  > chr${chromosome}.filtered.ibd
          else
              cp \
                  chr${chromosome}.raw.ibd \
                  chr${chromosome}.filtered.ibd
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
 * Final chromosome-level pair summaries.
 */
process PER_PAIR_CHROMOSOME {
    tag "final chr${chromosome}"

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
 * Final genome-wide pair summary used for graph clustering.
 */
process PER_PAIR_GENOMEWIDE {
    tag 'final genome-wide summary'

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
 * Apply recursive Louvain or Leiden clustering to the final unrelated-cohort
 * weighted IBD graph.
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
 * Prepare a common, autosomal, biallelic and LD-pruned PGEN dataset for FST.
 *
 * Input VCFs have already been restricted to the unrelated cohort.
 */
process PREPARE_FST_DATA {
    tag 'prepare unrelated-cohort FST data'

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
        echo "ERROR: FST LD pruning retained no variants." >&2
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

    if [[ ! -s fst.pruned.pgen ]]
    then
        echo "ERROR: PLINK 2 did not create fst.pruned.pgen." >&2
        exit 1
    fi

    if [[ ! -s fst.pruned.pvar.zst ]]
    then
        echo "ERROR: PLINK 2 did not create fst.pruned.pvar.zst." >&2
        exit 1
    fi

    if [[ ! -s fst.pruned.psam ]]
    then
        echo "ERROR: PLINK 2 did not create fst.pruned.psam." >&2
        exit 1
    fi

    sample_count=\$(
        awk '
            NR > 1 {
                count++
            }

            END {
                print count + 0
            }
        ' fst.pruned.psam
    )

    variant_count=\$(
        wc -l \
            < fst.prune.in
    )

    printf 'samples\\tvariants\\tmin_maf\\tmax_variant_missingness\\thwe_pvalue\\tprune_window_kb\\tprune_step_variants\\tprune_r2\\n' \
        > fst.prepare.summary.tsv

    printf '%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\t%s\\n' \
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
 * Iteratively merge terminal clusters using Hudson FST.
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
    tuple val(method),
          path("${method}.final_membership.tsv.gz"),
          emit: final_membership
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


/*
 * Summarize within- and between-cluster IBD sharing after FST merging.
 * Undetected pairs are included as zero through analytical denominators.
 */
process FINAL_IBD_SHARING {
    tag "${method} final IBD sharing"

    cpus params.resources.final_ibd.cpus
    memory params.resources.final_ibd.memory
    time params.resources.final_ibd.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    source ${params.python_venv}/bin/activate
    """

    publishDir {
        "${params.outdir}/clustering/${method}/ibd_sharing"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(method),
          path(membership)

    path pair_summary
    path analysis_script

    output:
    path "${method}.ibd_sharing.tsv"
    path "${method}.within_cluster_distribution.tsv"
    path "${method}.ibd_mean_matrix.tsv"
    path "${method}.ibd_heatmap.png"
    path "${method}.within_ibd_boxplot.png"

    script:
    """
    set -euo pipefail

    python3 ${analysis_script} \
        --edges ${pair_summary} \
        --membership ${membership} \
        --weight-column ${params.clustering.weight_column} \
        --output-prefix ${method}
    """
}


/*
 * Recompute Hudson FST from the post-merging assignments and plot the
 * final pairwise matrix.
 */
process FINAL_FST_HEATMAP {
    tag "${method} final FST heatmap"

    cpus params.resources.final_fst.cpus
    memory params.resources.final_fst.memory
    time params.resources.final_fst.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    module load ${params.plink2_module}
    source ${params.python_venv}/bin/activate
    """

    publishDir {
        "${params.outdir}/clustering/${method}/final_fst"
    },
        mode: 'copy',
        overwrite: true

    input:
    tuple val(method),
          path(membership)

    tuple path(pgen),
          path(pvar),
          path(psam)

    path heatmap_script

    output:
    path "${method}.final_fst.tsv"
    path "${method}.final_fst_matrix.tsv"
    path "${method}.final_fst_heatmap.png"

    script:
    """
    set -euo pipefail

    python3 ${heatmap_script} \
        --pgen ${pgen} \
        --pvar ${pvar} \
        --psam ${psam} \
        --membership ${membership} \
        --plink2 plink2 \
        --threads ${task.cpus} \
        --memory-mb ${task.memory.mega as int} \
        --output-prefix ${method}
    """
}







/*
 * Estimate recent effective population-size histories for eligible final
 * clusters. The driver streams the final Hap-IBD segments once, retains only
 * within-cluster segments, and runs cluster-level IBDNe analyses concurrently.
 */
process RUN_IBDNE_FINAL_CLUSTERS {
    tag "${method} final-cluster IBDNe"

    cpus params.resources.ibdne.cpus
    memory params.resources.ibdne.memory
    time params.resources.ibdne.time

    errorStrategy 'terminate'

    beforeScript """
    module load python/3.14.2
    module load java/25.36
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

    path ibd_segments
    path genetic_maps
    path ibdne_jar
    path ibdne_runner_script

    output:
    tuple val(method),
          path('ibdne'),
          emit: results

    script:
    """
    set -euo pipefail

    if [[ -z "\${SLURM_TMPDIR:-}" ]]
    then
        echo "ERROR: SLURM_TMPDIR is not defined." >&2
        exit 1
    fi

    temporary_directory="\${SLURM_TMPDIR}/${method}.ibdne.\${SLURM_JOB_ID:-task}"
    mkdir -p "\${temporary_directory}" ibdne

    python3 ${ibdne_runner_script} \
        --segments ${ibd_segments.join(' ')} \
        --maps ${genetic_maps.join(' ')} \
        --membership ${membership} \
        --ibdne-jar ${ibdne_jar} \
        --method ${method} \
        --minimum-cluster-size ${params.ibdne.min_cluster_size} \
        --mincm ${params.ibdne.mincm} \
        --nits ${params.ibdne.nits} \
        --nboots ${params.ibdne.nboots} \
        --filtersamples ${params.ibdne.filtersamples} \
        --seed ${params.ibdne.seed} \
        --threads ${task.cpus} \
        --parallel-clusters ${params.ibdne.parallel_clusters} \
        --java-heap-gb ${params.ibdne.java_heap_gb} \
        --temporary-directory "\${temporary_directory}" \
        --output-directory ibdne
    """
}


workflow {
    /*
     * ---------------------------------------------------------------
     * General configuration validation
     * ---------------------------------------------------------------
     */

    if (!params.input_pattern) {
        error 'Set params.input_pattern.'
    }

    if (!params.input_index_pattern) {
        error 'Set params.input_index_pattern.'
    }

    if (!(params.chromosomes as List)) {
        error 'params.chromosomes must contain at least one chromosome.'
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

    if (!params.gap_file) {
        error 'Set params.gap_file.'
    }

    if (!(params.chr_in_chrom_field instanceof Boolean)) {
        error 'params.chr_in_chrom_field must be true or false.'
    }


    /*
     * ---------------------------------------------------------------
     * Relatedness mode
     * ---------------------------------------------------------------
     *
     * discover_apply:
     *   Run KING and preliminary full-cohort Hap-IBD, generate the
     *   unrelated list, and then run the final analyses.
     *
     * reuse:
     *   Skip relatedness discovery and use an existing unrelated list.
     *
     * disabled:
     *   Run final analyses on the complete cohort.
     */

    def relatednessMode = params.relatedness.mode as String

    if (!(relatednessMode in ['discover_apply', 'reuse', 'disabled'])) {
        error """
        params.relatedness.mode must be one of:
        - discover_apply
        - reuse
        - disabled
        """.stripIndent().trim()
    }

    if (
        relatednessMode == 'discover_apply' &&
        !params.relatedness_selector_script
    ) {
        error 'Set params.relatedness_selector_script for discover_apply mode.'
    }

    if (
        relatednessMode == 'reuse' &&
        !params.relatedness.sample_list
    ) {
        error 'Set params.relatedness.sample_list for reuse mode.'
    }

    if (
        (relatednessMode == 'discover_apply' || params.fst.enabled) &&
        !params.plink2_module
    ) {
        error 'Set params.plink2_module when relatedness discovery or FST is enabled.'
    }


    /*
     * ---------------------------------------------------------------
     * Variant-QC validation
     * ---------------------------------------------------------------
     */

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


    /*
     * ---------------------------------------------------------------
     * Relatedness-discovery validation
     * ---------------------------------------------------------------
     */

    if (relatednessMode == 'discover_apply') {
        if (!(params.relatedness.removal_mode in ['optimal', 'remove_all'])) {
            error "params.relatedness.removal_mode must be 'optimal' or 'remove_all'."
        }

        if (
            params.relatedness.king_cutoff <= 0 ||
            params.relatedness.king_cutoff >= 0.5
        ) {
            error 'params.relatedness.king_cutoff must be greater than 0 and smaller than 0.5.'
        }

        if (params.relatedness.ibd_total_cm_cutoff <= 0) {
            error 'params.relatedness.ibd_total_cm_cutoff must be greater than zero.'
        }

        if (
            params.relatedness.min_maf < 0 ||
            params.relatedness.min_maf > 0.5
        ) {
            error 'params.relatedness.min_maf must be between 0 and 0.5.'
        }

        if (
            params.relatedness.max_variant_missingness < 0 ||
            params.relatedness.max_variant_missingness > 1
        ) {
            error 'params.relatedness.max_variant_missingness must be between 0 and 1.'
        }

        if (params.relatedness.prune_window_variants < 2) {
            error 'params.relatedness.prune_window_variants must be at least two.'
        }

        if (params.relatedness.prune_step_variants < 1) {
            error 'params.relatedness.prune_step_variants must be at least one.'
        }

        if (
            params.relatedness.prune_r2 <= 0 ||
            params.relatedness.prune_r2 >= 1
        ) {
            error 'params.relatedness.prune_r2 must be greater than 0 and smaller than 1.'
        }

        if (!(params.relatedness.indep_order in [1, 2])) {
            error 'params.relatedness.indep_order must be 1 or 2.'
        }

        if (params.relatedness_hapibd.min_seed_cm <= 0) {
            error 'params.relatedness_hapibd.min_seed_cm must be greater than zero.'
        }

        if (
            params.relatedness_hapibd.min_extend_cm <= 0 ||
            params.relatedness_hapibd.min_extend_cm >
                params.relatedness_hapibd.min_seed_cm
        ) {
            error """
            params.relatedness_hapibd.min_extend_cm must be greater than zero
            and no greater than relatedness_hapibd.min_seed_cm.
            """.stripIndent().trim()
        }

        if (params.relatedness_hapibd.min_output_cm <= 0) {
            error 'params.relatedness_hapibd.min_output_cm must be greater than zero.'
        }

        if (params.relatedness_hapibd.min_markers < 1) {
            error 'params.relatedness_hapibd.min_markers must be at least one.'
        }

        if (params.relatedness_hapibd.min_mac < 1) {
            error 'params.relatedness_hapibd.min_mac must be at least one.'
        }

        if (params.relatedness_hapibd.max_gap_bp < -1) {
            error 'params.relatedness_hapibd.max_gap_bp must be at least -1.'
        }

        if (params.relatedness_summary.segment_threshold_cm < 0) {
            error 'params.relatedness_summary.segment_threshold_cm must be at least zero.'
        }
    }


    /*
     * ---------------------------------------------------------------
     * Final Hap-IBD validation
     * ---------------------------------------------------------------
     */

    if (params.hapibd.min_seed_cm <= 0) {
        error 'params.hapibd.min_seed_cm must be greater than zero.'
    }

    if (
        params.hapibd.min_extend_cm <= 0 ||
        params.hapibd.min_extend_cm > params.hapibd.min_seed_cm
    ) {
        error """
        params.hapibd.min_extend_cm must be greater than zero and no greater
        than params.hapibd.min_seed_cm.
        """.stripIndent().trim()
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
        error 'params.summary.segment_threshold_cm must be at least zero.'
    }


    /*
     * ---------------------------------------------------------------
     * Clustering validation
     * ---------------------------------------------------------------
     */

    if (
        params.n_louvain < 0 ||
        params.n_leiden < 0
    ) {
        error 'Clustering refinement depths must be at least zero.'
    }

    if (
        (params.Louvain || params.Leiden) &&
        !params.clustering_script
    ) {
        error 'Set params.clustering_script when clustering is enabled.'
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
        error 'params.clustering.min_modularity_gain must be at least zero.'
    }

    if (
        params.Leiden &&
        params.clustering.leiden_resolution <= 0
    ) {
        error 'params.clustering.leiden_resolution must be greater than zero.'
    }


    /*
     * ---------------------------------------------------------------
     * FST validation
     * ---------------------------------------------------------------
     */

    if (
        params.fst.enabled &&
        !(params.Louvain || params.Leiden)
    ) {
        error 'Enable Louvain and/or Leiden when FST clumping is enabled.'
    }

    if (
        params.fst.enabled &&
        !params.fst_clump_script
    ) {
        error 'Set params.fst_clump_script when FST clumping is enabled.'
    }

    if (
        params.fst.enabled &&
        !params.final_ibd_script
    ) {
        error 'Set params.final_ibd_script when FST clumping is enabled.'
    }

    if (
        params.fst.enabled &&
        !params.final_fst_heatmap_script
    ) {
        error 'Set params.final_fst_heatmap_script when FST clumping is enabled.'
    }

    if (params.ibdne.enabled && !params.fst.enabled) {
        error 'Enable params.fst.enabled when final-cluster IBDNe is enabled.'
    }

    if (params.ibdne.enabled && !params.ibdne_jar) {
        error 'Set params.ibdne_jar when final-cluster IBDNe is enabled.'
    }

    if (params.ibdne.enabled && !params.ibdne_runner_script) {
        error 'Set params.ibdne_runner_script when final-cluster IBDNe is enabled.'
    }

    if (params.ibdne.enabled && params.ibdne.min_cluster_size < 2) {
        error 'params.ibdne.min_cluster_size must be at least two.'
    }

    if (params.ibdne.enabled && params.ibdne.mincm <= 0) {
        error 'params.ibdne.mincm must be greater than zero.'
    }

    if (
        params.ibdne.enabled &&
        (params.ibdne.nits < 1 || params.ibdne.nboots < 1)
    ) {
        error 'params.ibdne.nits and params.ibdne.nboots must be positive.'
    }

    if (
        params.ibdne.enabled &&
        (params.ibdne.parallel_clusters < 1 || params.ibdne.java_heap_gb < 1)
    ) {
        error 'IBDNe parallel_clusters and java_heap_gb must be positive.'
    }

    if (
        params.fst.enabled &&
        (
            params.fst.min_maf < 0 ||
            params.fst.min_maf > 0.5
        )
    ) {
        error 'params.fst.min_maf must be between 0 and 0.5.'
    }

    if (
        params.fst.enabled &&
        (
            params.fst.max_variant_missingness < 0 ||
            params.fst.max_variant_missingness > 1
        )
    ) {
        error 'params.fst.max_variant_missingness must be between 0 and 1.'
    }

    if (
        params.fst.enabled &&
        (
            params.fst.hwe_pvalue <= 0 ||
            params.fst.hwe_pvalue > 1
        )
    ) {
        error 'params.fst.hwe_pvalue must be greater than 0 and no greater than 1.'
    }

    if (
        params.fst.enabled &&
        params.fst.prune_window_kb < 1
    ) {
        error 'params.fst.prune_window_kb must be at least one.'
    }

    if (
        params.fst.enabled &&
        params.fst.prune_step_variants < 1
    ) {
        error 'params.fst.prune_step_variants must be at least one.'
    }

    if (
        params.fst.enabled &&
        (
            params.fst.prune_r2 <= 0 ||
            params.fst.prune_r2 >= 1
        )
    ) {
        error 'params.fst.prune_r2 must be greater than 0 and smaller than 1.'
    }

    if (
        params.fst.enabled &&
        params.fst.clump_threshold < 0
    ) {
        error 'params.fst.clump_threshold must be at least zero.'
    }

    if (
        params.fst.enabled &&
        params.fst.min_cluster_size < 2
    ) {
        error 'params.fst.min_cluster_size must be at least two.'
    }

    if (
        params.fst.enabled &&
        !params.fst.cluster_column
    ) {
        error 'Set params.fst.cluster_column.'
    }


    /*
     * ---------------------------------------------------------------
     * Resolve scripts and fixed files
     * ---------------------------------------------------------------
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

    if (params.Louvain || params.Leiden) {
        clusteringScript = file(
            params.clustering_script,
            checkIfExists: true
        )
    }

    if (params.fst.enabled) {
        fstClumpScript = file(
            params.fst_clump_script,
            checkIfExists: true
        )

        finalIbdScript = file(
            params.final_ibd_script,
            checkIfExists: true
        )

        finalFstHeatmapScript = file(
            params.final_fst_heatmap_script,
            checkIfExists: true
        )
    }

    if (params.ibdne.enabled) {
        ibdneJar = file(
            params.ibdne_jar,
            checkIfExists: true
        )

        ibdneRunnerScript = file(
            params.ibdne_runner_script,
            checkIfExists: true
        )
    }

    if (relatednessMode == 'discover_apply') {
        relatednessSelectorScript = file(
            params.relatedness_selector_script,
            checkIfExists: true
        )
    }

    if (relatednessMode == 'reuse') {
        suppliedUnrelatedSamples = file(
            params.relatedness.sample_list,
            checkIfExists: true
        )
    }


    /*
     * ---------------------------------------------------------------
     * Genetic-map path convention
     * ---------------------------------------------------------------
     */

    def chrInChromField = params.chr_in_chrom_field

    def geneticMapDirectory = chrInChromField \
        ? file("${projectDir}/assets/chr_in_chrom_field")
        : file("${projectDir}/assets/no_chr_in_chrom_field")

    def geneticMapPrefix = chrInChromField \
        ? 'plink.chrchr'
        : 'plink.chr'


    /*
     * ---------------------------------------------------------------
     * Construct chromosome input tuples
     * ---------------------------------------------------------------
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
     * ---------------------------------------------------------------
     * Full-cohort variant QC
     * ---------------------------------------------------------------
     */

    qc = QC_VCF(
        inputs
    )

    fullCohortSamples = VALIDATE_SAMPLES(
        qc.samples
            .map { chromosome, sampleFile ->
                sampleFile
            }
            .collect()
    )


    /*
     * ---------------------------------------------------------------
     * Relatedness discovery or sample-list reuse
     * ---------------------------------------------------------------
     */

    if (relatednessMode == 'discover_apply') {
        relatednessHapIBD = HAP_IBD_RELATEDNESS(
            qc.vcfs,
            hapibdJar
        )

        relatednessChromosomeSummaries =
            PER_PAIR_CHROMOSOME_RELATEDNESS(
                relatednessHapIBD.segments,
                summaryScript
            )

        relatednessGenomewideSummary =
            PER_PAIR_GENOMEWIDE_RELATEDNESS(
                relatednessChromosomeSummaries
                    .map { chromosome, summaryFile ->
                        summaryFile
                    }
                    .collect(),
                summaryScript
            )

        kingResults = KING_RELATEDNESS(
            qc.vcfs
                .toSortedList { first, second ->
                    (first[0] as Integer) <=>
                    (second[0] as Integer)
                }
                .map { chromosomeTuples ->
                    chromosomeTuples.collect { chromosomeTuple ->
                        chromosomeTuple[1]
                    }
                }
        )

        relatednessSelection = SELECT_UNRELATED_SAMPLES(
            fullCohortSamples,
            kingResults.pairs,
            relatednessGenomewideSummary,
            relatednessSelectorScript
        )

        /*
         * collect() converts the selector output into a reusable value
         * channel so the same list is supplied to every chromosome.
         */
        selectedSamples = relatednessSelection.unrelated
            .collect()
            .map { sampleLists ->
                sampleLists[0]
            }
    }
    else if (relatednessMode == 'reuse') {
        selectedSamples = suppliedUnrelatedSamples
    }


    /*
     * ---------------------------------------------------------------
     * Select the cohort used by all final analyses
     * ---------------------------------------------------------------
     */

    if (relatednessMode != 'disabled') {
        unrelatedQc = APPLY_UNRELATED_SAMPLES(
            qc.vcfs,
            selectedSamples
        )

        analysisVcfs = unrelatedQc.vcfs
        analysisSampleFiles = unrelatedQc.samples
    }
    else {
        analysisVcfs = qc.vcfs
        analysisSampleFiles = qc.samples
    }

    analysisSamples = VALIDATE_ANALYSIS_SAMPLES(
        analysisSampleFiles
            .map { chromosome, sampleFile ->
                sampleFile
            }
            .collect()
    )


    /*
     * ---------------------------------------------------------------
     * Final unrelated-cohort Hap-IBD analysis
     * ---------------------------------------------------------------
     */

    finalHapIBD = HAP_IBD(
        analysisVcfs,
        hapibdJar
    )

    finalChromosomeSummaries = PER_PAIR_CHROMOSOME(
        finalHapIBD.segments,
        summaryScript
    )

    finalGenomewideSummary = PER_PAIR_GENOMEWIDE(
        finalChromosomeSummaries
            .map { chromosome, summaryFile ->
                summaryFile
            }
            .collect(),
        summaryScript
    )


    /*
     * ---------------------------------------------------------------
     * Prepare the unrelated-cohort FST genotype dataset
     * ---------------------------------------------------------------
     */

    if (params.fst.enabled) {
        fstDataset = PREPARE_FST_DATA(
            analysisVcfs
                .toSortedList { first, second ->
                    (first[0] as Integer) <=>
                    (second[0] as Integer)
                }
                .map { chromosomeTuples ->
                    chromosomeTuples.collect { chromosomeTuple ->
                        chromosomeTuple[1]
                    }
                }
        )
    }


    /*
     * ---------------------------------------------------------------
     * Louvain and Leiden clustering
     * ---------------------------------------------------------------
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

    if (clusteringMethods) {
        clusteringResults = CLUSTER_GRAPH(
            Channel.fromList(clusteringMethods),
            finalGenomewideSummary,
            analysisSamples,
            clusteringScript
        )

        /*
         * FST clumping uses the same unrelated sample set because both the
         * membership and PGEN files descend from analysisVcfs.
         */
        if (params.fst.enabled) {
            fstClumpResults = FST_CLUMP(
                clusteringResults.membership,
                fstDataset.dataset,
                fstClumpScript
            )

            FINAL_IBD_SHARING(
                fstClumpResults.final_membership,
                finalGenomewideSummary,
                finalIbdScript
            )

            FINAL_FST_HEATMAP(
                fstClumpResults.final_membership,
                fstDataset.dataset,
                finalFstHeatmapScript
            )

            if (params.ibdne.enabled) {
                finalIbdSegmentsForIbdne = finalHapIBD.segments
                    .toSortedList { first, second ->
                        (first[0] as Integer) <=>
                        (second[0] as Integer)
                    }
                    .map { chromosomeTuples ->
                        chromosomeTuples.collect { chromosomeTuple ->
                            chromosomeTuple[1]
                        }
                    }

                finalMapsForIbdne = analysisVcfs
                    .toSortedList { first, second ->
                        (first[0] as Integer) <=>
                        (second[0] as Integer)
                    }
                    .map { chromosomeTuples ->
                        chromosomeTuples.collect { chromosomeTuple ->
                            chromosomeTuple[3]
                        }
                    }

                RUN_IBDNE_FINAL_CLUSTERS(
                    fstClumpResults.final_membership,
                    finalIbdSegmentsForIbdne,
                    finalMapsForIbdne,
                    ibdneJar,
                    ibdneRunnerScript
                )
            }
        }
    }
}
