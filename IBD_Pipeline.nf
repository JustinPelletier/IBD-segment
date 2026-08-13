/*
 * Hap-IBD segment-detection pipeline
 *
 * Author: Justin Pelletier
 * Version: 2.1
 */

nextflow.enable.dsl = 2


/*
 * Filter variants according to MAF and variant missingness.
 *
 * The process also verifies that all retained genotypes are phased
 * and non-missing, as required by Hap-IBD.
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

    input:
    tuple val(chromosome),
          path(vcf),
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

    script:
    """
    set -euo pipefail

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

    bcftools query \
        -l chr${chromosome}.qc.vcf.gz \
        > chr${chromosome}.samples.txt

    if [[ ! -s chr${chromosome}.samples.txt ]]
    then
        echo "ERROR: no samples were found in the chromosome ${chromosome} VCF." >&2
        exit 1
    fi

    number_of_variants=\$(
        bcftools index \
            -n chr${chromosome}.qc.vcf.gz
    )

    if [[ \${number_of_variants} -eq 0 ]]
    then
        echo "ERROR: no variants remained after QC on chromosome ${chromosome}." >&2
        exit 1
    fi

    # Hap-IBD requires identical chromosome identifiers in the VCF
    # and the first column of the PLINK genetic map.
    vcf_chromosome=\$(
        bcftools index -s chr${chromosome}.qc.vcf.gz |
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

    if [[ -z "\${vcf_chromosome}" || -z "\${map_chromosome}" ]]
    then
        echo "ERROR: unable to read the chromosome identifier from the VCF or genetic map." >&2
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

    # Hap-IBD requires every genotype to be phased and non-missing.
    #
    # AWK scans the complete input instead of exiting immediately.
    # This prevents a SIGPIPE error from bcftools under pipefail.
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
        echo "Hap-IBD requires all genotypes to be phased and non-missing." >&2
        exit 1
    fi

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

    python3 ${summary_script} \
        --input ${ibd_file} \
        --output chr${chromosome}.per_pair.tsv.gz \
        --threshold-cm ${params.summary.segment_threshold_cm}
    """
}


/*
 * Calculate genome-wide pair-level IBD summary statistics across
 * all chromosome-specific Hap-IBD segment files.
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
    path ibd_files
    path summary_script

    output:
    path 'genomewide.per_pair.tsv.gz'

    script:
    """
    set -euo pipefail

    python3 ${summary_script} \
        --input ${ibd_files.join(' ')} \
        --output genomewide.per_pair.tsv.gz \
        --threshold-cm ${params.summary.segment_threshold_cm}
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
    path "${method}.membership.tsv.gz"
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


workflow {
    /*
     * Validate configuration parameters.
     */
    if (!params.input_pattern) {
        error 'Set params.input_pattern.'
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

            def geneticMap = file(
                "${geneticMapDirectory}/${geneticMapPrefix}${chromosomeString}.GRCh38.map",
                checkIfExists: true
            )

            tuple(
                chromosome,
                vcf,
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
     * Detect IBD segments using the official GRCh38 PLINK maps.
     */
    hapIBD = HAP_IBD(
        qc.vcfs,
        hapibdJar
    )


    /*
     * Generate chromosome-specific and genome-wide pair summaries.
     */
    PER_PAIR_CHROMOSOME(
        hapIBD.segments,
        summaryScript
    )

    genomewideSummary = PER_PAIR_GENOMEWIDE(
        hapIBD.segments
            .map { chromosome, segmentFile ->
                segmentFile
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

        CLUSTER_GRAPH(
            Channel.fromList(clusteringMethods),
            genomewideSummary,
            validatedSamples,
            clusteringScript
        )
    }
}
