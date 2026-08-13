/*
 * Hap-IBD segment-detection pipeline
 *
 * Author: Justin Pelletier
 * Version: 2.0
 */

nextflow.enable.dsl = 2


/*
 * Filter variants according to MAF and variant missingness.
 * The process also confirms that all retained genotypes are phased
 * and contain no missing alleles, as required by Hap-IBD.
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
    tuple val(chromosome), path(vcf), path(map), path(gap)

    output:
    tuple val(chromosome),
          path("chr${chromosome}.qc.vcf.gz"),
          path("chr${chromosome}.qc.vcf.gz.tbi"),
          path(map),
          path(gap),
          emit: vcfs

    tuple val(chromosome),
          path("chr${chromosome}.samples.txt"),
          emit: samples

    path "chr${chromosome}.qc.stats.txt",
          emit: stats

    script:
    """
    set -o pipefail

    bcftools +fill-tags ${vcf} -Ou -- -t MAF,F_MISSING | \
        bcftools view \
            -i 'MAF>=${params.qc.min_maf} && F_MISSING<=${params.qc.max_variant_missingness}' \
            -Oz \
            -o chr${chromosome}.qc.vcf.gz

    tabix -f -p vcf chr${chromosome}.qc.vcf.gz

    bcftools query \
        -l chr${chromosome}.qc.vcf.gz \
        > chr${chromosome}.samples.txt

    if [[ ! -s chr${chromosome}.samples.txt ]]; then
        echo "ERROR: no samples were found in the chromosome ${chromosome} VCF." >&2
        exit 1
    fi

    if [[ \$(bcftools index -n chr${chromosome}.qc.vcf.gz) -eq 0 ]]; then
        echo "ERROR: no variants remained after QC on chromosome ${chromosome}." >&2
        exit 1
    fi

    # Hap-IBD requires every genotype to be phased and non-missing.
    if bcftools query -f '[%GT\\n]' chr${chromosome}.qc.vcf.gz | \
       awk '
           index(\$0, "/") || index(\$0, ".") {
               invalid = 1
               exit
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
 * Confirm that the sample IDs and their ordering are identical across
 * all chromosome-specific VCF files.
 *
 * The resulting cohort.samples.txt file is also used to include isolated
 * participants in the clustering graph.
 */
process VALIDATE_SAMPLES {
    tag 'validate sample lists'

    cpus 1
    memory '1 GB'
    time '30m'

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
    cp ${sample_files[0]} cohort.samples.txt

    for sample_file in ${sample_files.join(' ')}
    do
        if ! cmp -s cohort.samples.txt \${sample_file}
        then
            echo "ERROR: sample IDs or sample ordering differ between chromosomes." >&2
            exit 1
        fi
    done
    """
}


/*
 * Interpolate genetic positions for the variants retained after QC.
 *
 * The generated map contains exactly the markers present in the
 * chromosome-specific QC-filtered VCF.
 */
process GENETIC_MAP {
    tag "chr${chromosome}"

    cpus params.resources.map.cpus
    memory params.resources.map.memory
    time params.resources.map.time

    beforeScript """
    module load StdEnv/2020
    module load plink/1.9b_6.21-x86_64
    """

    errorStrategy 'terminate'

    publishDir "${params.outdir}/maps",
        pattern: '*.custom.map',
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chromosome),
          path(vcf),
          path(index),
          path(map),
          path(gap)

    output:
    tuple val(chromosome),
          path(vcf),
          path(index),
          path("chr${chromosome}.custom.map"),
          path(gap)

    script:
    """
    plink \
        --vcf ${vcf} \
        --cm-map ${map} ${chromosome} \
        --double-id \
        --make-bed \
        --out chr${chromosome}.custom

    awk '
        BEGIN {
            OFS = "\\t"
        }
        {
            print \$1, ".", \$3, \$4
        }
    ' chr${chromosome}.custom.bim \
        > chr${chromosome}.custom.map
    """
}


/*
 * Detect IBD segments using Hap-IBD.
 *
 * If gap removal is enabled, any complete IBD segment overlapping
 * at least one gap interval is removed.
 */
process HAP_IBD {
    tag "chr${chromosome}"

    cpus params.resources.hapibd.cpus
    memory params.resources.hapibd.memory
    time params.resources.hapibd.time

    beforeScript """
    module load java/25.36
    module load bcftools
    module load bedtools
    """

    errorStrategy 'retry'
    maxRetries 1

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
          path(index),
          path(map),
          path(gap)

    path hapibd

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
              bcftools query \
                  -f '%CHROM\\n' \
                  ${vcf} |
              head -n 1
          )
    
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
              ' ${gap} \
              > chr${chromosome}.normalized.gaps.bed
    
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
          """ \
        : """
          cp \
              chr${chromosome}.raw.ibd \
              chr${chromosome}.filtered.ibd
          """
    """
    set -o pipefail

    java \
        -Xmx${javaGb}g \
        -jar ${hapibd} \
        gt=${vcf} \
        map=${map} \
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
 * Generate a pair-level IBD summary for each chromosome.
 */
process PER_PAIR_CHROMOSOME {
    tag "chr${chromosome}"

    cpus params.resources.summary.cpus
    memory params.resources.summary.memory
    time params.resources.summary.time

    beforeScript """
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/per_pair/by_chromosome",
        pattern: '*.per_pair.tsv.gz',
        mode: 'copy',
        overwrite: true

    input:
    tuple val(chromosome), path(ibd_file)
    path summary_script

    output:
    tuple val(chromosome),
          path("chr${chromosome}.per_pair.tsv.gz")

    script:
    """
    python3 ${summary_script} \
        --input ${ibd_file} \
        --output chr${chromosome}.per_pair.tsv.gz \
        --threshold-cm ${params.summary.segment_threshold_cm}
    """
}


/*
 * Combine the chromosome-specific Hap-IBD results and calculate
 * genome-wide pair-level summary statistics.
 */
process PER_PAIR_GENOMEWIDE {
    tag 'genome-wide summary'

    cpus params.resources.summary.cpus
    memory params.resources.summary.memory
    time params.resources.summary.time

    beforeScript """
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
    tag "${method} recursive clustering"

    cpus params.resources.clustering.cpus
    memory params.resources.clustering.memory
    time params.resources.clustering.time

    beforeScript """
    source ${params.python_venv}/bin/activate
    """

    publishDir "${params.outdir}/clustering/${method}",
        mode: 'copy',
        overwrite: true

    input:
    tuple val(method), val(max_levels)
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
     * Validate required configuration parameters.
     */
    if (!params.input_pattern) {
        error 'Set params.input_pattern.'
    }

    if (!params.genetic_map_pattern) {
        error 'Set params.genetic_map_pattern.'
    }

    if (params.remove_gaps && !params.gap_file) {
        error 'Set params.gap_file when remove_gaps=true.'
    }

    if (!(params.chromosomes as List)) {
        error 'params.chromosomes must contain at least one chromosome.'
    }

    if (params.qc.min_maf < 0 || params.qc.min_maf > 0.5) {
        error 'params.qc.min_maf must be between 0 and 0.5.'
    }

    if (
        params.qc.max_variant_missingness < 0 ||
        params.qc.max_variant_missingness > 1
    ) {
        error 'params.qc.max_variant_missingness must be between 0 and 1.'
    }

    if (params.n_louvain < 0 || params.n_leiden < 0) {
        error 'Clustering refinement depths must be greater than or equal to zero.'
    }


    /*
     * Construct one input tuple per chromosome.
     */
    inputs = Channel
        .fromList(params.chromosomes as List)
        .map { chromosome ->
            def chromosomeString = chromosome as String

            def vcf = file(
                params.input_pattern.replace('{chr}', chromosomeString),
                checkIfExists: true
            )

            def map = file(
                params.genetic_map_pattern.replace('{chr}', chromosomeString),
                checkIfExists: true
            )

            def gap = file(
                params.gap_file,
                checkIfExists: true
            )

            tuple(
                chromosome,
                vcf,
                map,
                gap
            )
        }


    /*
     * Perform variant QC and validate the sample lists.
     */
    qc = QC_VCF(inputs)

    validatedSamples = VALIDATE_SAMPLES(
        qc.samples
            .map { chromosome, samples -> samples }
            .collect()
    )


    /*
     * Generate marker-matched maps when requested.
     */
    if (params.generate_custom_maps) {
        customMaps = GENETIC_MAP(qc.vcfs)
        mapsForHapIBD = customMaps.out
    } else {
        mapsForHapIBD = qc.vcfs
    }


    /*
     * Detect IBD segments.
     */
    hapIBD = HAP_IBD(
        mapsForHapIBD,
        file(
            params.hapibd_jar,
            checkIfExists: true
        )
    )


    /*
     * Generate chromosome-specific and genome-wide summaries.
     */
    summaryScript = file(
        params.per_pair_script,
        checkIfExists: true
    )

    PER_PAIR_CHROMOSOME(
        hapIBD.segments,
        summaryScript
    )

    genomewideSummary = PER_PAIR_GENOMEWIDE(
        hapIBD.segments
            .map { chromosome, segmentFile -> segmentFile }
            .collect(),
        summaryScript
    )


    /*
     * Build the list of requested clustering analyses.
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
     * first() converts the single genome-wide summary and sample-list
     * outputs into reusable value channels. This allows both Louvain and
     * Leiden jobs to consume the same files.
     */
    if (clusteringMethods) {
        CLUSTER_GRAPH(
            Channel.fromList(clusteringMethods),
            genomewideSummary.out.first(),
            validatedSamples.out.first(),
            file(
                params.clustering_script,
                checkIfExists: true
            )
        )
    }
}
