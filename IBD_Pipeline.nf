/*
* AUTHOR: Justin Pelletier, MSc <justin.pelletier2@mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/

process geneticMap {
   time = "1h"
   memory = "5GB"
   cpus 1

   errorStrategy "finish"
   beforeScript "module load plink/1.9b_6.21-x86_64"


   input:
   tuple val(chromosome), path(genotype), path(map),path(gap), val(removegap)

   output:
   tuple val(chromosome), path(genotype), path("*.custom.map"), path(gap), val(removegap)

   """
   #with a genetic map need to be the exact same variants than the input vcf file (interpolation)
   plink --vcf ${genotype} --cm-map ${map} ${chromosome} --double-id --make-bed --out ${genotype.getBaseName()}.custom
   awk '{print \$1" . "\$3" "\$4}' ${genotype.getBaseName()}.custom.bim > ${genotype.getBaseName()}.custom.map

   """
}




process HapIBD {
   errorStrategy "retry"
   maxRetries 1
   publishDir 'test_Results_200bp', pattern: '*.hapibd.header.ibd.gz', mode: "copy"
   beforeScript 'module load bcftools'

   input:
   tuple val(chromosome), path(genotype), path(map), path(gap), val(removegap)
   path(hapibd)
   val(min_size)
   val(min_markers)

   output:
   tuple (val(chromosome), path("${genotype.getBaseName()}.ibd"), path("header.tmp"), emit: hapfiles)
   path("*.hapibd.header.ibd.gz")
   //tuple val(chromosome), path("*.hapibd.header.ibd.gz")


   script:
   """
   java -Xmx100g -jar ${hapibd} gt=${genotype} min-output=${min_size} min-markers=${min_markers} out=${genotype.getBaseName()} map=${map} nthreads=$task.cpus
   gunzip ${genotype.getBaseName()}.ibd.gz

   #add header to the output file
   echo "SAMPLE1 HAP_INDEX1      SAMPLE2 HAP_INDEX2      CHROM   START   END     GEN_LENGTH"  > header.tmp

   #remove IBD overlapping gaps in the genome
   if [ ${removegap} == "true" ]
   then
           FILENAME=${gap}
           IFS=\$'\t'
           while read CHOM START END TYPE; do
                awk -v a=\$START -v b=\$END '{ if (!( (\$6<=a && \$7>=a) || (\$6<=b && \$7>=b) || (\$6>=a && \$7<=b) )) { print }}' ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.ibd.tmp
                mv ${genotype.getBaseName()}.ibd.tmp ${genotype.getBaseName()}.ibd
           done < \$FILENAME

           #add header to the ibd output file
           cat header.tmp ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.nogap.hapibd.header.ibd
           bgzip -c ${genotype.getBaseName()}.nogap.hapibd.header.ibd > ${genotype.getBaseName()}.nogap.hapibd.header.ibd.gz

   else
           #add header to the ibd output file
           cat header.tmp ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.hapibd.header.ibd
           bgzip -c ${genotype.getBaseName()}.hapibd.header.ibd > ${genotype.getBaseName()}.hapibd.header.ibd.gz

   fi
   """
}



process PhaseIBD {
   cpus 1

   errorStrategy "retry"
   maxRetries 1
   publishDir 'test_Results_200bp', pattern: '*.phaseibd.header.ibd.gz', mode: "copy"
   beforeScript "source ${params.virtualenv} ; module load bcftools"


   input:
   tuple val(chromosome), path(genotype), path(map), path(gap), val(removegap)
   path(phaseibd)
   val(min_size)
   val(min_markers)

   output:
   tuple (val(chromosome), path("${genotype.getBaseName()}.ibd"), path("${genotype.getBaseName()}.ibd.header"), emit: phasefiles)
   path("*.phaseibd.header.ibd.gz")


   """

   echo "index,VCF_ID" | sed 's/,/\t/g' > ${genotype}.id
   bcftools query -l ${genotype} | awk '{print int((NR-1)) " " \$0}' | sed 's/ /\t/g' >> ${genotype}.id

   #Unzip the vcf for phaseIBD to run
   gunzip -c ${genotype} > ${genotype.getBaseName()}

   python3 ${phaseibd} ${genotype.getBaseName()} ${chromosome} ${genotype.getBaseName()} ${map} ${genotype}.id ${min_size} ${min_markers}


   #remove IBD overlapping gaps in the genome
   if [ ${removegap} == "true" ]
   then
           head -1 ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.ibd.header
           sed -i '1d' ${genotype.getBaseName()}.ibd

           FILENAME=${gap}
           IFS=\$'\t'
           while read CHOM START END TYPE; do
                awk -v a=\$START -v b=\$END '{ if (!( (\$10<=a && \$11>=a) || (\$10<=b && \$11>=b) || (\$10>=a && \$11<=b) )) { print }}' ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.ibd.tmp
                mv ${genotype.getBaseName()}.ibd.tmp ${genotype.getBaseName()}.ibd
           done < \$FILENAME

           cat ${genotype.getBaseName()}.ibd.header ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.nogap.phaseibd.header.ibd
          bgzip -c ${genotype.getBaseName()}.nogap.phaseibd.header.ibd > ${genotype.getBaseName()}.nogap.phaseibd.header.ibd.gz

   else
           mv ${genotype.getBaseName()}.ibd ${genotype.getBaseName()}.phaseibd.header.ibd
           bgzip -c ${genotype.getBaseName()}.phaseibd.header.ibd > ${genotype.getBaseName()}.phaseibd.header.ibd.gz
   fi
   """
}



process PerPairHapIBD {
   cpus 1

   publishDir 'test_Results_200bp', pattern: '*.per_pair.ibd', mode: "copy"
   beforeScript "source ${params.virtualenv} ; module load bcftools"


   input:
   tuple val(chromosome), path(ibd_file), path(header)
   val(method)
   path(script_per_pair)

   output:
   tuple val(chromosome), path("${method}.${chromosome}.per_pair.ibd"), path(header)


   """
   python3 ${script_per_pair} ${ibd_file} ${method}.${chromosome}.per_pair.ibd $chromosome $method
   #bgzip -c ${method}.${chromosome}.per_pair.ibd > ${method}.${chromosome}.per_pair.ibd.gz
   """
}


process PerPairPhaseIBD {
   cpus 1

   publishDir 'test_Results_200bp', pattern: '*.per_pair.ibd', mode: "copy"
   beforeScript "source ${params.virtualenv} ; module load bcftools"


   input:
   tuple val(chromosome), path(ibd_file), path(header)
   val(method)
   path(script_per_pair)

   output:
   tuple val(chromosome), path("${method}.${chromosome}.per_pair.ibd"), path(header)


   """
   python3 ${script_per_pair} ${ibd_file} ${method}.${chromosome}.per_pair.ibd $chromosome $method
   #bgzip -c ${method}.${chromosome}.per_pair.ibd > ${method}.${chromosome}.per_pair.ibd.gz
   """
}


process MergeIBD {
   cpus 1

   publishDir 'test_Results_200bp', pattern: '*.merged.ibd.gz', mode: "copy"
   beforeScript "module load bcftools"


   input:
   tuple val(chromosome), path(ibd_files), path(header)

   output:
   tuple val(chromosome), path("$ibd_files.getSimpleName().merged.ibd.gz")


   """
   cat $header > $ibd_files.getSimpleName().merged.ibd
   for f in ${ibd_files}; do cat \${f} >> $ibd_files.getSimpleName().merged.ibd ; done
   bgzip -c $ibd_files.getSimpleName().merged.ibd > $ibd_files.getSimpleName().merged.ibd.gz
   """
}






workflow {

   geneticMaps_out = Channel.from(params.chromosomes).map { chr -> [ "${chr}" , params.genoFile + ".chr${chr}.vcf.gz",  params.geneticMap + ".chr${chr}." + params.assembly +".gmap",  params.gapfile + params.assembly + ".chr${chr}.gap.bed", params.removeGaps] } | geneticMap

   hapIBD_out =  HapIBD(geneticMaps_out, params.hapibd, params.minimun_size, params.minimum_markers)

   phaseIBD_out = PhaseIBD(geneticMaps_out, params.phaseibd, params.minimun_size, params.minimum_markers)

   IBD_pair_hapIBD = PerPairHapIBD(hapIBD_out.hapfiles, "HapIBD", params.per_pair)
   IBD_pair_phaseIBD = PerPairPhaseIBD(phaseIBD_out.phasefiles, "PhaseIBD", params.per_pair)

   //out_merge_hapIBD = MergeIBD(IBD_pair_HapIBD.groupTuple())
   //out_merge_phaseIBD = MergeIBD(IBD_pair_PhaseIBD.groupTuple(by: [0, 1]))


}

