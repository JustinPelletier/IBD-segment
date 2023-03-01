/*
* AUTHOR: Justin Pelletier, MSc <justin.pelletier2@mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/




process geneticMap {
   errorStrategy "finish"
   beforeScript "module load plink/1.9b_6.21-x86_64"


   input:
   tuple val(chromosome), path(genotype), path(map), path(hapibd),path(gap), val(removegap)

   output:
   tuple val(chromosome), path(genotype), path("*.custom.map"), path(hapibd), path(gap), val(removegap)

   """
   #with a genetic map need to be the exact same variants than the input vcf file (interpolation)
   plink --vcf ${genotype} --cm-map ${map} ${chromosome} --make-bed --out ${genotype.getBaseName()}.custom
   awk '{print \$1" . "\$3" "\$4}' ${genotype.getBaseName()}.custom.bim > ${genotype.getBaseName()}.custom.map

   """
}



process HapIBD {
   errorStrategy "finish"
   publishDir 'Results', pattern: '*.hapibd.header.ibd.gz', mode: "copy"
   beforeScript 'module load bcftools'

   input:
   tuple val(chromosome), path(genotype), path(map), path(hapibd), path(gap), val(removegap)

   output:
   tuple val(chromosome), path("*.hapibd.header.ibd.gz")

   script:
   """
   java -Xmx100g -jar ${hapibd} gt=${genotype} out=${genotype.getBaseName()} map=${map} nthreads=$task.cpus
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
           bgzip -f ${genotype.getBaseName()}.nogap.hapibd.header.ibd

   else
        #add header to the ibd output file
        cat header.tmp ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.hapibd.header.ibd
        bgzip -f ${genotype.getBaseName()}.hapibd.header.ibd

   fi
   """
}



process PhaseIBD {
   cpus 1
   time = "1h"
   memory = "5GB"

   errorStrategy "finish"
   publishDir 'Results', pattern: '*.phaseibd.header.ibd.gz', mode: "copy"
   beforeScript "source ${params.virtualenv} ; module load plink/1.9b_6.21-x86_64 ; module load bcftools"


   input:
   tuple val(chromosome), path(genotype), path(map), path(phaseibd), path(gap), val(removegap)

   output:
   tuple val(chromosome), path("chr*.phaseibd.header.ibd.gz")


   """
   #with a genetic map need to be the exact same variants than the input vcf file (interpolation)
   #plink --vcf ${genotype} --cm-map ${map} ${chromosome} --make-bed --out ${genotype.getBaseName()}.custom
   #awk '{print \$1" . "\$3" "\$4}' ${genotype.getBaseName()}.custom.bim > ${genotype.getBaseName()}.custom.map

   echo "index,VCF_ID" | sed 's/,/\t/g' > ${genotype}.id
   bcftools query -l ${genotype} | awk '{print int((NR-1)) " " \$0}' | sed 's/ /\t/g' >> ${genotype}.id

   #Unzip the vcf for phaseIBD to run
   gunzip -c ${genotype} > ${genotype.getBaseName()}

   python3 ${phaseibd} ${genotype.getBaseName()} ${chromosome} ${genotype.getBaseName()} ${map} ${genotype}.id


   #remove IBD overlapping gaps in the genome
   if [ ${removegap} == "true" ]
   then
           head -1 ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.ibd.header
           sed -i '1d' ${genotype.getBaseName()}.ibd

           FILENAME=${gap}
           IFS=\$'\t'
           while read CHOM START END TYPE; do
                awk -v a=\$START -v b=\$END '{ if (!( (\$6<=a && \$7>=a) || (\$6<=b && \$7>=b) || (\$6>=a && \$7<=b) )) { print }}' ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.ibd.tmp
                mv ${genotype.getBaseName()}.ibd.tmp ${genotype.getBaseName()}.ibd
           done < \$FILENAME

           cat ${genotype.getBaseName()}.ibd.header ${genotype.getBaseName()}.ibd > ${genotype.getBaseName()}.nogap.phaseibd.ibd
           bgzip -f ${genotype.getBaseName()}.nogap.phaseibd.ibd

   else
        mv ${genotype.getBaseName()}.ibd ${genotype.getBaseName()}.phaseibd.ibd
        bgzip -f ${genotype.getBaseName()}.phaseibd.ibd


   fi

   """
}



workflow {

   //genetic_map = Channel.from(params.chromosomes).map { chr -> [ "${chr}" , params.genoFile + ".chr${chr}.vcf.gz",  params.geneticMap + ".chr${chr}." + params.assembly +".gmap", params.phaseibd,  params.gapfile, params.removeGaps] } | geneticMap | HapIBD

   genetic_map = Channel.from(params.chromosomes).map { chr -> [ "${chr}" , params.genoFile + ".chr${chr}.vcf.gz",  params.geneticMap + ".chr${chr}." + params.assembly +".gmap", params.phaseibd,  params.gapfile, params.removeGaps] } | geneticMap | PhaseIBD
   genetic_map.view()



   //paf = genetic_map.flatten()
   //paf = genetic_map.groupTuple().map { file -> [ ${file}, params.hapibd ] }
   //paf.view()

   //hapibd_out = genetic_map.groupTuple().map { file -> [ $file, params.hapibd ] } | HapIBD
   //hapibd_out = Channel.from(params.chromosomes).map { chr -> [ "chr${chr}" , params.genoFile + ".chr${chr}.vcf.gz", params.geneticMap + ".chr${chr}." + params.assembly +".map", params.hapibd, params.gapfile, params.removeGaps] } | HapIBD
   //hapibd_out.view()

   //phaseibd_out = Channel.from(params.chromosomes).map { chr -> [ "${chr}" , params.genoFile + ".chr${chr}.vcf.gz", genetic_map , params.phaseibd, params.gapfile, params.removeGaps] } | PhaseIBD
   //phaseibd_out = Channel.from(params.chromosomes).map { chr -> [ "${chr}" , params.genoFile + ".chr${chr}.vcf.gz", params.geneticMap + ".chr${chr}." + params.assembly +".map", params.phaseibd, params.gapfile, params.removeGaps] } | PhaseIBD
   //phaseibd_out.view()



}

