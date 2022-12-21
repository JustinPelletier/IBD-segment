!/usr/bin/env nextflow

/*
* AUTHOR: Justin Pelletier, MSc <justin.pelletier2@mcgill.ca>
* VERSION: 1.0
* YEAR: 2022
*/


process Phase {  
   cpus 8
   memory "16 GB"
   time "4h"
   errorStrategy "finish"
   cache "lenient"
   scratch true

   
   input:
   path(genotype)
   path(map)
   path(beagle)
   each chromosome

   output:
   path("*.vcf.gz")
   
   publishDir "BeaglePhased/", pattern: "*.vcf.gz", mode: "copy"
   
   
   """
   java -jar ${beagle} gt=${genotype} out=${genotype}.phased map=${map} nthreads=8 
   """
     
}


//process MergeCHR { 
   cpus 1
   memory "20 GB"
   time "1h"
   errorStrategy "finish"
   cache "lenient"
   scratch true

   
   input:
   path(genotype)
   val(chromosome)
   output:
   path "${genotype}.merge.vcf.gz"

   """
   for chr in $(seq ${chromosome})
   do
   	if [ $chr == 1 ]
	then
   		zgrep "^#" ${genotype} > ${genotype}.merge.vcf
	fi
	zgrep -v "^#" ${genotype} >> ${genotype}.merge.vcf
   done
   bgzip -f ${genotype}.merge.vcf
   """
}


process HapIBD {  
   cpus 8
   memory "16 GB"
   time "4h"
   errorStrategy "finish"
   cache "lenient"
   scratch true

   
   input:
   path(hapibd)
   path(genotype)
   path(map)
   path(beagle)
   path(gap)
   each chromosome

   output:
   path("${genotype}.nogap.header.ibd")
   
   script:
   if (param.phased == FLASE) {
      """
      #launch Phase process
      """
   }
   //then launch the hap-ibd program
   """   
   java -Xmx100g -jar ${hapibd} gt=${genotype} out=${genotype} map=${map} nthreads=8 
   gunzip ${genotype}.ibd.gz
   
   #remove IBD overlapping gaps in the genome
   FILENAME=${gap}
   IFS=$'\t'
   
   while read CHOM START END TYPE; do
	awk -v a=$START -v b=$END '{ if (!( ($6<=a && $7>=a) || ($6<=b && $7>=b) || ($6>=a && $7<=b) )) { print }}' ${genotype}.ibd > ${genotype}.ibd.tmp
	mv ${genotype}.ibd.tmp ${genotype}.ibd
   done <$FILENAME
   
   #add header to the ibd output file
   echo "SAMPLE1 HAP_INDEX1      SAMPLE2 HAP_INDEX2      CHROM   START   END     GEN_LENGTH"  > header.tmp
   cat header.tmp ${genotype}.ibd > ${genotype}.nogap.ibd
   bgzip -f ${genotype}.nogap.header.ibd
   """   
}



process PhaseIBD {  
   cpus 8
   memory "16 GB"
   time "4h"
   errorStrategy "finish"
   cache "lenient"
   scratch true

   
   input:
   path(phaseibd)
   path(genotype)
   path(map)
   path(beagle)
   each chromosome

   output:
   path("${chromosome}.nogap.header.ibd.gz")
   
   script:
   if (param.phased == FLASE) {
      """
      #launch Phase process
      """
   }
   //then launch the hap-ibd program
   """   
   #with a genetic map need to be the exact same variants than the input vcf file (interpolation)
   plink --vcf ${genotype} --cm-map ${map} ${chromsome} --make-bed --out ${TMPDIR}.${chromosome}.custom
   awk '{print $1" . "$3" "$4}' ${TMPDIR}.${chromosome}.custom.bim > ${map}.${chromosome}.custom.map
   
   echo "index,VCF_ID" | sed 's/,/\t/g' > ${genotype}.id
   bcftools query -l ${genotype} | awk '{print int((NR-1)) " " $0}' | sed 's/ /\t/g' >> ${genotype}.id
   
   python3 ${phaseibd} ${genotype} ${chromosome} ${TMPDIR} ${map}.${chromosome}.custom.map ${genotype}.id
   head -1 ${chromosome}.ibd > ${chromosome}.ibd.header
   sed -i '1d' ${chromosome}.ibd
   
   #remove IBD overlapping gaps in the genome
   FILENAME=${gap}
   IFS=$'\t'
   while read CHOM START END TYPE; do
   	awk -v a=$START -v b=$END '{ if (!( ($10<=a && $11>=a) || ($10<=b && $11>=b) || ($10>=a && $11<=b) )) { print }}' ${chromosome}.ibd > ${chromosome}.ibd.tmp
	mv ${chromosome}.ibd.tmp ${chromosome}.ibd
   done <$FILENAME
   
   #add header to the ibd output file
   cat ${chromosome}.ibd.header ${chromosome}.ibd > ${chromosome}.nogap.header.ibd
   bgzip -f ${chromosome}.nogap.header.ibd
   """   
}




process RemoveGaps {  
   cpus 1
   memory "8 GB"
   time "4h"
   scratch true
   
   input:
   path(phaseibd_result)
   path(hapibd_result)
   path(gap)
   path(beagle)
   each chromosome

   output:
   path("*.")
   
   script:
   if (param.gaps == TRUE) {
   """
   #remove IBD overlapping gaps in the genome
   FILENAME=${gap}
   IFS=$'\t'
   
   while read CHOM START END TYPE; do
	   #echo "CHOM=$CHOM START=$START END=$END TYPE=$TYPE"

	   awk -v a=$START -v b=$END '{ if (!( ($6<=a && $7>=a) || ($6<=b && $7>=b) || ($6>=a && $7<=b) )) { print }}' $path_gap/cp_nogap_hapIBD_cag_chrCHR.ibd > $path_gap/tmp_nogap_hapIBD_cag_chrCHR.ibd
	   mv $path_gap/tmp_nogap_hapIBD_cag_chrCHR.ibd $path_gap/cp_nogap_hapIBD_cag_chrCHR.ibd
	   OUTPUT=$(wc -l $path_gap/cp_nogap_hapIBD_cag_chrCHR.ibd | cut -d " " -f1)

	   awk -v a=$START -v b=$END '{ if (!( ($10<=a && $11>=a) || ($10<=b && $11>=b) || ($10>=a && $11<=b) )) { print }}' $path_gap/cp_nogap_phaseIBD_cag_chrCHR.ibd > $path_phaseIBD/tmp_nogap_phaseIBD_cag_chrCHR.ibd
	   mv $path_phaseIBD/tmp_nogap_phaseIBD_cag_chrCHR.ibd $path_phaseIBD/cp_nogap_phaseIBD_cag_chrCHR.ibd
	   TEST=$(wc -l $path_phaseIBD/cp_nogap_phaseIBD_cag_chrCHR.ibd | cut -d " " -f1)

   done <$FILENAME
   
   """
   }else{
   """
   echo "Gaps not removed"
   """   
   }   
}
	
	
	
workflow {
   
   chromosomes = Channel.from(params.chromosomes)
   
   phased_geno = Phase()
}
