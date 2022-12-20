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
   each chromosome

   output:
   path("*.header.ibd.gz")
   
   script:
   if (param.phased == FLASE) {
      """
      #launch Phase process
      """
   }
   //then launch the hap-ibd program
   """   
   java -Xmx100g -jar ${hapibd} gt=${genotype} out=${genotype}.phased map=${map} nthreads=8 
   
   #add header to the ibd output file
   echo "SAMPLE1 HAP_INDEX1      SAMPLE2 HAP_INDEX2      CHROM   START   END     GEN_LENGTH" | bgzip -f > header.tmp.gz
   cat header.tmp.gz ${genotype}.phased.ibd.gz > {genotype}.phased.header.ibd.gz 
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
   path("*.")
   
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
   
   """   
   
   
}


process RemoveGaps {  
   cpus 1
   memory "8 GB"
   time "4h"
   scratch true
   
   
   
   script:
   if (param.gaps == TRUE) {
   """
   
   """
   }else{
   """
   
   
   """   
   }
   
   
   
}
	
	
	
workflow {
   
   chromosomes = Channel.from(params.chromosomes)
   
   phased_geno = Phase(path(genotype), path(map), path(beagle), chromosomes)
}
