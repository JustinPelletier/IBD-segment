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
   val phasestep
   tuple path(genotype), path(map), path(beagle), path(out)
   each chromosome

   output:
   tuple val(chromosome), path("${out}/${genotype.getSimpleName()}.chr${chromosome}.phased.vcf.gz"), path(beagle), path(map), path(out)
      
   script:
   if( CHR == 1 && phasestep == 1)
   	"""
   	java -jar ${beagle} gt=${genotype} out=${out}/${genotype.getSimpleName()}.chr${chromosome}.phased map=${map}/plink.${chromosome}.GRCh38.map nthreads=8 
   	"""
   if(CHR == 0 && phasestep == 1)
    	"""
	java -jar ${beagle} gt=${genotype} out=${out}/${genotype.getSimpleName()}.chr${chromosome}.phased map=${map}/no_chr_plink.${chromosome}.GRCh38.map nthreads=8
	"""
   else //don't phase, just output the tuple while changing the input vcf name
   	"""
	cp ${genotype} ${out}/${genotype.getSimpleName()}.chr${chromosome}.phased.vcf.gz
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
   tuple val(chromosome), path(genotype), path(beagle), path(map), path(out)
   path(hapibd)
   path(gap)
   val removegap
   
   output:
   path("$$out/chr*.header.ibd.gz")
   
   script:
   if (param.phased == FLASE) {
      """
      #launch Phase process
      """
   }
   //then launch the hap-ibd program
   """   
   java -Xmx100g -jar ${hapibd} gt=${genotype} out=$out/chr${chromosome} map=${map} nthreads=8 
   gunzip $out/chr${chromosome}.ibd.gz
   
   echo "SAMPLE1 HAP_INDEX1      SAMPLE2 HAP_INDEX2      CHROM   START   END     GEN_LENGTH"  > header.tmp
   
   #remove IBD overlapping gaps in the genome and AccessibilityMask
   if [ removegap == 1 ]
   then
	   FILENAME=${gap}
	   IFS=$'\t'
	   while read CHOM START END TYPE; do
		awk -v a=$START -v b=$END '{ if (!( ($6<=a && $7>=a) || ($6<=b && $7>=b) || ($6>=a && $7<=b) )) { print }}' $out/chr${chromosome}.ibd > $out/chr${chromosome}.ibd.tmp
		mv $out/chr${chromosome}.ibd.tmp $out/chr${chromosome}.ibd
	   done <$FILENAME
	   
	   #add header to the ibd output file
   	   cat header.tmp $out/chr${chromosome}.ibd > $out/chr${chromosome}.nogap.ibd
   	   bgzip -f $out/chr${chromosome}.nogap.header.ibd
	   
   else
   	#add header to the ibd output file
   	cat header.tmp $out/chr${chromosome}.ibd > $out/chr${chromosome}.ibd
   	bgzip -f $$out/chr${chromosome}.header.ibd
   
   fi
   """   
}



process PhaseIBD {  
   cpus 8
   memory "16 GB"
   time "4h"
   errorStrategy "finish"
   cache "lenient"
   scratch true

   beforeScript "source ${params.virtualenv}"
   
   input:
   tuple val(chromosome), path(genotype), path(beagle), path(map), path(out)
   path(phaseibd)
   path(gap)
   val removegap

   output:
   path("${out}/chr*.header.ibd.gz")
   
 
   """   
   #with a genetic map need to be the exact same variants than the input vcf file (interpolation)
   plink --vcf ${genotype} --cm-map ${map} ${chromsome} --make-bed --out ${TMPDIR}.${chromosome}.custom
   awk '{print $1" . "$3" "$4}' ${TMPDIR}.${chromosome}.custom.bim > ${map}.${chromosome}.custom.map
   
   echo "index,VCF_ID" | sed 's/,/\t/g' > ${genotype}.id
   bcftools query -l ${genotype} | awk '{print int((NR-1)) " " $0}' | sed 's/ /\t/g' >> ${genotype}.id
   
   python3 ${phaseibd} ${genotype} ${chromosome} ${out} ${map}.${chromosome}.custom.map ${genotype}.id
   
   
   #remove IBD overlapping gaps in the genome and AccessibilityMask
   if [ removegap == 1 ]
   then
   	head -1 ${out}/${chromosome}.ibd > ${out}/${chromosome}.ibd.header
   	sed -i '1d' ${out}/${chromosome}.ibd
	
   	FILENAME=${gap}
	IFS=$'\t'
	while read CHOM START END TYPE; do
		awk -v a=$START -v b=$END '{ if (!( ($10<=a && $11>=a) || ($10<=b && $11>=b) || ($10>=a && $11<=b) )) { print }}' ${out}/${chromosome}.ibd > ${out}/${chromosome}.ibd.tmp
		mv ${out}/${chromosome}.ibd.tmp ${out}/${chromosome}.ibd
	done <$FILENAME
	
	#add header to the ibd output file
   	cat ${out}/${chromosome}.ibd.header ${out}/${chromosome}.ibd > ${out}/chr${chromosome}.nogap.header.ibd
   	bgzip -f ${out}/chr${chromosome}.nogap.header.ibd
   else
   	#add header to the ibd output file
	mv {out}/${chromosome}.ibd {out}/chr${chromosome}.header.ibd
   	bgzip -f ${out}/chr${chromosome}.header.ibd
   fi
   """   
}


	
	
workflow {
   chromosomes = Channel.from(params.chromosomes)
   phasingStep = Channel.from(params.phasingStep)
   removegaps = Channel.from(params.removeGaps)
   
   genotypes = Channel.fromPath(params.genoFile)
   workdir = Channel.fromPath(params.workingDir)
   beagle =  Channel.fromPath(params.beagleDir)
   phaseIBD = Channel.fromPath(params.phaseibd)
   hapIBD = Channel.fromPath(params.hapibd)
   geneticmap = Channel.fromPath(params.geneticMap)
   gaps = Channel.fromPath(params.gapfile)
   
   
   phased_geno = Phase(phasingStep, [genotypes, geneticmap, beagle, workdir], chromosomes ) //val(chromosome), path("${out}/${genotype.getSimpleName()}.chr${chromosome}.phased.vcf.gz"), path(beagle), path(map), path(out)
   hapIBD_segments = HapIBD(phased_geno, hapIBD, gaps, removegaps)
   phaseIBD_segments = PhaseIBD(phased_geno, phaseIBD, gaps, removegaps)
}
