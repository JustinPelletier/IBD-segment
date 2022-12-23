# IBD-segment
Identical By Descent (IBD) segment detection



INTRODUCTION
------------

CERC-Genomic-Medecine (https://github.com/CERC-Genomic-Medicine)

This pipeline identifies IBD segments in genetic data (sequencing or genotyping).

Run Beagle phasing if input data is unphased (https://faculty.washington.edu/browning/beagle/beagle.html#download)

Run hap-ibd (https://github.com/browning-lab/hap-ibd)

Run phaseibd (https://github.com/23andMe/phasedibd) 



REQUIREMENTS
------------

Nextflow (tested with version 21)
module load tabix/0.2.6

Python 3.9.6

Installation of hap-ibd (https://github.com/browning-lab/hap-ibd)
> wget https://faculty.washington.edu/browning/hap-ibd.jar


Installation of phaseibd as a python module(https://github.com/23andMe/phasedibd) 
> make
> python setup.py install 
> python tests/unit_tests.py



RUNNING
------------

    Clone this repository to the directory where you will run the pipeline:

    git clone https://github.com/

    Modify nextflow.config configuration file.
        params.workingDir -- path to the working output directory
        params.virtualenv -- path to the virtual environment where to run phaseibd in python
        params.genoFile -- path to your VCF file. You can use glob expressions to selecect multiple files.
        params.beagleDir -- path to beagle executable java (.jar) file
        params.phaseibd -- path to phaseibd python custom script (phaseibd_chr.py) file that is this github repository
        params.hapibd -- path to hap-ibd executable java (.jar) file
        params.geneticMap -- path to the genetic map of your choice (GRCh36, GRCh37 or GRCh38) available in this repository
        params.gapfile -- path to bed file with intervals to remove from the IBD segments
        params.phasingStep -- Specify if your data is phased (0) or if you want to phase it (1)
        params.chromosomePrefix -- Specify if your genotype file (genoFile) has "chr" as chromosome prefix (1) or not (0)
        params.removeGaps -- Specify if you want the gaps to be removed from the IBD segment (1) or kept (0)

        

    Run pipeline:

    module load nextflow
    module load singularity
    nextflow run IBD_Pipeline.nf -w ~/scratch/workingDirectory

    Important: when working on Compute Canada HPC, set working directory to ~/scratch/<new directory name>. This will speed up IO and also save space on your project partition. After the execution, if there were no errors and you are happy with the results, you can remove this working directory.



OPTIONS
-----------

List of possible options;



OUTPUTS	
-----------


chr*.hapibd.header.ibd.gz -- hap-ibd IBD segment (with or without gaps)
chr*.nogap.hapibd.header.ibd.gz

chr*.phaseibd.header.ibd.gz -- phaseibd IBD segment detected (with or without gaps)
chr*.nogap.phaseibd.header.ibd.gz


EXAMPLES
-----------

> 




AUTHOR
-----------
PELLETIER Justin (https://www.genomic-medicine-cerc.online/current-team)

email: justin.pelletier2@mcgill.ca
