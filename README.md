# IBD-segment
Identical By Descent (IBD) segment detection



INTRODUCTION
------------

CERC-Genomic-Medecine (https://github.com/CERC-Genomic-Medicine)

This pipeline identifies IBD segments in genetic data (sequencing or genotyping).

Implements phasing if necessary

Run hapIBD

Run phaseibd



REQUIREMENTS
------------

Nextflow (tested with version 21)

Python 3.9.6

Required python packages are listed in: requirements.txt 
To install them, use : ```pip install -r requirements.txt```

List of packages used:
-pandas
-numpy
-re



RUNNING
------------

    Clone this repository to the directory where you will run the pipeline:

    git clone https://github.com/

    Modify nextflow.config configuration file.
        params.vcfs -- path to your VCF/BCF file(s). You can use glob expressions to selecect multiple files.
        params.assembly -- set to "GRCh37" or "GRCh38".
        params.vep_cache -- full path to your local vep_cache directory.
        params.vep_flags -- flags you want to pass to VEP.
        

    Run pipeline:

    module load nextflow
    module load singularity
    nextflow run Annotation.nf -w ~/scratch/work_directory

    Important: when working on Compute Canada HPC, set working directory to ~/scratch/<new directory name>. This will speed up IO and also save space on your project partition. After the execution, if there were no errors and you are happy with the results, you can remove this working directory.



OPTIONS
-----------

List of possible options;

-f / --vcf : Input VCF file

-g / --gtf :  Input Gencode GTF file

-i / --genes-in  :  Generate the list of overlapping genes

-a / --genes-around :  Generate the list of genes within +/-200000 base pairs

-n / --genes-nearest : Get the nearest gene 

-o / --output : Prefix of the outputs files. If not specified, output will be redirect as the input VCF file

-h / --help  : Show help message and exit


OUTPUTS	
-----------


InputFile.vcf.gz : Input VCF file with its INFO field modified with additional information required by the tool (GENES_IN, GENES_200KB, GENE_NEAREST)



EXAMPLES
-----------

> python3 CERC_vcf_annotation.py -f <VCF_filename>.vcf.gz -g <GENCODE_filename>.gtf.gz -o "output"




AUTHOR
-----------
PELLETIER Justin (https://www.genomic-medicine-cerc.online/current-team)

email: justin.pelletier2@mcgill.ca
