# NGSclean  
A pipeline to trim the reads and remove ribosomal RNA from RNA-Seq data  


This pipeline trims the adaptor from the ends of reads, and move plant Ribosomal reads from RNAseq reads.  
The process is based on two tools Trimmomatic and STAR for trimming and mapping respectively. First of all  
please make sure how to run these two tool on your system.  

On Sapelo/UGA, it is like this:
Trrimmomatic
```
module load ava/jdk1.8.0_20 
/usr/local/apps/trimmomatic/0.33/trimmomatic-0.33.jar 
```
STAR
```
modle load java/jdk1.7.0_67
/usr/local/apps/star/latest/bin/STAR
```

Dependency:
Trimmomatic  http://www.usadellab.org/cms/?page=trimmomatic
STAR   https://github.com/alexdobin/STAR


# 0 Copy the script and generate index file for rRNA sequences
copy NGSclean directory to your system, for example
NGSclean=/lustre1/lxue/NGSclean

prepare rRNA reference for STAR
'''
cd $NGSclean
mkdir rRNA_ref
module load java/jdk1.7.0_67
/usr/local/apps/star/2.4.2a/bin/STAR \
  --runThreadN 2  \
  --runMode genomeGenerate  \
  --genomeDir rRNA_ref  \
  --genomeChrBinNbits  5 \
  --genomeFastaFiles rRNA_only_NR.fas  
'''



# 1. Prepare design file
'''
cd /lustre1/lxue/RNAseq/01clean
NGSclean=/lustre1/lxue/NGSclean
python $NGSclean/generate_design_file.py -f /lustre1/lxue/RNAseq/00reads -d RNAseq_design.txt -p 

Check the design file. Revise the sample names if necessary, it will be used as a prefix for fastq files
more RNAseq_design.txt 
'''

# 2. Trim and Clean

```
cd /lustre1/lxue/RNAseq/01clean
settings
NGSclean=/lustre1/lxue/NGSclean
trimmoFull=/usr/local/apps/trimmomatic/0.33/trimmomatic-0.33.jar 
starFull=/usr/local/apps/star/latest/bin/STAR
adaptor=/usr/local/apps/trimmomatic/latest/adapters/TruSeq3-PE.fa
trimmo_module=java/jdk1.8.0_20 
star_module=java/jdk1.7.0_67
```

generate shell files  

-s       how to treat the singleton reads:merge, keep , discard  
-q       queue to run the jobs: queue(defualt), inter(run interactive jobs)  

```
python $NGSclean/trim_and_clean.py -d RNAseq_design.txt -t 8  -s merge \
  --run_trimmomatic $trimmoFull --load_trimmo_module $trimmo_module  --adaptor $adaptor  \
  --run_star $starFull --load_star_module $star_module 

python $NGSclean/trim_and_clean.py -d RNAseq_design.txt -t 8  -s keep \
  --run_trimmomatic $trimmoFull --load_trimmo_module $trimmo_module  --adaptor $adaptor  \
  --run_star $starFull --load_star_module $star_module 

python $NGSclean/trim_and_clean.py -d RNAseq_design.txt -t 8  -s discard -q inter \
  --run_trimmomatic $trimmoFull --load_trimmo_module $trimmo_module  --adaptor $adaptor  \
  --run_star $starFull --load_star_module $star_module 
```



