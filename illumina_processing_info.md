# Working codes for processing illumina sequencing files for Daniele's group:'
### This is what I used to process these files.

#### 1. Run fastqc and multiqc on raw reads ###############################################

This is done in Colonial One and assumes you are in the directory with all the sequencing directories (they start with L4##). You need to be in an interactive node to run fastqc, like we have talked about before.

```
module load fastqc
fastqc L4*/*.fastq.gz

module unload python
module load miniconda3/4.3.27.1
conda activate multiqc
multiqc L4*/*fastqc.zip
```

#### 2. Run fexbar first, then prinseq ####################################################
Run flexbar with illumina adapters (NexteraPE-PE.fa).

```sbatch flexbar.sh```

Run flexbar with IonTorrent adapters (adapters_SE.fa). You need to alter the script to change the adapter file.

```sbatch flexbar.sh```

Run prinseq.

```sbatch prinseq.sh```

Remove bad files and singletons.

```
rm L4*/*bad*.fastq
rm L4*/*singletons*.fastq
```

Run fastqc.

```
module load fastqc
fastqc L4*/fbIT_*prinseq_good*.fastq
```

Run multiqc.

```
module unload python
module load miniconda3/4.3.27.1
conda activate multiqc
multiqc L4*/fbIT_*prinseq_good*fastqc.zip
```

#### 2. Align cleaned reads to hg19 with bowtie2 ##########################################

Run bowtie2.

```sbatch bowtie2.sh```

The base of the script looks like this:
```
for name in $(cat samps.txt); do
     sampname=$(ls ${name}/L4*_R1_*.fastq.gz | cut -d"-" -f3);
     bowtie2 -p $(nproc) \
 	--sensitive \
 	--rg PL:MiSeq \
 	--rg-id ${sampname} \
 	--rg SM:${sampname} \
 	-x refs/genome \
 	-1 $(ls ${name}/flexbar_1_prinseq_good_*.fastq) \
 	-2 $(ls ${name}/flexbar_2_prinseq_good_*.fastq) \
 	-S $name/bowtie2/aligned.sam 2> $name/bowtie2/bt2_log.txt;
 	echo $name
done
```

Convert sam to bam.

```
module load samtools
for f in L4*/bowtie2; do samtools view -Sb $f/aligned.sam > $f/aligned.bam && echo $f; done
```

Sort bam file

```for f in L4*/bowtie2; do samtools sort -o $f/sorted.bam $f/aligned.bam && echo $f; done```

Index bam file.

```for f in L4*/bowtie2; do samtools index $f/sorted.bam && echo $f; done```

#### 3. Pull files from colonial one to computer ##########################################
```
rsync -av kmgibson@login4@colonialone.gwu.edu:/lustre/groups/cbi/kgibson/FORENSICS/illumina/20190222_Podini_rerun/samps.txt .
for f in $(cat samps.txt); do mkdir $f; done
for f in L*; do
	rsync -av kmgibson@login4.colonialone.gwu.edu:/lustre/groups/cbi/kgibson/FORENSICS/illumina/20190222_Podini_rerun/$f/bowtie2 $f/
done

for f in L*; do
    rsync -av $f/bowtie2/sorted.bam* ~/Box/FONTANA/20190222_Podini_illumina/$f/
    echo $f
done
```

#### 4. Installing samtools on your local computer to view the bam files if need be ##############

Download the tar file for samtools at http://www.htslib.org/download/

```
cd Downloads/samtools-1.9
./configure --prefix=/where/to/install #(I did /User/keyliegibson)
make
make install
export PATH=/where/to/install/bin:$PATH
```

Then you should be able to call samtools as such:

```samtools view sorted.bam | head```
