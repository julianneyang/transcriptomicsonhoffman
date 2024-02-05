# transcriptomicsonhoffman

After you log in to Hoffman2 and request a computational node: 

## Preprocessing the data 

Assuming your fastq files are in your current working directory:

1. Install FastQC (here we create a new conda env to install fastqc)
```bash
conda create -n fastqc fastqc
```

```bash
conda activate fastqc
```

2. Run FastQC on your file (I will later edit this to become a job submission script) :
```bash
mkdir 1_QC_output/
fastqc *.fastq.gz -o 1_QC_output/
```

3. Aggregate quality reports for all samples by using multiQC (note: for some reason I had issues with forcing multiqc to use python 3.10 so I had to use the below workaround. MultiQC takes as input a directory full of report.html files.

For downloading MultiQC, do not use conda, it downloads an outdated version. Instead I used pip to install the development version, and I also forced installed to $PROJECT which has enough space as opposed to the default $HOME installation

```bash
pip install --upgrade --force-reinstall git+https://github.com/MultiQC/MultiQC.git -t /u/project/jpjacobs/jpjacobs/rna_seq/
```
You may need to find the exact filepath to multiqc via the following command:
```bash
which multiqc
```
Replace ~/.local/bin/multiqc with the exact filepath:

```bash
python ~/.local/bin/multiqc ./
``` 

4. Copy the .html report over to your local directory with `scp` or push to Github from Hoffman. open report.html in a browser. For help interpreting multiqc results, see the following resoureces:

5. Trim adapters and low-quality reads with Trimmomatic. Since we already have trimmomatic installed in the kneaddata env, we are going to activate the kneaddata env:
```bash
conda activate kneaddata
```
```bash
trimmomatic PE JJ1715_393_S43_R1_001.fastq.gz JJ1715_393_S43_R2_001.fastq.gz output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:/u/home/j/jpjacobs/project-jpjacobs/software_rna_seq/Trimmomatic/trimmomatic-0.39/adapters/TruSeq3-PE.fa:2:30:10:2:True LEADING:3 TRAILING:3 MINLEN:36
```

## References: 
https://bookdown.org/jean_souza/PreProcSEQ/quality-control.html#fastqc-1 
https://github.com/hbctraining/Intro-to-rnaseq-hpc-gt/blob/master/lessons/08_rnaseq_workflow.md
