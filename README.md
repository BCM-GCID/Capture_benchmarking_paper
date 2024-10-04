# Complete Genomic Characterization of Global Pathogens, Respiratory Syncytial Virus (RSV), and Human Norovirus (HuNoV) Using Probe-based Capture Enrichment.
 - _Preprint: https://www.biorxiv.org/content/10.1101/2024.09.16.613242v1_

<details>
 <summary><h2>1. Viral reads from pre & post capture data sets were were assembled and read stats were generated using [Virmap](https://github.com/cmmr/virmap) </h2></summary>

</details>

<details>
 <summary><h2>2. Figure 2 deptics viral read recovery efficiency.</h2></summary>

</details>

<details>
 <summary><h2>3. Figure 3 compares genome coverage post-capture vs pre-capture.</h2></summary>

</details>

<details>
 <summary><h2>4. Figure 4 shows the breadth of coverage (20x) for each sample. </h2></summary>

 
 <h3>Breadth of 20x coverage</h3>

To calculate the breadth of coverage, we first align the reads to a given reference genome (see below), and then use `samtools depth` to calculate the coverage at each base across the genome.

For the alignments, we used `bwa mem` and different reference genomes depending on the virus. For RSV, we used the RSV/A and RSV/B reference genomes that were recently published by our group, which can be found [here](https://academic.oup.com/ve/article/10/1/vead086/7503540). For Norovirus, we used the assembled genome from each sample (assembled using capture probes) as a reference. To ensure quality, we applied a filter for a minimum mapping quality of 20 Phred scores (`-q 20`) when calculating the coverage.

Here’s the code we used for the alignment and coverage calculation:

 ```
# Performing alignment for each sample
bwa mem -t 4 -T 0 reference read1 read2 | samtools view -hb - | samtools sort -o $outputdir/${name}.sorted.bam -

# Calculating the breadth of coverage for 20x and 30x
cov20=$(samtools depth -q 20 $outputdir/${name}.sorted.bam | awk '$3 >= 20 {count++} END {print count}')
cov30=$(samtools depth -q 20 $outputdir/${name}.sorted.bam | awk '$3 >= 30 {count++} END {print count}')
```
Where:
`reference`: is the reference genome
`read1`: the fastq file containing reads 1
`read2`: the fastq file containing reads 2
`outputdir`: the output directory
`name`: the sample name

Info about the RSV reference genomes here: https://doi.org/10.1093/ve/vead086
</details>
   
