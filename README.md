<h1>Complete Genomic Characterization of Global Pathogens, Respiratory Syncytial Virus (RSV), and Human Norovirus (HuNoV) Using Probe-based Capture Enrichment.</h1>
 <h3>Preprint: https://www.biorxiv.org/content/10.1101/2024.09.16.613242v1</h3>

<details>
 <summary><h2>1. Pre & post capture viral assemblies</h2></summary>
 To evaluate the capability of the capture methodology to assemble full-length genomes, the VirMAP:https://github.com/cmmr/virmap pipeline was used to reconstruct RSV and HuNoV genomes. Trimmed and host-filtered reads were processed through VirMAP (24) to assemble complete RSV or HuNoV genomes. The VirMAP summary statistics include information on reconstructed genome length, the number of reads mapped to the reconstruction, and the average coverage across the genome.
</details>

<details>
 <summary><h2>2. Figure 2 depicts viral read recovery efficiency.</h2></summary>
Viral read recovery efficiency. Percent of trimmed, non-human sequence reads (post-processing) that mapped to the target viral genome in pre-capture (circles) and post-capture (triangles) libraries. CT value range of samples: ‘CT <20’ (red), ‘CT 20 to 30’ (light blue), ‘CT > 30’ (green) & ND (not detected) (pink). A: Viral reads mapping to RSV genomes, split by two subtypes. B: Viral reads mapping to HuNoV genomes, split by genotypes (GI.1, GII.4, Other GII).
</details>

<details>
 <summary><h2>3. Figure 3 compares genome coverage post-capture vs pre-capture.</h2></summary>
Average genome coverage obtained in post-capture (triangles) and pre-capture (circles) samples. Genome reconstruction was classified as follows: ‘complete’ (within expected length range, >90% completeness & >20x coverage), ‘complete with low coverage’ (within expected length range, >90% completeness & <20x coverage), or ‘incomplete’ (below expected length range, <90% completeness & <20x coverage). CT value range of samples: ‘CT <20’ (red), ‘CT 20 to 30’ (light blue), ‘CT > 30’ (green) & ‘ND’ (pink). A: RSV samples split by RSV-A or RSV-B genotype. B: HuNoV samples split by five genotypes (GI.1, GII.4, Other GII).
</details>

<details>
 <summary><h2>4. Figure 4 shows the breadth of coverage (20x) for each sample. </h2></summary>

 
 <h3>Breadth of 20x coverage</h3>

To calculate the breadth of coverage, we first align the reads to a given reference genome (see below), and then use `samtools depth` to calculate the coverage at each base across the genome.

For the alignments, we used `bwa mem` and different reference genomes depending on the virus. For RSV, we used the RSV/A and RSV/B reference genomes that were recently published by our group, which can be found [here](https://academic.oup.com/ve/article/10/1/vead086/7503540). For Norovirus, we used the assembled genome from each sample (assembled using capture probes) as a reference. To ensure quality, we applied a filter for a minimum mapping quality of 20 Phred scores (`-q 20`) when calculating the coverage.

Here’s the code we used for the alignment and coverage calculation:

 ```
# Performing alignment for each sample. The samtools commands will convert the output to bam and immediatelly sort the output into the final sorted file.
bwa mem -t 4 -T 0 reference read1 read2 | samtools view -hb - | samtools sort -o $outputdir/${name}.sorted.bam -

# Calculating the breadth of coverage for 20x and 30x
cov20=$(samtools depth -q 20 $outputdir/${name}.sorted.bam | awk '$3 >= 20 {count++} END {print count}')
cov30=$(samtools depth -q 20 $outputdir/${name}.sorted.bam | awk '$3 >= 30 {count++} END {print count}')
```
Where:
`reference`: is the reference genome ;
`read1`: the fastq file containing reads 1 ;
`read2`: the fastq file containing reads 2 ;
`outputdir`: the output directory ;
`name`: the sample name.

<h3> Plotting the breadth of coverage (20x)</h3>
The plots were created using the R script uploaded in the `scripts/fig4` subfolder of this Github page. The table containing the calculated breadth fof coverage (20x) is also found in that folder. For calcuating the percentages of genomes covered by 20x coverage, we used the reference genome lengths for RSV (15243bp), or an average of the assembled genome lengths for all samplesin the case of NoV (7526.148148bp).


Info about the RSV reference genomes here: https://doi.org/10.1093/ve/vead086
</details>


<details>
 <summary><h2>5. Figure 5 </h2></summary>

</details>
