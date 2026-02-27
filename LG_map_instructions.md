
# Linkage group map instructions

## Here I go over how to group loci according to linkage order, including our method for binning loci together. This approach did not work well for the leatherback turtles, probably because we only had 28 individuals total (a small map) and because background linkage among chromosomes was high in this population, but future users might find this approach helpful.

### We are starting with large, whole-genome paired-end fastq files and trimming them with the program [fastp](https://github.com/OpenGene/fastp)

    fastp --trim_poly_g --length_required 100 -i R1.fastq.gz -I R2.fastq.gz -o R1_trimmed.fq.gz -O R2_trimmed.fq.gz

Before trimming: 117565414 reads R1, sequence length 150
After trimming:  108362002 reads R1, sequence length 100-150

### Then we index the Dermochelys coriacea genome and map short reads with the program [minimap2](https://github.com/lh3/minimap2)

    minimap2 -t 8 -d Dcor_genome_ref.mmi GCF_009764565.3_rDerCor1.pri.v4_genomic.fna

    minimap2 -t 8 -ax sr Dcor_genome_ref.mmi X_R1.fastq.gz X_R2.fastq.gz > X.sam

### Here's a batch script in Python that you can use, where the batchF and batchR are the R1 and R2 batch files respectively. I use something similar for fastp as well.

    with open("batchF.txt","r") as forwards, open("batchR.txt","r") as reverses:
       for f,r in zip(forwards,reverses):
          os.system("minimap2 -t 8 -ax sr Dcor_genome_ref.mmi %s %s > %s.sam" %(f.rstrip("\n"), r.rstrip("\n"), f.lstrip("/pathway to the directory holding fastq files/").rstrip("_R1.fastq.gz\n")))

### The files are big, so even with multithreading (-t) this can take a long time.

### remove unmapped reads and convert to bam

    samtools view --threads 10 -b -F 4 X.sam > X.bam

### then you need to run the command samtools fixmate with the -m option (this is necessary for removing PCR duplicate reads)

    samtools fixmate -m -@ 14 X.bam X_fxm.bam

### In my experience some samtools functions don't work well with standard out, such as fixmate and markdup, so it is hard to do these in a batch

### Sort the .bam 

    samtools sort --threads 8 X_fxm.bam > X_sortedFxm.bam

### Remove PCR duplicates, using the samtools markdup 

    samtools markdup -@ 14 -r -s -d 100 X_sortedFxm.bam X_rmdup.bam

### Check the coverage for all individuals across chromosomes

    samtools coverage -m X_rmdup.bam

    NC_050068.2 (354.45Mbp)  # Chromosome 1
    >  89.70% │▅▅▆▆▇▇█▇▇▇███▇█ █▇█▇▇█▇█████▇███████▇███▇▇██████▆████▇███▅▇█████████▇████████████████▆▇▃▇█████▅▆█▇▆█▇██████▇▇██▇█▃█████████▇███▇█▇▇███████▇▃▆▆│ Number of reads: 29775497
    >  79.73% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ 
    >  69.77% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Covered bases:   350.17Mbp
    >  59.80% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Percent covered: 98.79%
    >  49.83% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Mean coverage:   12.5x
    >  39.87% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Mean baseQ:      42.2
    >  29.90% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Mean mapQ:       56.2
    >  19.93% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ 
    >   9.97% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo bin width: 2.50Mbp
    >   0.00% │██████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████████│ Histo max bin:   99.666%
              1       24.96M    49.92M    74.88M    99.85M   124.81M   149.77M   174.73M   199.69M   224.65M   249.61M   274.58M   299.54M   324.50M     354.45M

#####  Mean coverage across turtles ranged from 8x to 35x with most falling in the 12-18x range

### Now use bcftools to call SNPs, piping the "mpileup" and "call" scripts together

    bcftools mpileup \
    --skip-indels \
    -q 10 -Q10 \
    --annotate FORMAT/AD,FORMAT/DP \
    -f GCF_009764565.3_rDerCor1.pri.v4_genomic.fna \
    --bam-list the_bams.txt \
    --threads 10 | \
    bcftools call \
    -o BCF2_raw.vcf \
    -O z -mv \
    --threads 10

##### All .bam files are listed in "the_bams.txt"

### We have 3,761,592 raw SNPs, or about 1 SNP for every 585 bp on average.

    vcftools --vcf BCF_raw.vcf
    After filtering, kept 28 out of 28 Individuals
    After filtering, kept 3761529 out of a possible 3761529 Sites
        
    > 2200000000 / 3760000
    [1] 585.1064

### Filter for missing data in [vcftools](https://vcftools.github.io/index.html)

    vcftools --vcf BCF_raw.vcf --max-missing 0.9 --recode --out BCF2_filt1
    After filtering, kept 3402269 out of a possible 3761529 Sites

##### For some projects you might be want to tolerate more missing data, but here we only have 28 individuals, so I am only going to allow sites to miss data 10% of the time, or in no more than 2 individuals.

### Now examine variants for depth

    vcftools --vcf BCF2_filt1.recode.vcf --site-mean-depth --out BCF2_depth

### In R, examine the distribution of mean read depth

    > hist(d$VAR_DEPTH, breaks=50)
    > hist(d$VAR_DEPTH[which(d$VAR_DEPTH < 100)], breaks=50)
    > hist(d$VAR_DEPTH[which(d$VAR_DEPTH < 40)], breaks=50)
    > hist(d$VAR_DEPTH[which(d$VAR_DEPTH < 10)], breaks=50)
    
    > quantile(d$VAR_DEPTH[which(d$VAR_DEPTH < 100)],0.95)
        95% 
    45.4339 
    > quantile(d$VAR_DEPTH[which(d$VAR_DEPTH < 100)],0.05)
         5% 
    4.74074  

##### The distribution here looks good, with a normal read depth ranging between 5 and 45

### Filter for depth in vcftools

    vcftools --vcf BCF2_filt1.recode.vcf --min-meanDP 5 --max-meanDP 45 --recode --out BCF2_filt2
    After filtering, kept 3032210 out of a possible 3402269 Sites

### Let's take a quick look at how the SNPs are distributed across chromosomes

cut -f 1 BCF2_filt2.recode.vcf | sort | uniq -c

     479118 NC_050068.2
     356688 NC_050069.1
     276187 NC_050070.1
     196492 NC_050071.1
     183217 NC_050072.1
     181312 NC_050073.1
     175123 NC_050074.1
     146282 NC_050075.1
     136866 NC_050076.1
     113625 NC_050077.1
     106096 NC_050078.2
      62086 NC_050079.1
      59814 NC_050080.1
      81656 NC_050081.1
      45800 NC_050082.1
      38583 NC_050083.1
      37520 NC_050084.1
      32800 NC_050085.1
      30136 NC_050086.2
      47322 NC_050087.2
      31649 NC_050088.1
      30793 NC_050089.1
      37190 NC_050090.1
      32789 NC_050091.1
      26365 NC_050092.1
      28118 NC_050093.1
      28816 NC_050094.1
      26194 NC_050095.1
        462 NW_025091561.1
         14 NW_025091562.1
          1 NW_025091563.1
          3 NW_025091564.1
       1306 NW_025091565.1
        152 NW_025091566.1
        348 NW_025091567.1
        103 NW_025091568.1
        967 NW_025091569.1
        114 NW_025091570.1
         48 NW_025091571.1
         55 NW_025091572.1

##### The 28 autosomal chromosomes all start with the prefix "NC_" anything with a prefix of "NW_" is a non-localized scaffold. When I began this project, I had hoped to be able to place them using linkage, improving the genome assembly, but none of the loci on these unlocalized scaffolds were very informative.

### Before we proceed, we should calculate relatedness in vcftools.

    vcftools --vcf BCF2_filt2.recode.vcf --relatedness2 --out STX2_rel

##### All previously inferred relationships are confirmed using the high-density genomic SNPs.


### Now for [Lep-MAP3](https://sourceforge.net/p/lep-map3/wiki/LM3%20Home/)

### This next line of code formats a .vcf file for input into lep-MAP3 and runs the module ParentCall2 (All Lep-MAP3 modules are java-based)

    cat BCF2_filt2.recode.vcf | gawk '/^#/{print}/^[^#]/{m=0;for (i=10;i<=NF;++i) if ($i ~ /\.\/\./) ++m; if (m<=26) print}'| \
    java -cp lepmap3/bin ParentCall2 data=ped_STXf2.txt vcfFile=- removeNonInformative=1 ignoreParentOrder=1|gzip > data.call.gz
    Number of called markers = 1652486 (1652486 informative)

##### On some computers awk is the default, but on mine gawk is what works best on my operating system.

##### We already filtered for missing data, but the m variable in the awk script will remove a certain count (those that match the pattern ./.) of missing data. In this case, over 26/28 missing individual genotypes.

##### The ParentCall2 module is removing non-informative markers, formatting the data for the next step, and imputing genotypes for missing individuals in the pedigree file (ped_STXf2.txt). In this case we are imputing all the paternal genotypes from offspring alleles, because we only have nesting mother samples.

##### Below is the pedigree file. If you copy and paste this as a template for your own run, be sure that it is tab-delimited in your text file, otherwise the program will throw an error.

    CHROM	POS	STX1	STX1	STX1	STX1	STX1	STX1	STX1	STX1	STX1	STX1	STX1	STX1	STX2	STX2	STX2	STX2	STX2	STX2	STX2	STX2	STX2	STX2	STX2	STX2	STX3	STX3	STX3	STX3	STX3	STX3	STX3	STX3	STX3	STX3
    CHROM	POS	P1	14395	P2	24210	96249	96239	96240	96245	96246	96247	96251	96252	P1	14395	P3	23560	96540	96541	96542	96543	96544	96545	96546	96547	P4	88837	85952	85953	85954	85955	85956	85960	85961	86024
    CHROM	POS	0	0	0	P1	P2	P2	P2	P2	P2	P2	P2	P2	0	0	0	P1	P3	P3	P3	P3	P3	P3	P3	P3	0	0	P4	P4	P4	P4	P4	P4	P4	P4
    CHROM	POS	0	0	0	14395	24210	24210	24210	24210	24210	24210	24210	24210	0	0	0	14395	23560	23560	23560	23560	23560	23560	23560	23560	0	0	88837	88837	88837	88837	88837	88837	88837	88837
    CHROM	POS	1	2	1	2	0	0	0	0	0	0	0	0	1	2	1	2	0	0	0	0	0	0	0	0	1	2	0	0	0	0	0	0	0	0
    CHROM	POS	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0

##### The first line of the file is the family. Here we have three families (STX1, STX2, STX3) and note that the grandparents (P1, and 14395) are present in both STX1 and 2. That is, STX1 and 2 are three-generation families, where the mothers are full sisters.

##### STX3 is a two-generation family with an unrelated mother, but curiously she seems to be the aunt of 24210's offspring (i.e. P2's sister). That doesn't affect what we are trying to do with the linkage map.
 
##### The second line is the individual, starting with the grandparents (P1 and 14395), then parents (e.g. P2, 24210), then the offspring (e.g. 96249).
##### The third line is the father of each individual
##### The fourth line is the mother of each individual
##### The fifth line is the sex of each individual (male=1, female=2, unknown-0).
##### The sixth line is for adding a phenotype. Zero if none.

### If you don't want to bin your SNPs together, then you can proceed to the SeparateChromosomes 2 module and run it normally, according to the software documentation.

### Binning SNPs together can be done as follows.

    zcat data.call.gz|awk 'BEGIN{print "#binned markers"}(NR>7){if (prev==0 || $2-prev >= 10000 || $1 != prev1) {++n;prev=$2;prev1=$1}; print n}' > bin.txt

 ##### The bin width is 10,000 bp. Or in other words, that is the maximum distance between two SNPs in the same bin.

### Now use the OrderMarkers2 module to bin the loci

    zcat data.call.gz|java -cp lepmap3/bin OrderMarkers2 map=bin.txt data=- phasingIterations=3 outputPhasedData=3 improveOrder=0 recombination1=0 recombination2=0 2>///dev/null > om_binned.txt

##### Don't be confused by the pathway "lepmap3/bin" that is a directory where the java files are, and has nothing to do with binning.

### Next, we use a custom awk script (all awk scripts herein are written by Pasi Rastas) to complete the binning process

    awk -vpedigree=1 -f order2data.awk om_binned.txt | gzip > data_binned.gz

##### This script is found in the same Github repository as these instructions.

### With the SNPs binned, we are now ready for the SeparateChromosomes2 module

    zcat data_binned.gz | java -cp lepmap3/bin SeparateChromosomes2 data=- lodLimit=6.5 numThreads=12 lod3Mode=1 distortionLod=1 > binned_map6d.txt

##### Setting lod3Mode=1 considers all binned SNPs as a block.

##### The "lodLimit" parameter sets the statistical association threshold (the logarithm of the odds) between loci.
##### In a straightforward scenario, an lod score of three means that the odds of two loci being linked are 1,000 to 1, but setting a lower lod threshold for linkage mapping requires some exploration of the data.

##### A useful rule of thumb is that the lod score should be at least one tenth the number of offspring in your analysis. Here we have 24 offspring so I wouldn't go any lower than lod=2.5 here.

##### Once you've run SeparateChromosomes2, use the following shell script to see how markers are sorting into linkage groups

    sort binned_map5.txt|uniq -c|sort -r -n | less -S
     115757 0
      53206 1
        143 2
         77 3
         23 4
         18 5
         11 7
         11 6
         10 8
          9 9

##### Here's a bad example. 
##### 115,757 SNPs were singles in the 0 group that didn't associate with any linkage groups at all.
##### 53,206 SNPs were in linkage group 1, and only 143 in linkage group 2. 
##### This suggests that the lod threshold was too low, and that a lot of SNPs are being agglomerated into a single big linkage group.

##### This is a better example...

     127244 0
       6876 1          
       5774 2          
       3572 5
       3499 4
       3145 6
       2756 7
       2545 3
       2486 10
       2393 8          
       2386 9
       1799 11
       1110 14
       1048 13
        871 12
        789 16
        637 15
        540 20
        528 22
        472 19
        414 17
        401 23
        400 27        
        369 18
        353 21
        303 24
        299 26
        255 25
         74 28

##### We still have a lot of unassociated SNPs that didn't fall into any linkage group, but the SNP counts for the following linkage groups are more even, and correspond better to the size of each Chromosome in BP (keep in mind that we have microchromosomes in this species). 

##### After linkage group 28 the loci numbers drop dramatically, but that's what we would expect anyways, because there are only 28 chromosomes.

### Once the LOD threshold has been fine tuned, you will run the JoinSingles2All module, which will attempt to add more loci to the 28 linkage groups from the pool of unassociated SNPs

    zcat data_binned.gz | java -cp ../lepmap3/bin JoinSingles2All data=- lodLimit=4 lod3Mode=1 map=binned_map6d.txt numThreads=4 lodDifference=1 distortionLod=1 > map_js4.txt

##### For this module you will want to relax the lodLimit to allow more loci to associate
##### It can be a good idea to run JoinSingles2All multiple times using the output (map_js4.txt) from the previous run.

### Next you want to order the markers in the linkage groups with the OrderMarkers2 module.

##### This is best done in a loop that iterates across the number of linkage groups that correspond to the correct number of chromosomes.

    for i in {1..28}; do zcat data_binned.gz | java -cp lepmap3/bin/ OrderMarkers2 map=map_js4.txt data=- chromosome=$i > order$i.txt; done

##### This is another module that can, and maybe should, be run more than once. Subsequent runs look like this...

    for i in `ls | grep order`; do zcat data_binned.gz | java -cp lepmap3/bin/ OrderMarkers2 numThreads=8 evaluateOrder=$i improveOrder=1 data=- > B_${i}; done

##### Instead of looping through a set of numbers I am looping through the output files (with "order" in the file name) from the last run, and adding a "B_" onto the front of the output file name to differentiate it from the original.

##### Sometimes it works better to type: `evaluateOrder=${i}`

### To assess ordering, we generate .dot format files

    for i in {1..28}; do java -cp lepmap3/bin/ LMPlot order$i.txt > order$i.dot; done

### We read these in R and plot with the function grViz

    library(DiagrammeR)
    grViz("order1.dot")
    
##### This will give you a little chain-like diagram. Everything should look linear, but if there are small (or large) reticulations in the chain then you may need to order the markers some more, or remove the problematic loci.

### The process of separating and ordering chromosomes doesn't preserve the scaffold and base-pair position of each SNP, so we have to go back and add that information to the results

    zcat data.call.gz|cut -f 1,2|awk '(NR>=7)' > snps.txt

    paste snps.txt bin.txt|awk '{if (!($3 in d)) {d[$3]; print}}' > snps_binned.txt

    for i in {1..28}; do awk -vchr=$i -f map2genotypes.awk order$i.txt | awk -vcolumn=1 -f map.awk pass=1 <(awk '{print NR-1"\t"$0}' snps_binned.txt) pass=2 - > ord$i.mapped; done

##### The two awk scripts (map.awk and map2genotypes.awk) can be found in the main github repository

### Next we want to purify the linkage groups

##### In theory, if you have a chromosome-level genome assembly, you should be able to look through your .mapped files and see only the same chromosome for all loci.

##### However, for the leatherback linkage map the results were a bit noisy, and the linkage groups were not pure chromosomes. For example, if we look at chromosome 1.

    cat ord1.mapped | grep "NC_050068.2" > pure1.mapped
    wc ord1.mapped
      10629  573966 1481839 ord1.mapped
    wc pure1.mapped
       9962  537948 1388586 pure1.mapped

##### Using grep to find only loci on NC_050068.2 removed 667 SNPs from the original .mapped file.

##### All linkage groups need to modified in this way to obtain pure chromosome linkage groups.

### How many loci are in our pure map?

      for i in {1..28}; do wc pure$i.mapped; done
       
     9962
     8955
     3657
     5465
     5716
     4987
     4483
     3739
     3451
     3682
     2493
     1348
     1638
     1630
     999
     1023
     643
     698
     670
     777
     452
     735
     542
     437
     356
     392
     550
     113

##### That's 69,593 loci total

##### This is what the first six columns of the .mapped file lines looks like

    NC_050070.1	2249	49394	3	0.000	0.000
##### The first column is the chromosome
##### The second is the bp position
##### The third is the SNP index
##### The fourth is the linkage group number
##### The fifth and sixth are the male and female map cM values––both zeros here because this is the beginning of the file.

##### Sometimes the bp are inverted relative to cM in the output, like here on linkage group 4

    NC_050071.1	144665679	77809	4	0.000	0.000
...

    NC_050071.1	10304876	66989	4	58.333	100.000	


 ##### That's not a problem. What really matters are the cM distances between loci. The descending or ascending order can be reversed.

    new_cM = max(cM) - cM

##### Something else to notice. The female map is longer than the male map. This is typical and is a good sign. Look at the comparisons of the total lengths between the male and female maps across all chromosomes...

    > sum(mal_cM)
    [1] 1999.998
    > sum(fem_cM)
    [1] 2516.667

##### This pattern is consistent across sexually reproducing organisms. This is even true for species like sequentially hermaphroditic fish with individuals that reproduce as both males and females. The male map always has less recombination.

##### But, for most applications you'll want to average the two recombination rates together, to get the sex-averaged rate, and sex-averaged cM positions of SNPs.

### The last thing to do is to plot the cM values against the bp positions. This is called a Marey map.

![Marey Maps for Linkage groups 9 and 10](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/Ch9Ch10.pdf)

##### cM values on the y-axis. Bp values on the x-axis

##### Viewing the raw marey maps will help you trouble shoot any issues that your map may have. In this case the recombination landscape of the chromosomes is a bit fuzzy and needs some more polishing.

##### There are tools that can help remove unwanted noise, but outliers can also be removed manually when needed.












