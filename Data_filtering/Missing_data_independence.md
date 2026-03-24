# Filtering loci based on the independence of missing data

### For contemporary Ne analysis, some amount of missing data can be tolerated, but only if the missing data are independent of allele frequencies.

### One way to assess this is to test loci for Hardy-Weinberg equilibrium at different missing data thresholds.

### The data consists of 1088 SNPs, and has already undergone some missing data filtering, where loci with more than 15% missing data have been excluded. We are now going to compare *G*is and p-values between a 20% and 30% missing data threshold for individual turtles.

### I do this for each temporal sample set, from the 1990s and the 2010s.

### There are many population genetics packages that you can use to compute *G*is, and calculate a HWE p-value, but in this instance I used [GenoDive](https://www.patrickmeirmans.com/software/Home.html). 

### After I've calculated all my statistics from the different missing data treatments I've compiled a [table](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/data/HWE.txt) that I read in R like this...

    head(H)
             Loci    N30    T30    O30  N30p  T30p  O30p    N20    T20  N20p  T20p
    1 01_11028227  0.209  0.281  0.236 0.251 0.217 0.071  0.189  0.277 0.276 0.228
    2 01_23416158 -0.070  0.370  0.132 0.793 0.119 0.309 -0.088  0.364 0.755 0.151
    3 01_23416385  0.000  0.000  0.000 1.000 1.000 1.000  0.000  0.000 1.000 1.000
    4 01_24898306  1.000  0.000  0.728 0.023 1.000 0.023     NA     NA    NA    NA
    5 01_28115838 -0.059  0.359  0.166 0.539 0.033 0.101 -0.107  0.305 0.410 0.068
    6 01_38259617     NA -0.014 -0.013    NA 0.985 0.985     NA -0.016    NA 0.985

### N columns are the 90s turtles, while T columns are the 2010s turtles, and O columns are combined.

### 30 refers to 30% missing data. 20 refers to 20% missing data. And column names with a p are the p-values.

### Next we create linear models for each temporal sample set, and plot.

    HWE90s_model <- lm(N30 ~ N20, data = H)
    HWE10s_model <- lm(T30 ~ T20, data = H)
    plot(HWE90s_model)
    plot(HWE10s_model)

### Then we isolate the residuals

    HWE90s_res <- rstandard(HWE90s_model)
    HWE10s_res <- rstandard(HWE10s_model)

### Looking at the plots it appears that plus or minus two standard residuals is an appropriate threshold for exclusion.

### How many loci exceed the threshold in each temporal sample?

    > length(which(abs(HWE90s_res) > 2))
    [1] 36
    
    > length(which(abs(HWE10s_res) > 2))
    [1] 33

### And how many loci overlap the two sets?

   

    bad_loci90s <- names(HWE90s_res[which(abs(HWE90s_res) > 2)])
    bad_loci10s <- names(HWE10s_res[which(abs(HWE10s_res) > 2)])
    
    intersect(bad_loci90s, bad_loci10s)
        [1] "26_13486177"

### Only 1? I thought there'd be more than that, but we're still going to remove all these loci from both data sets.

    all_bad_loci <- unique(c(bad_loci90s,bad_loci10s))
    newH <- H[!row.names(H) %in% all_bad_loci,]

### Now, what did that do? What this even necessary?
### Let's do before and after plots to assess the [results](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/Data_filtering/Gis_before_and_after.pdf).

### There most certainly was a subset of SNPs that were missing data that were not fully independent. Overall, we only removed 68 out of 1088 loci and have cleaner data. 
