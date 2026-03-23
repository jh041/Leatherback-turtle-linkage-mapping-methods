# Predicting cM from bp position

### With cM values for all linkage mapping SNPs, you can examine the recombination landscape of chromosomes using [marey maps](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/Chromosome_linkage_mapping/marey.png)

### Patterns of recombination in the genome are informative to evolutionary biology and molecular ecology, but not all SNPs from your population- or species-level studies are going to be present in the map. Therefore, it is beneficial to be able to predict cM values for any new polymorphic variants.

### This can be accomplished from the marey maps by using locally weighted scatterplot smoothing (LOWESS) to model the relationship between cM and bp.

### I am going to do this in R

##### Step 1, load the "[map.txt](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/data/map.txt)" file resulting from the Chromosome-informed map into R

    > map <- read.table("map.txt", header=F)
    
    > head(map)
    
    V1 V2 V3 V4 V5
    1 NC_050068.2  10267  1  0  0
    2 NC_050068.2  50317  1  0  0
    3 NC_050068.2  86728  1  0  0
    4 NC_050068.2 111121  1  0  0
    5 NC_050068.2 135771  1  0  0
    6 NC_050068.2 148195  1  0  0

##### Columns are: chromosome, bp, linkage group, male cM, and female cM, respectively

### I am going to compute the sex-averaged cM values and use them as the basis for my model.

    sa_vals <- cbind.data.frame(map$V3, map$V2, rowMeans(map[, c("V4", "V5")], na.rm = TRUE))
    names(sa_vals) <- c("Chr","bp","cM")


    > head(sa_vals)
    
    Chr bp cM
    1 1  10267  0
    2 1  50317  0
    3 1  86728  0
    4 1 111121  0

### Make individual data frames for each chromosome 1-28

    sa_chr01 <- sa_vals[which(sa_vals$Chr==1),]
    ...
    sa_chr28 <- sa_vals[which(sa_vals$Chr==28),]

### Before we make each model, it's a good idea to plot it against the marey map, to see how well it fits.

    > library(ggplot2)
    
    > model_ch01_plot <- ggplot(sa_chr01, aes(bp, cM)) +
    geom_point(color="blue", alpha=0.1) +
    geom_line(data=as.data.frame(lowess(sa_chr01$bp, sa_chr01$cM, f = 0.05, iter = 1)), aes(x,y), color="orange", linewidth=1) + labs(title="Chromosome 1", x=NULL, y="cM") + theme_classic() + theme(plot.title = element_text(size = 12))

##### you can see that plot [here](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/Chromosome1_lowess.pdf).

##### you can adjust the fit of the model using the f and iter parameters.

##### You don't want to overfit the model too much, but the bigger risk is to underfit the model, especially at the chromosomal extremes.

##### In addition, the map ends often need to be edited,  as these regions typically have problems in read mapping and the support for the first/last crossover is always weaker than in the middle crossovers.

##### You can calculate how many markers are supporting the first and last crossovers (number of markers with the corresponding map position). 

##### If on average there are X markers at each map position, removing crossover with < 0.05X markers could be ok... (chance of removing a correct marker is 5%). This is especially important when ordering markers using linkage, as there are frequently SNPs at the extreme ends that stand apart from the rest of the map.  

### Once you are satisfied with your map, you can use that lowess() function, native to R, to create models for each chromosome.

    ch01_bp2cM_model <- lowess(sa_chr01$bp, sa_chr01$cM, f = 0.05, iter = 1)

...

    ch28_bp2cM_model <- lowess(sa_chr28$bp, sa_chr28$cM, f = 0.05, iter = 3)

##### Each model is a list object, with two elements $x (bp) and $y (cM)

##### Some of these models have negative elements.
##### This happens when the model is underfit at the front end and bp positions below the bottom of the model are extrapolated downwards.

##### So all nexative cM values in the models need to be converted to zero. In my experience this is only a few loci. If not, try to fit the model to the data better.

    ch01_bp2cM_model$y[which(ch01_bp2cM_model$y < 0)] <- 0

...

    ch28_bp2cM_model$y[which(ch28_bp2cM_model$y < 0)] <- 0

### Now that we've done that, save the models to file

    saveRDS(ch01_bp2cM_model, file = "ch01_bp2cM_model.rds")
...

    saveRDS(ch28_bp2cM_model, file = "ch28_bp2cM_model.rds")

##### To load these files again use this command...

    ch01_bp2cM_model <- readRDS("ch01_bp2cM_model.rds")

### Next, we can create a predictor function to predict values based on bp in other data sets

    ch01_bp2cM_predict <- approxfun(ch01_bp2cM_model)

...

    ch28_bp2cM_predict <- approxfun(ch28_bp2cM_model)

### Assume, for example, that you have some population genetic SNPs that you want to place on the linkage map. To predict cM positions first load the chromosome and bp coordinates into R.

    pop_loci <- read.table("/Users/john.horne/Documents/projects/Dc_linkage_map/map_v2/LinkNe/pop_bp_chr.txt", header=F)
    
    > head(pop_loci)
      V1       V2
    1  1 11028227
    2  1 23416158
    3  1 23416385
    4  1 24898306
    5  1 28115838
    6  1 38259617

### Then you can feed that data into the models, chromosome by chromosome

    ch01_pop_cM <- ch01_bp2cM_predict(pop_loci$V2[which(pop_loci$V1==1)])
    [1] 23.05080 25.00000 25.00000 25.00000 25.00000 27.48837

##### These cM are then available to be used for other applications.
##### But sometimes if SNPs have a bp position that is outside the model range, then it returns an NA value. That can happen if the SNP is near one of the extreme ends of the chromosome. In which case, just delete those loci.
