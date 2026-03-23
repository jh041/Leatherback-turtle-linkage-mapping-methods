# Chromosome linkage map instructions

### Here I go over how to estimate cM values for SNP loci when arranged in their base-pair order.

### This process has fewer steps than linkage grouping and uses less of the functionality of Lep-MAP3, so for a more in-depth exploration of the software see the [linkage grouping notes](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/Linkage_grouping/LG_map_instructions.md).

### Data processing follows as described in the linkage grouping notes, and we start here with the data.call.gz file which is the output of the ParentCall2 module in Lep-MAP3.

### First, we are going to thin this down, so that no SNP is within 10 bp of another polymorphism.

    zcat data.call.gz | cut -f 1,2 | awk '(NR>=7)' > snps.txt

    awk 'BEGIN{print "#"}NR>1{if (/NC/ && !($1 in d)) d[$1]=++n; if ($1!=p1 || $2-p>10) {print (d[$1]+0);p=$2} else print 0; p1=$1}' snps.txt > map.txt

##### We are now dealing with 368,409 SNPs and we are going to make an index for them

    seq 1 368409 > phys.txt

    for i in {1..28}; do echo "zcat data.call.gz | java -cp bin OrderMarkers2 numMergeIterations=3 data=- map=map.txt calculateIntervals=orders/order$i.int chromosome=$i evaluateOrder=phys.txt improveOrder=0 proximityScale=100 > orders/order$i.txt 2>orders/order$i.err"; done | parallel --jobs 4

##### As with the linkage grouping process, "bin" is the directory where the java executables are found.

##### And orders/ is a folder where we are sending the output

##### Here, the module is doing something different than in the linkage grouping exercise. We already know the order, and are using the option "evaluateOrder" to 

##### The output of OrderMarkers2 looks like this...

    #marker_number	male_position	female_position	( parental_phase )[	phased data]
    1	0.000	0.000	( AC,AC;--,-- )	1011101001000101 0 0	0110111111111101 0 0
    2	0.000	0.000	( --,CA;CA,-- )	1011101001000101 0 0	0110111111111101 0 0
    3	0.000	0.000	( CA,--;AC,-- )	1011101001000101 0 0	0110111111111101 0 0

##### But we want the SNP positions so we can make marey maps.

### Add the SNP positions

    for i in {1..28}; do awk -vn=$i '/^[^#]/{print $1,n,$2,$3}' orders/order$i.txt; done|awk -vOFS="\t" 'NR==FNR{m[NR-1]=$0}NR!=FNR{$1=m[$1];print}' snps.txt -> map.txt

##### Now we have this...

    NC_050095.1	6698553	28	56.250	68.750
    NC_050095.1	6698575	28	56.250	68.750
    NC_050095.1	6698836	28	56.250	68.750
    NC_050095.1	6698849	28	56.250	68.750

##### All chromosomes are in a single file, with the: chromosome (28), bp position, linkage group (28), and male and female cM values respectively.

##### Click here to view [marey maps](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/Chromosome_linkage_mapping/marey.png) for each of the 28 chromosomes.

### For the paper, this map was chosen over the linkage grouping map because it retains more markers and eliminates noise from background linkage that we later learned was a strong signal in this genome. With a good chromosome-level genome assembly this method is straightforward and relatively easy. However, without the high-quality reference genome, it is far better to let the linkage groups determine marker order.

### Let's look at some stats for this linkage map in R

    LM_stats <- function(M, Ch) {
    Chromosome <- M[which(M$V3==Ch),]
    bp <- max(Chromosome$V2)
    male <- max(Chromosome$V4)
    female <- max(Chromosome$V5)
    sex_av <- (male + female)/2
    sa_rate <- bp / sex_av
    return(c(Ch, format(bp, scientific=F), format(male, scientific=F), format(female, scientific=F), format(sex_av, scientific=F), format(sa_rate, scientific=F)))
    }
    
    LM_stats(map, 1)
    [1] "1"         "354445513" "125"       "131.25"    "128.125"   "2766404"  
    [1] "2"         "272701501" "131.25"    "187.5"     "159.375"   "1711068"  
    [1] "3"         "212153107" "143.75"    "200"       "171.875"   "1234345"  
    [1] "4"         "146519289" "100"       "150"       "125"       "1172154"  
    [1] "5"         "137567980" "81.25"     "137.5"     "109.375"   "1257764"  
    [1] "6"         "130998374" "118.75"    "181.25"    "150"       "873322.5" 
    [1] "7"         "127639529" "87.5"      "175"       "131.25"    "972491.6" 
    [1] "8"         "109268961" "81.25"     "125"       "103.125"   "1059578"  
    [1] "9"         "105169864" "62.5"      "181.25"    "121.875"   "862932.2" 
    [1] "10"       "86472350" "87.5"     "156.25"   "121.875"  "709516.7"
    [1] "11"       "79988231" "87.5"     "137.5"    "112.5"    "711006.5"
    [1] "12"       "44227750" "62.5"     "125"      "93.75"    "471762.7"
    [1] "13"       "41137519" "75"       "93.75"    "84.375"   "487555.8"
    [1] "14"       "40018190" "56.25"    "100"      "78.125"   "512232.8"
    [1] "15"       "33183484" "56.25"    "93.75"    "75"       "442446.5"
    [1] "16"       "26399810" "56.25"    "75"       "65.625"   "402282.8"
    [1] "17"       "25426300" "37.5"     "75"       "56.25"    "452023.1"
    [1] "18"       "23658492" "31.25"    "62.5"     "46.875"   "504714.5"
    [1] "19"       "20011669" "31.25"    "43.75"    "37.5"     "533644.5"
    [1] "20"       "19189351" "68.75"    "87.5"     "78.125"   "245623.7"
    [1] "21"       "18857373" "68.75"    "81.25"    "75"       "251431.6"
    [1] "22"       "18744624" "43.75"    "68.75"    "56.25"    "333237.8"
    [1] "23"       "17220591" "75"       "56.25"    "65.625"   "262409"  
    [1] "24"       "16917659" "56.25"    "68.75"    "62.5"     "270682.5"
    [1] "25"       "16451672" "62.5"     "93.75"    "78.125"   "210581.4"
    [1] "26"       "16414810" "43.75"    "75"       "59.375"   "276460"  
    [1] "27"       "16292564" "56.25"    "62.5"     "59.375"   "274401.1"
    [1] "28"       "6698849"  "56.25"    "68.75"    "62.5"     "107181.6"

    bp <-             c(354445513,272701501,212153107,146519289,137567980,130998374,127639529,109268961,105169864,86472350,79988231,44227750,41137519,40018190,33183484,26399810,25426300,23658492,20011669,19189351,18857373,18744624,17220591,16917659,16451672,16414810,16292564,6698849)
    sum(bp)
    [1] 2163775406

    male_cM <- c(125,131.25,143.75,100,81.25,118.75,87.5,81.25,62.5,87.5,87.5,62.5,75,56.2,56.2,56.2,37.5,31.25,31.25,68.75,68.75,43.75,75,56.25,62.5,43.75,56.25,56.25)
    sum(male_cM)
    [1] 2043.6

    female_cM <- c(131.25,187.5,200,150,137.5,181.25,175,125,181.25,156.25,137.5,125,93.75,100,93.75,75,75,62.5,43.75,87.5,81.25,68.75,56.25,68.75,93.75,75,62.5,68.75)
    sum(female_cM)
    [1] 3093.75

    rowMeans(cbind(male_cM, female_cM))
    > sex_av_cM <-     c(128.125,159.375,171.875,125,109.375,150,131.250,103.125,121.875,121.875,112.5,93.75,84.375,78.1,74.975,65.6,56.25,46.875,37.5,78.125,75,56.25,65.625,62.5,78.125,59.375,59.375,62.5)
    sum(sex_av_cM)
    [1] 2568.675

    2163775406/2568.675
    [1] 842370.3          ### 1 cM for every 842 kbp

  
