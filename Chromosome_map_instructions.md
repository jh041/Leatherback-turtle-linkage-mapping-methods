# Chromosome linkage map instructions

### Here I go over how to estimate cM values for SNP loci when arranged in their base-pair order.

### This process has fewer steps than linkage grouping and uses less of the functionality of Lep-MAP3, so for a more in-depth exploration of the software see the [linkage grouping notes](https://github.com/jh041/Leatherback-turtle-linkage-mapping-methods/blob/main/LG_map_instructions.md).

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

> Written with [StackEdit](https://stackedit.io/).

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

  