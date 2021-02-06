# Analyzing Hi-C and HiChIP data with HiCDCPlus

Merve Sahin

02/05/2021

 A necessary task in the analysis of HiC or HiChIP count data is the detection of statistically significant and differential genomic interactions. 
  The count data are available as a table which reports, with regions typically as genomic regions binned uniformly or across restriction enzyme fragments, the number of interactions between pairs of genomic regions. The package HiCDCPlus
  provides methods to determine significant and differential chromatin interactions by use of a
  negative binomial generalized linear model, as well as implementations for TopDom to call topologically associating domains (TADs), and Juicer eigenvector to find the A/B compartments. This vignette explains the use of
  the package and demonstrates typical workflows on HiC and HiChIP data.
  HiCDCPlus package version: 0.99.4
output:
  html_document:
    keep_md: true





**Note:** if you use HiCDCPlus in published research, please cite:

> Sahin, M., Wong, W., Zhan, Y., Van Deyze, K., Koche, R., and Leslie, C. S.
> (2020)
> HiC-DC+: systematic 3D interaction calls and differential analysis 
> for Hi-C and HiChIP
> *BioRxiv*, **335273**.
> [10.1101/2020.10.11.335273](http://dx.doi.org/10.1101/2020.10.11.335273)

# Installation

To install this package, start R and enter:


```r
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("HiCDCPlus")
```

# Standard workflow

## Overview {#overview}

`HiCDCPlus` can take the outputs of the popular Hi-C pre-processing tools such as .hic (from Juicebox), .matrix, and .allValidPairs (from HiC-Pro). It can also be used with HTClist objects (from Bioconductor package HiTC). 

In the standard workflow, one first needs to generate genomic features present in the `HiCDCPlus` model (GC content, mappability, effective length) using the `construct_features` function (see [Creating genomic feature files](#bintolen)). This can be either done for uniformly or multiple restriction fragment binned data. 

`HiCDCPlus` stores counts and features in a memory-efficient way using what we name as a `gi_list` instance(see [The `gi_list` instance](#gi_list)). One next feeds the genomic features in the form of a `gi_list` instance using `generate_bintolen_gi_list` function. Then, counts can be added to this `gi_list` instance using dedicated functions for each input Hi-C file format (`add_hic_counts`, `add_hicpro_matrix_counts`,`add_hicpro_allvalidpairs.counts`).  

Before modeling, 1D features from the `gi_list` instance coming from the bintolen file must be expanded to 2D using `expand_1D_features` function. Different transformations can be applied to combine genomic features derived for each anchor. 

At the core of `HiCDCPlus` is an efficient implementation of the [HiC-DC](https://www.nature.com/articles/ncomms15454) negative binomial count model for normalization and removal of biases
(see ?HiCDCPlus). A platform-agnostic parallelizable implementation is also available in the `HiCDCPlus_parallel` function for efficient interaction 
calling across chromosomes. The `HiCDCPlus` (or `HiCDCPlus_parallel`)
function outputs the significance of each interaction (`pvalue` and FDR
adjusted p-value `qvalue`) along with following estimated from the model: 
1. `mu`: expected interaction frequency estimated from biases, 
2. `sdev`: the standard deviation of expected interaction frequencies. 

Once results are obtained, they can be output into text files using `gi_list_write` function or to a `.hic` file using the `hicdc2hic`function (where one can pass either raw counts, observed/expected normalized
counts, -log10 _P_-value, -log10 _P_-adjusted value, or 
negative binomial Z-score normalized counts: (counts-mu)/sdev to the `.hic` file

To detect differential significant interactions across conditions, `HiCDCPlus` also provides a modified implementation of
[DESeq2](https://bioconductor.org/packages/DESeq2/) using replicate Hi-C/HiChIP datasets `hicdcdiff`. This function requires a
(1) definition of the experimental setup (see ?hicdcdiff), (2) a filtered set of interactions to consider, as a text file containing columns `chr`, `startI`, and `startJ` (startI<=startJ) and (3)
count data per each condition and replicate either as `gi_list` instances or as output text files generated using the `gi_list_write` function that can be read as valid `gi_list` instances using `gi_list_read`.
The `hicdcdiff`
function performs the differential analysis and outputs genomic coordinates of
pairs of regions with corresponding logFC difference, _P_-value and BH adjusted
_P_-value (see the example in [Quickstart](#diff_int)).

We next demonstrate the standard workflow to detect significant as well as differential interactions. 

## Quickstart {#quickstart}

In this section we show a complete workflow for identifying significant
interactions and differential interactions from Hi-C data across replicate
experiments. For HiChIP, the functions used are the same, but the distance thresholds used are slightly reduced (recommended Dmax = 1.5e6).

### Finding Significant Interactions from Hi-C/HiChIP

Here, we identify significant interactions from HiC data at 50kb resolution across multiple chromosomes (in the example
below, across chromosomes 21 and 22). The following example code chunk assumes that
you have downloaded a `.hic` file from [GSE63525](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE63525) and 
also downloaded [Juicebox command line tools](https://github.com/aidenlab/juicer/wiki/Download). Example
below runs with [GSE63525_HMEC_combined.hic](http://hgdownload.soe.ucsc.edu/gbdb/hg19/bbi/hic/GSE63525_HMEC_combined.hic) and stores the path to it into the variable `hicfile_path` with features generated for restriction enzyme fragments with the pattern `"GATC"` in hg19 genome.


```r
hicfile_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic", package = "HiCDCPlus")
outdir<-tempdir(check=TRUE)
#generate features
construct_features(output_path=paste0(outdir,"/hg19_50kb_GATC"),
                   gen="Hsapiens",gen_ver="hg19",
                   sig="GATC",
                   bin_type="Bins-uniform",
                   binsize=50000,
                   chrs=c("chr21","chr22"))
```

```
## [1] "Using chr21 chr22and cut patterns GATC"
## [1] "chr21"
## [1] "chr22"
```

```
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/hg19_50kb_GATC_bintolen.txt.gz"
```

If you have a multiple enzyme cocktail used to generate Hi-C data, you can specify multiple patterns including `"N"` as string to this function (e.g., sig=c("GATC","GANTC")).
If you want to analyze data binned by multiple restriction enzyme fragments, you can change bin_type to "Bins-RE-sites", and binsize to the number of fragments that you would like to merge as bin (e.g., bin_type="Bins-RE-sites" and binsize=10 means 10 restriction fragment binning).


```r
#generate gi_list instance
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(outdir,"/hg19_50kb_GATC_bintolen.txt.gz"))
#add .hic counts
gi_list<-add_hic_counts(gi_list,hic_path = hicfile_path)
```

```
## [1] "Chromosome chr21 intrachromosomal counts processed."
## [1] "Chromosome chr22 intrachromosomal counts processed."
```

If you have HiC-Pro outputs instead, you can use either `add_hicpro_matrix_counts` or `add_hicpro_allvalidpairs.counts` depending on the file format. `add_hicpro_matrix_counts` function requires .bed output from HiC-Pro matrix generation step, together with count data in .matrix format.


```r
#expand features for modeling
gi_list<-expand_1D_features(gi_list)
#run HiC-DC+ on 2 cores
set.seed(1010) #HiC-DC downsamples rows for modeling
gi_list<-HiCDCPlus_parallel(gi_list,ncore=2)
head(gi_list)
```

```
## $chr21
## GInteractions object with 18026 interactions and 9 metadata columns:
##           seqnames1           ranges1     seqnames2           ranges2 |
##               <Rle>         <IRanges>         <Rle>         <IRanges> |
##       [1]     chr21   9450000-9500000 ---     chr21  9950000-10000000 |
##       [2]     chr21   9450000-9500000 ---     chr21 10100000-10150000 |
##       [3]     chr21   9450000-9500000 ---     chr21 10150000-10200000 |
##       [4]     chr21   9450000-9500000 ---     chr21 11000000-11050000 |
##       [5]     chr21   9450000-9500000 ---     chr21 11100000-11150000 |
##       ...       ...               ... ...       ...               ... .
##   [18022]     chr21 47900000-47950000 ---     chr21 47900000-47950000 |
##   [18023]     chr21 47900000-47950000 ---     chr21 47950000-48000000 |
##   [18024]     chr21 47900000-47950000 ---     chr21 48000000-48050000 |
##   [18025]     chr21 47900000-47950000 ---     chr21 48050000-48100000 |
##   [18026]     chr21 48000000-48050000 ---     chr21 48000000-48050000 |
##                   D    counts        gc       map       len        mu      sdev
##           <integer> <numeric> <numeric> <numeric> <numeric> <numeric> <numeric>
##       [1]    500000         0 -0.415451         0 0.2861372  36.62682  20.85646
##       [2]    650000         0 -0.767405         0 0.2689031  29.08634  16.74247
##       [3]    700000         0 -0.461107         0 0.2347835  25.74290  14.91746
##       [4]   1550000         0 -0.795873         0 0.0742543   7.95680   5.17270
##       [5]   1650000         0  0.299047         0 0.2230742   8.34009   5.38471
##       ...       ...       ...       ...       ...       ...       ...       ...
##   [18022]         0      2200  0.633416         0 0.4511124  1789.121  975.8642
##   [18023]     50000       545  0.752876         0 0.1295933   537.613  293.8783
##   [18024]    100000       260  0.725338         0 0.3141704   276.445  151.5580
##   [18025]    150000       271  0.570906         0 0.0480836   125.598   69.3535
##   [18026]         0      2709  0.817260         0 0.1772284  1484.744  809.9996
##              pvalue    qvalue
##           <numeric> <numeric>
##       [1]         1         1
##       [2]         1         1
##       [3]         1         1
##       [4]         1         1
##       [5]         1         1
##       ...       ...       ...
##   [18022] 0.2837514  0.857618
##   [18023] 0.4184339  0.922653
##   [18024] 0.4722949  0.942038
##   [18025] 0.0386625  0.387830
##   [18026] 0.0811332  0.562719
##   -------
##   regions: 963 ranges and 3 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $chr22
## GInteractions object with 18241 interactions and 9 metadata columns:
##           seqnames1           ranges1     seqnames2           ranges2 |
##               <Rle>         <IRanges>         <Rle>         <IRanges> |
##       [1]     chr22 16050000-16100000 ---     chr22 17650000-17700000 |
##       [2]     chr22 16050000-16100000 ---     chr22 17700000-17750000 |
##       [3]     chr22 16050000-16100000 ---     chr22 17800000-17850000 |
##       [4]     chr22 16050000-16100000 ---     chr22 17950000-18000000 |
##       [5]     chr22 16050000-16100000 ---     chr22 18000000-18050000 |
##       ...       ...               ... ...       ...               ... .
##   [18237]     chr22 51000000-51050000 ---     chr22 51100000-51150000 |
##   [18238]     chr22 51000000-51050000 ---     chr22 51150000-51200000 |
##   [18239]     chr22 51050000-51100000 ---     chr22 51050000-51100000 |
##   [18240]     chr22 51050000-51100000 ---     chr22 51150000-51200000 |
##   [18241]     chr22 51150000-51200000 ---     chr22 51150000-51200000 |
##                   D    counts         gc       map       len        mu
##           <integer> <numeric>  <numeric> <numeric> <numeric> <numeric>
##       [1]   1600000         1  0.0731079         0 0.1565584   6.52084
##       [2]   1650000         3  0.0445727         0 0.0809363   5.78804
##       [3]   1750000         3 -0.2745115         0 0.0120789   4.57491
##       [4]   1900000         0 -0.6320482         0 0.0357296   3.31760
##       [5]   1950000         1 -0.2308528         0 0.1088195   3.23818
##       ...       ...       ...        ...       ...       ...       ...
##   [18237]    100000       477   1.350818         0 0.1560164   414.965
##   [18238]    150000       269   1.113756         0 0.3827365   246.312
##   [18239]         0      3903   0.440179         0 0.1683074  2610.342
##   [18240]    100000       350   0.540858         0 0.0930645   361.689
##   [18241]         0      2665   0.641536         0 0.0178216  2415.709
##                sdev    pvalue    qvalue
##           <numeric> <numeric> <numeric>
##       [1]   4.31807  0.974889         1
##       [2]   3.91674  0.795263         1
##       [3]   3.24700  0.704750         1
##       [4]   2.54088  1.000000         1
##       [5]   2.49563  0.899118         1
##       ...       ...       ...       ...
##   [18237]   222.523  0.329603  1.000000
##   [18238]   132.462  0.366718  1.000000
##   [18239]  1394.840  0.163675  0.950785
##   [18240]   194.074  0.453822  1.000000
##   [18241]  1290.907  0.357891  1.000000
##   -------
##   regions: 1027 ranges and 3 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

```r
#write normalized counts (observed/expected) to a .hic file
hicdc2hic(gi_list,hicfile=paste0(outdir,'/GSE63525_HMEC_combined_result.hic'),
          mode='normcounts',gen_ver='hg19')
```

```
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/GSE63525_HMEC_combined_result.hic"
```

```r
#write results to a text file
gi_list_write(gi_list,fname=paste0(outdir,'/GSE63525_HMEC_combined_result.txt.gz'))
```

```
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/GSE63525_HMEC_combined_result.txt.gz"
```

`HiCDCPlus` results can be converted into .hic using `hicdc2hic` function. Values that should be supplied as "mode" into the `hicdc2hic` function for the corresponding score stored in the .hic file are: 'pvalue' for -log10 significance p-value, 'qvalue' for -log10 FDR corrected p-value, 'normcounts' for raw counts/expected counts, 'zvalue' for standardized counts (raw counts-expected counts)/modeled standard deviation of expected counts and 'raw' to pass-through raw counts. 

.hic files can be further converted into .cool format using hic2cool software and be visualized using HiCExplorer. 

## Finding Differential Interactions {#diff_int}

Suppose we're interested in finding differential interactions on `chr21` 
and `chr22` at 50kb between
NSD2 and NTKO/TKO cells given the following `.hic` files available in
[GSE131651](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE131651):
`GSE131651_NSD2_LOW_arima.hic`, `GSE131651_NSD2_HIGH_arima.hic`,
`GSE131651_TKOCTCF_new.hic`, `GSE131651_NTKOCTCF_new.hic`. We first find 
significant interactions in each, and save results to a file:


```r
#generate features
construct_features(output_path=paste0(outdir,"/hg38_50kb_GATC"),
                   gen="Hsapiens",gen_ver="hg38",
                   sig="GATC",bin_type="Bins-uniform",
                   binsize=50000,
                   chrs=c("chr21","chr22"))
```

```
## [1] "Using chr21 chr22and cut patterns GATC"
## [1] "chr21"
## [1] "chr22"
```

```
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/hg38_50kb_GATC_bintolen.txt.gz"
```

```r
#add .hic counts
hicfile_paths<-c(
system.file("extdata", "GSE131651_NSD2_LOW_arima_example.hic", package = "HiCDCPlus"),
system.file("extdata", "GSE131651_NSD2_HIGH_arima_example.hic", package = "HiCDCPlus"),
system.file("extdata", "GSE131651_TKOCTCF_new_example.hic", package = "HiCDCPlus"),
system.file("extdata", "GSE131651_NTKOCTCF_new_example.hic", package = "HiCDCPlus"))
indexfile<-data.frame()
for(hicfile_path in hicfile_paths){
output_path<-paste0(outdir,'/',
                    gsub("^(.*[\\/])", "",gsub('.hic','.txt.gz',hicfile_path)))
#generate gi_list instance
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(outdir,"/hg38_50kb_GATC_bintolen.txt.gz"),
  gen="Hsapiens",gen_ver="hg38")
gi_list<-add_hic_counts(gi_list,hic_path = hicfile_path)
#expand features for modeling
gi_list<-expand_1D_features(gi_list)
#run HiC-DC+ on 2 cores
set.seed(1010) #HiC-DC downsamples rows for modeling
gi_list<-HiCDCPlus(gi_list,ssize=0.1)
for (i in seq(length(gi_list))){
indexfile<-unique(rbind(indexfile,
  as.data.frame(gi_list[[i]][gi_list[[i]]$qvalue<=0.05])[c('seqnames1',
                                                           'start1','start2')]))
}
#write results to a text file
gi_list_write(gi_list,fname=output_path)
}
```

```
## Warning in add_2D_features(gi_list[[chrom]], count_matrix): Bins and counts mismatch. This will slow down the
##             performance of counts integration.Check if genome and/or bin size of
##             counts data is aligned with the GenomicInteractions object.
```

```
## [1] "Chromosome chr21 intrachromosomal counts processed."
```

```
## Warning in add_2D_features(gi_list[[chrom]], count_matrix): Bins and counts mismatch. This will slow down the
##             performance of counts integration.Check if genome and/or bin size of
##             counts data is aligned with the GenomicInteractions object.
```

```
## [1] "Chromosome chr22 intrachromosomal counts processed."
## [1] "Chromosome chr21 complete."
## [1] "Chromosome chr22 complete."
## [1] "Chromosome chr21 intrachromosomal counts processed."
## [1] "Chromosome chr22 intrachromosomal counts processed."
## [1] "Chromosome chr21 complete."
## [1] "Chromosome chr22 complete."
## [1] "Chromosome chr21 intrachromosomal counts processed."
## [1] "Chromosome chr22 intrachromosomal counts processed."
## [1] "Chromosome chr21 complete."
## [1] "Chromosome chr22 complete."
## [1] "Chromosome chr21 intrachromosomal counts processed."
## [1] "Chromosome chr22 intrachromosomal counts processed."
## [1] "Chromosome chr21 complete."
## [1] "Chromosome chr22 complete."
```

```r
#save index file---union of significants at 50kb
colnames(indexfile)<-c('chr','startI','startJ')
data.table::fwrite(indexfile,
            paste0(outdir,'/GSE131651_analysis_indices.txt.gz'),
            sep='\t',row.names=FALSE,quote=FALSE)
```

We next get the union set of significant interactions and save it as the index file, and then run `hicdcdiff`.


```r
#Differential analysis using modified DESeq2 (see ?hicdcdiff)
hicdcdiff(input_paths=list(NSD2=c(paste0(outdir,'/GSE131651_NSD2_LOW_arima_example.txt.gz'),
                 paste0(outdir,'/GSE131651_NSD2_HIGH_arima_example.txt.gz')),
TKO=c(paste0(outdir,'/GSE131651_TKOCTCF_new_example.txt.gz'),
paste0(outdir,'/GSE131651_NTKOCTCF_new_example.txt.gz'))),
filter_file=paste0(outdir,'/GSE131651_analysis_indices.txt.gz'),
output_path=paste0(outdir,'/diff_analysis_example/'),
fitType = 'mean',
binsize=50000,
diagnostics=TRUE)
```

```
## $deseq2paths
## NULL
## 
## $outputpaths
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/diff_resTKOoverNSD2_chr21.txt.gz"
## [2] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/diff_resTKOoverNSD2_chr22.txt.gz"
## 
## $plotpaths
##  [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/sizefactors_chr21.pdf"        
##  [2] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/geomean_sizefactors_chr21.pdf"
##  [3] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/plotMA_TKOoverNSD2_chr21.pdf" 
##  [4] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/diff_chr21_PCA.pdf"           
##  [5] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/dispersionplot.pdf"           
##  [6] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/sizefactors_chr22.pdf"        
##  [7] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/geomean_sizefactors_chr22.pdf"
##  [8] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/plotMA_TKOoverNSD2_chr22.pdf" 
##  [9] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/diff_chr22_PCA.pdf"           
## [10] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/diff_analysis_example/dispersionplot.pdf"
```

```r
#Check the generated plots as well as DESeq2 results
```

Suppose you provide multiple conditions in input_paths such as input_paths=list(A="..",B="..",C=".."), then the pairwise comparisons reported by `hicdcdiff` will be B over A, C over B, C over A.

### ICE normalization using HiTC {#ice}
To find TADs, we use ICE normalized Hi-C data. If you use HiC-Pro to process counts, we suggest feeding ICE normalized .matrix files into a `gi_list` instance. 


```r
gi_list<-generate_binned_gi_list(50000,chrs=c("chr21","chr22"))
gi_list<-add_hicpro_matrix.counts(gi_list,absfile_path,matrixfile_path,chrs=c("chr21","chr22")) #add paths to iced absfile and matrix files here
```

If you have .hic file instead, then you can perform ICE normalization with our HiTC wrapper as follows:


```r
hic_path<-system.file("extdata", "GSE63525_HMEC_combined_example.hic", package = "HiCDCPlus")
gi_list=hic2icenorm_gi_list(hic_path,binsize=50e3,chrs=c('chr21','chr22'))
```

You can also output a ICE normalized .hic file from `hic2icenorm_gi_list` function if you set file_out=TRUE with the name of gsub(".hic","_icenorm.hic",hic_path).

### Finding TADs using TopDom {#topdom}

`HiCDCPlus` converts the gi_list instance with ICE normalized counts into TAD annotations through an implementation of TopDom v0.0.2 (https://github.com/HenrikBengtsson/TopDom) adapted as TopDom. We recommend call TADs with to ICE normalized counts at 50kb resolution with window.size 10 in TopDom. 


```r
tads<-gi_list_topdom(gi_list,chrs=c("chr21","chr22"),window.size = 10)
```

```
## [1] "#########################################################################"
## [1] "Step 0 : File Read "
## [1] "#########################################################################"
## [1] "-- Done!"
## [1] "Step 0 : Done !!"
## [1] "#########################################################################"
## [1] "Step 1 : Generating binSignals by computing bin-level contact frequencies"
## [1] "#########################################################################"
## [1] "Step 1 Running Time :  0.00400000000000489"
## [1] "Step 1 : Done !!"
## [1] "#########################################################################"
## [1] "Step 2 : Detect TD boundaries based on binSignals"
## [1] "#########################################################################"
## [1] "Process Regions from  1 to 665"
## [1] "Step 2 Running Time :  0.00499999999999545"
## [1] "Step 2 : Done !!"
## [1] "#########################################################################"
## [1] "Step 3 : Statistical Filtering of false positive TD boundaries"
## [1] "#########################################################################"
## [1] "-- Matrix Scaling...."
## [1] "-- Compute p-values by Wilcox Ranksum Test"
## [1] "Process Regions from  1 to 665"
## [1] "-- Done!"
## [1] "-- Filtering False Positives"
## [1] "-- Done!"
## [1] "Step 3 Running Time :  0.298000000000002"
## [1] "Step 3 : Done!"
## [1] "Done!!"
## [1] "Job Complete !"
## [1] "#########################################################################"
## [1] "Step 0 : File Read "
## [1] "#########################################################################"
## [1] "-- Done!"
## [1] "Step 0 : Done !!"
## [1] "#########################################################################"
## [1] "Step 1 : Generating binSignals by computing bin-level contact frequencies"
## [1] "#########################################################################"
## [1] "Step 1 Running Time :  0.00400000000000489"
## [1] "Step 1 : Done !!"
## [1] "#########################################################################"
## [1] "Step 2 : Detect TD boundaries based on binSignals"
## [1] "#########################################################################"
## [1] "Process Regions from  1 to 675"
## [1] "Step 2 Running Time :  0.00499999999999545"
## [1] "Step 2 : Done !!"
## [1] "#########################################################################"
## [1] "Step 3 : Statistical Filtering of false positive TD boundaries"
## [1] "#########################################################################"
## [1] "-- Matrix Scaling...."
## [1] "-- Compute p-values by Wilcox Ranksum Test"
## [1] "Process Regions from  1 to 675"
## [1] "-- Done!"
## [1] "-- Filtering False Positives"
## [1] "-- Done!"
## [1] "Step 3 Running Time :  0.385000000000005"
## [1] "Step 3 : Done!"
## [1] "Done!!"
## [1] "Job Complete !"
```

### Finding A/B compartment using Juicer {#comp}

`HiCDCPlus` can call Juicer eigenvector function to determine A/B compartments from .hic files. `extract_hic_eigenvectors` generates text files for each chromosome containing chromosome, start, end and compartment score values that may need to be flipped signs for each chromosome. File paths follow gsub('.hic','_<chromosome>_eigenvalues.txt',hicfile).


```r
extract_hic_eigenvectors(
  hicfile=system.file("extdata", "GSE63525_HMEC_combined_example.hic", package = "HiCDCPlus"),
  mode = "KR",
  binsize = 50e3,
  chrs = "chr22",
  gen = "Hsapiens",
  gen_ver = "hg19"
)
```

# Creating genomic feature files {#bintolen}
Genomic features can be generated using the `construct_features` function. 
This function finds all restriction enzyme cutsites of a given genome and genome
version and computes GC content, mappability (if a relevant 
`.bigWig` file is provided) and effective fragment length for
uniform bin or across specified multiples of restriction enzyme cutsites of
given pattern(s).

```r
#generate features
construct_features(output_path=paste0(outdir,"/hg19_50kb_GATC"),
                   gen="Hsapiens",gen_ver="hg19",
                   sig=c("GATC","GANTC"),bin_type="Bins-uniform",
                   binsize=50000,
                   wg_file=NULL, #e.g., 'hg19_wgEncodeCrgMapabilityAlign50mer.bigWig',
                   chrs=c("chr21","chr22"))
```

```
## [1] "Using chr21 chr22and cut patterns GATC GANTC"
## [1] "chr21"
## [1] "chr22"
```

```
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/hg19_50kb_GATC_bintolen.txt.gz"
```

```r
#read and print
bintolen<-data.table::fread(paste0(outdir,"/hg19_50kb_GATC_bintolen.txt.gz"))
tail(bintolen,20)
```

```
##                        bins     gc map   len
##  1: chr22-49850001-49900000 0.4833   0 49875
##  2: chr22-49900001-49950000 0.5126   0 47521
##  3: chr22-49950001-50000000 0.5139   0 46220
##  4: chr22-50000001-50050000 0.5472   0 48270
##  5: chr22-50050001-50100000 0.5241   0 49289
##  6: chr22-50100001-50150000 0.5014   0 49584
##  7: chr22-50150001-50200000 0.5466   0 48171
##  8: chr22-50200001-50250000 0.5232   0 47970
##  9: chr22-50250001-50300000 0.4675   0 49242
## 10: chr22-50300001-50350000 0.6117   0 41993
## 11: chr22-50350001-50400000 0.1997   0 49875
## 12: chr22-50400001-50450000 0.3623   0 47683
## 13: chr22-50450001-50500000 0.5526   0 49544
## 14: chr22-50500001-50550000 0.5126   0 49105
## 15: chr22-50550001-50600000 0.4790   0 49875
## 16: chr22-50600001-50650000 0.6103   0 49700
## 17: chr22-50650001-50700000 0.5927   0 49531
## 18: chr22-50700001-50750000 0.6473   0 49469
## 19: chr22-50750001-50800000 0.5409   0 49669
## 20: chr22-50800001-50818468 0.4792   0  7806
```

# The `gi_list` instance {#gi_list}

`HiCDCPlus` stores features and count data in a list of `InteractionSet` objects generated for each chromosome, what we name as a `gi_list` instance. 

A `gi_list` instance can be initialized through multiple ways. One can generate
a uniformly binsized `gi_list` instance using `generate_binned_gi_list`. One can
also generate a restriction enzyme fragment binning of the 
genome as a `data.frame` and ingest it as a `gi_list` instance (see 
?generate_df_gi_list) Third, one can generate
some genomic features (GC content, mappability, effective length) and
restriction enzyme fragment regions into as a `bintolen` file
(see [Creating bintolen files](#bintolen)) and generate a `gi_list` instance
from this `bintolen` file. Finally, one can read a `gi_list` instance from a
file generated by `gi_list_write` (see ?gi_list_read). 

## Uniformly binned `gi_list` instance {#uniform}
One can generate a uniform binsized `gi_list` instance for a genome using
`generate_binned_gi_list`:

```r
gi_list<-generate_binned_gi_list(binsize=50000,chrs=c('chr20','chr21'),
                                 gen="Hsapiens",gen_ver="hg19")
head(gi_list)
```

```
## $chr20
## GInteractions object with 52029 interactions and 1 metadata column:
##           seqnames1           ranges1     seqnames2           ranges2 |
##               <Rle>         <IRanges>         <Rle>         <IRanges> |
##       [1]     chr20           0-50000 ---     chr20           0-50000 |
##       [2]     chr20           0-50000 ---     chr20      50000-100000 |
##       [3]     chr20           0-50000 ---     chr20     100000-150000 |
##       [4]     chr20           0-50000 ---     chr20     150000-200000 |
##       [5]     chr20           0-50000 ---     chr20     200000-250000 |
##       ...       ...               ... ...       ...               ... .
##   [52025]     chr20 64300000-64350000 ---     chr20 64350000-64400000 |
##   [52026]     chr20 64300000-64350000 ---     chr20 64400000-64444167 |
##   [52027]     chr20 64350000-64400000 ---     chr20 64350000-64400000 |
##   [52028]     chr20 64350000-64400000 ---     chr20 64400000-64444167 |
##   [52029]     chr20 64400000-64444167 ---     chr20 64400000-64444167 |
##                   D
##           <integer>
##       [1]         0
##       [2]     50000
##       [3]    100000
##       [4]    150000
##       [5]    200000
##       ...       ...
##   [52025]     50000
##   [52026]     97083
##   [52027]         0
##   [52028]     47083
##   [52029]         0
##   -------
##   regions: 1289 ranges and 0 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $chr21
## GInteractions object with 37515 interactions and 1 metadata column:
##           seqnames1           ranges1     seqnames2           ranges2 |
##               <Rle>         <IRanges>         <Rle>         <IRanges> |
##       [1]     chr21           0-50000 ---     chr21           0-50000 |
##       [2]     chr21           0-50000 ---     chr21      50000-100000 |
##       [3]     chr21           0-50000 ---     chr21     100000-150000 |
##       [4]     chr21           0-50000 ---     chr21     150000-200000 |
##       [5]     chr21           0-50000 ---     chr21     200000-250000 |
##       ...       ...               ... ...       ...               ... .
##   [37511]     chr21 46600000-46650000 ---     chr21 46650000-46700000 |
##   [37512]     chr21 46600000-46650000 ---     chr21 46700000-46709983 |
##   [37513]     chr21 46650000-46700000 ---     chr21 46650000-46700000 |
##   [37514]     chr21 46650000-46700000 ---     chr21 46700000-46709983 |
##   [37515]     chr21 46700000-46709983 ---     chr21 46700000-46709983 |
##                   D
##           <integer>
##       [1]         0
##       [2]     50000
##       [3]    100000
##       [4]    150000
##       [5]    200000
##       ...       ...
##   [37511]     50000
##   [37512]     79991
##   [37513]         0
##   [37514]     29991
##   [37515]         0
##   -------
##   regions: 935 ranges and 0 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Restriction enzyme binned `gi_list` instance {#re_sites}
One can also generate an restriction enzyme fragment binning
(indeed, any arbitrary binning) of the 
genome containing columns named `chr` and `start` as a `data.frame` 
(e.g., a `data.frame` read from a `.bed` file) and use it
to generate a `gi_list` instance using `generate_df_gi_list`.


```r
df<-data.frame(chr='chr9',start=c(1,300,7867,103938))
gi_list<-generate_df_gi_list(df)
gi_list
```

```
## $chr9
## GInteractions object with 7 interactions and 1 metadata column:
##       seqnames1          ranges1     seqnames2          ranges2 |         D
##           <Rle>        <IRanges>         <Rle>        <IRanges> | <integer>
##   [1]      chr9            1-300 ---      chr9            1-300 |         0
##   [2]      chr9            1-300 ---      chr9         300-7867 |      3933
##   [3]      chr9            1-300 ---      chr9      7867-103938 |     55752
##   [4]      chr9         300-7867 ---      chr9         300-7867 |         0
##   [5]      chr9         300-7867 ---      chr9      7867-103938 |     51819
##   [6]      chr9      7867-103938 ---      chr9      7867-103938 |         0
##   [7]      chr9 103938-138394717 ---      chr9 103938-138394717 |         0
##   -------
##   regions: 4 ranges and 0 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

## Generating `gi_list` instance from a bintolen file
One can generate genomic features (gc, mappability, effective length)
and restriction enzyme
fragment regions as a `bintolen` file
(see [Creating bintolen files](#bintolen)) first and then generate 
a `gi_list` instance
from it. This instance will readily store genomic features of the
`bintolen` file.

```r
#generate features
construct_features(output_path=paste0(outdir,"/hg19_50kb_GATC"),
                   gen="Hsapiens",gen_ver="hg19",
                   sig="GATC",bin_type="Bins-uniform",
                   binsize=50000,
            wg_file=NULL, #e.g., 'hg19_wgEncodeCrgMapabilityAlign50mer.bigWig',
                   chrs=c("chr21","chr22"))
```

```
## [1] "Using chr21 chr22and cut patterns GATC"
## [1] "chr21"
## [1] "chr22"
```

```
## [1] "/var/folders/31/gz4s1kx132l4_hvgmfc1jl2m0000gn/T//RtmpUrJ86Z/hg19_50kb_GATC_bintolen.txt.gz"
```

```r
#generate gi_list instance
gi_list<-generate_bintolen_gi_list(
  bintolen_path=paste0(outdir,"/hg19_50kb_GATC_bintolen.txt.gz"))
head(gi_list)
```

```
## $chr21
## GInteractions object with 37515 interactions and 1 metadata column:
##           seqnames1           ranges1     seqnames2           ranges2 |
##               <Rle>         <IRanges>         <Rle>         <IRanges> |
##       [1]     chr21           0-50000 ---     chr21           0-50000 |
##       [2]     chr21           0-50000 ---     chr21      50000-100000 |
##       [3]     chr21           0-50000 ---     chr21     100000-150000 |
##       [4]     chr21           0-50000 ---     chr21     150000-200000 |
##       [5]     chr21           0-50000 ---     chr21     200000-250000 |
##       ...       ...               ... ...       ...               ... .
##   [37511]     chr21 46600000-46650000 ---     chr21 46650000-46700000 |
##   [37512]     chr21 46600000-46650000 ---     chr21 46700000-46709983 |
##   [37513]     chr21 46650000-46700000 ---     chr21 46650000-46700000 |
##   [37514]     chr21 46650000-46700000 ---     chr21 46700000-46709983 |
##   [37515]     chr21 46700000-46709983 ---     chr21 46700000-46709983 |
##                   D
##           <integer>
##       [1]         0
##       [2]     50000
##       [3]    100000
##       [4]    150000
##       [5]    200000
##       ...       ...
##   [37511]     50000
##   [37512]     79991
##   [37513]         0
##   [37514]     29991
##   [37515]         0
##   -------
##   regions: 935 ranges and 3 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
## 
## $chr22
## GInteractions object with 40877 interactions and 1 metadata column:
##           seqnames1           ranges1     seqnames2           ranges2 |
##               <Rle>         <IRanges>         <Rle>         <IRanges> |
##       [1]     chr22           0-50000 ---     chr22           0-50000 |
##       [2]     chr22           0-50000 ---     chr22      50000-100000 |
##       [3]     chr22           0-50000 ---     chr22     100000-150000 |
##       [4]     chr22           0-50000 ---     chr22     150000-200000 |
##       [5]     chr22           0-50000 ---     chr22     200000-250000 |
##       ...       ...               ... ...       ...               ... .
##   [40873]     chr22 50700000-50750000 ---     chr22 50750000-50800000 |
##   [40874]     chr22 50700000-50750000 ---     chr22 50800000-50818468 |
##   [40875]     chr22 50750000-50800000 ---     chr22 50750000-50800000 |
##   [40876]     chr22 50750000-50800000 ---     chr22 50800000-50818468 |
##   [40877]     chr22 50800000-50818468 ---     chr22 50800000-50818468 |
##                   D
##           <integer>
##       [1]         0
##       [2]     50000
##       [3]    100000
##       [4]    150000
##       [5]    200000
##       ...       ...
##   [40873]     50000
##   [40874]     84234
##   [40875]         0
##   [40876]     34234
##   [40877]         0
##   -------
##   regions: 1017 ranges and 3 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

# Using custom features with HiCDCPlus

HiCDCPlus allows modeling with user-defined 1D (genomic features for each bin) and 2D (features belonging to an interaction) features. 

Once a `gi_list` instance is at hand, one can ingest counts (and 2D features) using a sparse matrix format text file  containing `chr`, `startI`, `startJ` and `<featurename>` 
columns (see ?add_2D_features) for features you would like to add. `counts`
can be ingested this way as well provided you have a text file containing 
columns named `chr`, `startI` and `startJ`.


```r
df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6))
gi_list<-generate_df_gi_list(df,Dthreshold=500e3,chrs="chr9")
feats<-data.frame(chr='chr9',
startI=seq(1e6,10e6,1e6),
startJ=seq(1e6,10e6,1e6),
counts=rpois(20,lambda=5))
gi_list[['chr9']]<-add_2D_features(gi_list[['chr9']],feats)
gi_list
```

```
## $chr9
## GInteractions object with 10 interactions and 2 metadata columns:
##        seqnames1            ranges1     seqnames2            ranges2 |
##            <Rle>          <IRanges>         <Rle>          <IRanges> |
##    [1]      chr9    1000000-2000000 ---      chr9    1000000-2000000 |
##    [2]      chr9    2000000-3000000 ---      chr9    2000000-3000000 |
##    [3]      chr9    3000000-4000000 ---      chr9    3000000-4000000 |
##    [4]      chr9    4000000-5000000 ---      chr9    4000000-5000000 |
##    [5]      chr9    5000000-6000000 ---      chr9    5000000-6000000 |
##    [6]      chr9    6000000-7000000 ---      chr9    6000000-7000000 |
##    [7]      chr9    7000000-8000000 ---      chr9    7000000-8000000 |
##    [8]      chr9    8000000-9000000 ---      chr9    8000000-9000000 |
##    [9]      chr9   9000000-10000000 ---      chr9   9000000-10000000 |
##   [10]      chr9 10000000-138394717 ---      chr9 10000000-138394717 |
##                D    counts
##        <integer> <numeric>
##    [1]         0         9
##    [2]         0        10
##    [3]         0         4
##    [4]         0         8
##    [5]         0         4
##    [6]         0        12
##    [7]         0         9
##    [8]         0         6
##    [9]         0        10
##   [10]         0        10
##   -------
##   regions: 10 ranges and 0 metadata columns
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

One can also ingest 1D features using a sparse matrix format text file 
containing `chr`, `start` and `<featurename>` (see ?add_1D_features)
and broadcast 1D features to 2D for modeling using a user-specified function
(see ?expand_1D_features). Ingesting 1D features first and then expanding has
a better memory footprint compared to using `add_2D_features` directly.


```r
df<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6),end=seq(2e6,11e6,1e6))
gi_list<-generate_df_gi_list(df)
feats<-data.frame(chr='chr9',start=seq(1e6,10e6,1e6),gc=runif(10))
gi_list<-add_1D_features(gi_list,feats)
gi_list
```

```
## $chr9
## GInteractions object with 27 interactions and 1 metadata column:
##        seqnames1           ranges1     seqnames2           ranges2 |         D
##            <Rle>         <IRanges>         <Rle>         <IRanges> | <integer>
##    [1]      chr9   1000000-2000000 ---      chr9   1000000-2000000 |         0
##    [2]      chr9   1000000-2000000 ---      chr9   2000000-3000000 |   1000000
##    [3]      chr9   1000000-2000000 ---      chr9   3000000-4000000 |   2000000
##    [4]      chr9   2000000-3000000 ---      chr9   2000000-3000000 |         0
##    [5]      chr9   2000000-3000000 ---      chr9   3000000-4000000 |   1000000
##    ...       ...               ... ...       ...               ... .       ...
##   [23]      chr9   8000000-9000000 ---      chr9  9000000-10000000 |   1000000
##   [24]      chr9   8000000-9000000 ---      chr9 10000000-11000000 |   2000000
##   [25]      chr9  9000000-10000000 ---      chr9  9000000-10000000 |         0
##   [26]      chr9  9000000-10000000 ---      chr9 10000000-11000000 |   1000000
##   [27]      chr9 10000000-11000000 ---      chr9 10000000-11000000 |         0
##   -------
##   regions: 10 ranges and 1 metadata column
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```


```r
mcols(InteractionSet::regions(gi_list[['chr9']]))
```

```
## DataFrame with 10 rows and 1 column
##           gc
##    <numeric>
## 1   0.396247
## 2   0.478465
## 3   0.393194
## 4   0.568527
## 5   0.232949
## 6   0.316447
## 7   0.331840
## 8   0.799278
## 9   0.778489
## 10  0.280763
```


```r
gi_list<-expand_1D_features(gi_list)
gi_list
```

```
## $chr9
## GInteractions object with 27 interactions and 2 metadata columns:
##        seqnames1           ranges1     seqnames2           ranges2 |         D
##            <Rle>         <IRanges>         <Rle>         <IRanges> | <integer>
##    [1]      chr9   1000000-2000000 ---      chr9   1000000-2000000 |         0
##    [2]      chr9   1000000-2000000 ---      chr9   2000000-3000000 |   1000000
##    [3]      chr9   1000000-2000000 ---      chr9   3000000-4000000 |   2000000
##    [4]      chr9   2000000-3000000 ---      chr9   2000000-3000000 |         0
##    [5]      chr9   2000000-3000000 ---      chr9   3000000-4000000 |   1000000
##    ...       ...               ... ...       ...               ... .       ...
##   [23]      chr9   8000000-9000000 ---      chr9  9000000-10000000 |   1000000
##   [24]      chr9   8000000-9000000 ---      chr9 10000000-11000000 |   2000000
##   [25]      chr9  9000000-10000000 ---      chr9  9000000-10000000 |         0
##   [26]      chr9  9000000-10000000 ---      chr9 10000000-11000000 |   1000000
##   [27]      chr9 10000000-11000000 ---      chr9 10000000-11000000 |         0
##                gc
##         <numeric>
##    [1] -0.2104882
##    [2]  0.0814864
##    [3] -0.2224653
##    [4]  0.3734610
##    [5]  0.0695092
##    ...        ...
##   [23]   1.921865
##   [24]   0.342575
##   [25]   1.881054
##   [26]   0.301764
##   [27]  -1.277526
##   -------
##   regions: 10 ranges and 1 metadata column
##   seqinfo: 1 sequence from an unspecified genome; no seqlengths
```

# How to get help for HiCDCPlus

Any and all HiCDCPlus questions should be posted to the 
**Bioconductor support site**, which serves as a searchable knowledge
base of questions and answers:

<https://support.bioconductor.org>

Posting a question and tagging with "HiCDCPlus" or "HiC-DC+" will automatically
send an alert to the package authors to respond on the support site.  
You should **not** email your question to the package authors directly, 
as we will just reply that the question should be posted to the 
**Bioconductor support site** instead.

# Session info


```r
sessionInfo()
```

```
## R version 4.0.2 (2020-06-22)
## Platform: x86_64-apple-darwin17.0 (64-bit)
## Running under: macOS Catalina 10.15.7
## 
## Matrix products: default
## BLAS:   /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRblas.dylib
## LAPACK: /Library/Frameworks/R.framework/Versions/4.0/Resources/lib/libRlapack.dylib
## 
## locale:
## [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
## 
## attached base packages:
## [1] stats4    parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] BSgenome.Hsapiens.UCSC.hg38_1.4.3 BSgenome.Hsapiens.UCSC.hg19_1.4.3
##  [3] BSgenome_1.56.0                   rtracklayer_1.48.0               
##  [5] Biostrings_2.56.0                 XVector_0.28.0                   
##  [7] GenomicRanges_1.40.0              GenomeInfoDb_1.24.2              
##  [9] IRanges_2.22.2                    S4Vectors_0.26.1                 
## [11] BiocGenerics_0.34.0               HiCDCPlus_0.2.0                  
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.4-1            ellipsis_0.3.1             
##   [3] biovizBase_1.36.0           htmlTable_2.1.0            
##   [5] base64enc_0.1-3             dichromat_2.0-0            
##   [7] rstudioapi_0.11             farver_2.0.3               
##   [9] bit64_4.0.5                 AnnotationDbi_1.50.3       
##  [11] splines_4.0.2               R.methodsS3_1.8.1          
##  [13] geneplotter_1.66.0          knitr_1.30                 
##  [15] Formula_1.2-3               Rsamtools_2.4.0            
##  [17] annotate_1.66.0             cluster_2.1.0              
##  [19] dbplyr_1.4.4                HiTC_1.32.0                
##  [21] png_0.1-7                   R.oo_1.24.0                
##  [23] compiler_4.0.2              httr_1.4.2                 
##  [25] backports_1.1.10            assertthat_0.2.1           
##  [27] Matrix_1.2-18               lazyeval_0.2.2             
##  [29] htmltools_0.5.0             prettyunits_1.1.1          
##  [31] tools_4.0.2                 igraph_1.2.5               
##  [33] gtable_0.3.0                glue_1.4.2                 
##  [35] GenomeInfoDbData_1.2.3      dplyr_1.0.2                
##  [37] rappdirs_0.3.1              Rcpp_1.0.5                 
##  [39] Biobase_2.48.0              vctrs_0.3.4                
##  [41] xfun_0.18                   stringr_1.4.0              
##  [43] lifecycle_0.2.0             ensembldb_2.12.1           
##  [45] XML_3.99-0.5                InteractionSet_1.16.0      
##  [47] zlibbioc_1.34.0             MASS_7.3-53                
##  [49] scales_1.1.1                VariantAnnotation_1.34.0   
##  [51] GenomicInteractions_1.22.0  hms_0.5.3                  
##  [53] ProtGenerics_1.20.0         SummarizedExperiment_1.18.2
##  [55] AnnotationFilter_1.12.0     RColorBrewer_1.1-2         
##  [57] yaml_2.2.1                  curl_4.3                   
##  [59] memoise_1.1.0               gridExtra_2.3              
##  [61] ggplot2_3.3.2               biomaRt_2.44.1             
##  [63] rpart_4.1-15                latticeExtra_0.6-29        
##  [65] stringi_1.5.3               RSQLite_2.2.1              
##  [67] genefilter_1.70.0           checkmate_2.0.0            
##  [69] GenomicFeatures_1.40.1      BiocParallel_1.22.0        
##  [71] rlang_0.4.7                 pkgconfig_2.0.3            
##  [73] matrixStats_0.57.0          bitops_1.0-6               
##  [75] evaluate_0.14               lattice_0.20-41            
##  [77] purrr_0.3.4                 labeling_0.3               
##  [79] GenomicAlignments_1.24.0    htmlwidgets_1.5.2          
##  [81] bit_4.0.4                   tidyselect_1.1.0           
##  [83] magrittr_1.5                DESeq2_1.28.1              
##  [85] R6_2.4.1                    generics_0.0.2             
##  [87] Hmisc_4.4-1                 DelayedArray_0.14.1        
##  [89] DBI_1.1.0                   pillar_1.4.6               
##  [91] foreign_0.8-80              survival_3.2-7             
##  [93] RCurl_1.98-1.2              nnet_7.3-14                
##  [95] tibble_3.0.3                crayon_1.3.4               
##  [97] BiocFileCache_1.12.1        rmarkdown_2.4              
##  [99] jpeg_0.1-8.1                progress_1.2.2             
## [101] locfit_1.5-9.4              grid_4.0.2                 
## [103] data.table_1.13.0           blob_1.2.1                 
## [105] digest_0.6.25               xtable_1.8-4               
## [107] tidyr_1.1.2                 R.utils_2.10.1             
## [109] openssl_1.4.3               munsell_0.5.0              
## [111] Gviz_1.32.0                 askpass_1.1
```
