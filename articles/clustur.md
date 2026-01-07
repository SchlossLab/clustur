# clustur

## Background

clustur was developed to be similar to [mothur’s cluster
function](https://mothur.org/wiki/cluster/) that was written in C++. In
order to cluster your data, users need to provide their own sparse or
phylip-formatted distance matrix. They also need to provide a count
table that either comes from mothur or that they create in R. Once these
objects are built users can call the
[`cluster()`](http://www.schlosslab.org/clustur/reference/cluster.md)
function. We currently support 5 methods: opticlust (default) and
furthest, nearest, weighted, and average neighbor. The opticlust method
is [`cluster()`](http://www.schlosslab.org/clustur/reference/cluster.md)
and mothur’s default. The speed of the methods implemented in {clustur}
and mothur are comparable; {clustur} may even be faster! Below we will
show you how to create your sparse matrix and count table. If you do not
have a count table, clustur can produce one from you, but it will assume
the abundance of each sequence is one and it will only cluster the
sequences in the distance matrix. The output of running
[`clustur()`](http://www.schlosslab.org/clustur/reference/clustur-package.md)
includes what is typically provided in a mothur-formatted shared file.

## Starting Up

For the official release from CRAN you can use the standard
[`install.packages()`](https://rdrr.io/r/utils/install.packages.html)
function:

``` r
# install via cran
install.packages("clustur")
```

For the developmental version, you can use the `install_github()`
function from the {devtools} package

``` r
# install via github
devtools::install_github("SchlossLab/clustur")
```

Because {clustur}’s functions make use of a random number generator,
users are strongly encouraged to set the seed.

``` r
library(clustur)
set.seed(19760620)
```

## Read count files

clustur will produce the same output using either a sparse (default) or
full [count table](https://mothur.org/wiki/count_file/)

``` r
full_count_table <- read_count(example_path("amazon.full.count_table"))
sparse_count_table <- read_count(example_path("amazon.sparse.count_table"))
```

## Read distance matrix file

clustur will read both mothur’s [column/sparse distance
matrix](https://mothur.org/wiki/column-formatted_distance_matrix/) and
[Phylip-formatted
distance](https://mothur.org/wiki/phylip-formatted_distance_matrix/)
matrix formats.

``` r
column_distance <- read_dist(example_path("amazon_column.dist"), full_count_table, cutoff = 0.03)
```

or

``` r
phylip_distance <- read_dist(example_path("amazon_phylip.dist"), full_count_table, cutoff = 0.03)
```

The return value of `distance_data` will be a memory address. If you
want a data frame version of the distances, you can use
`get_distance_df(distance_data)`.

``` r
get_distance_df(column_distance)
#>    FirstName SecondName Distance
#> 1     U68638     U68618 0.020396
#> 2     U68638     U68620 0.020396
#> 3     U68638     U68658 0.027067
#> 4     U68618     U68620 0.000000
#> 5     U68618     U68658 0.022512
#> 6     U68620     U68658 0.022512
#> 7     U68641     U68667 0.000000
#> 8     U68641     U68673 0.018238
#> 9     U68667     U68673 0.018238
#> 10    U68636     U68631 0.006024
#> 11    U68680     U68615 0.003141
#> 12    U68679     U68663 0.020354
#> 13    U68679     U68665 0.017144
#> 14    U68679     U68688 0.009987
#> 15    U68665     U68663 0.012295
#> 16    U68665     U68688 0.008059
#> 17    U68663     U68688 0.020272
get_distance_df(phylip_distance)
#>    FirstName SecondName Distance
#> 1     U68615     U68680 0.003141
#> 2     U68618     U68620 0.000000
#> 3     U68618     U68638 0.020396
#> 4     U68618     U68658 0.022512
#> 5     U68620     U68638 0.020396
#> 6     U68620     U68658 0.022512
#> 7     U68631     U68636 0.006024
#> 8     U68638     U68658 0.027067
#> 9     U68641     U68667 0.000000
#> 10    U68641     U68673 0.018238
#> 11    U68663     U68665 0.012295
#> 12    U68663     U68679 0.020354
#> 13    U68663     U68688 0.020272
#> 14    U68665     U68679 0.017144
#> 15    U68665     U68688 0.008059
#> 16    U68667     U68673 0.018238
#> 17    U68679     U68688 0.009987
```

## Clustering the data

The default method for clustering in `cluster` is “opticlust”

``` r
cutoff <- 0.03
cluster_data <- cluster(column_distance, cutoff)
```

### Selecting different clustering methods

``` r

cluster_data <- cluster(column_distance, cutoff, method = "furthest")
cluster_data <- cluster(column_distance, cutoff, method = "nearest")
cluster_data <- cluster(column_distance, cutoff, method = "average")
cluster_data <- cluster(column_distance, cutoff, method = "weighted")
```

## Output data from clustering

All methods produce a list object with an indicator of the cutoff that
was used (`label`), as well as cluster composition (`cluster`) and
shared (`abundance`) data frames. The `clusters` data frame shows which
OTU (Operation Taxonomic Unit) each sequence was assigned to. The
`abundance` data frame contains columns indicating the `OTU` and
`sample` identifiers and the abundance of each OTU in each sample. The
OptiClust method also includes the `metrics` data frame, which describe
the optimization value for each iteration in the fitting process; the
data in `clusters` and `shared` are taken from the last iteration.
clustur provides getter functions, `get_label()`, `get_clusters()`,
`get_shared()`, and
[`get_metrics()`](http://www.schlosslab.org/clustur/reference/get_metrics.md),
which will be demonstrated below.

``` r
clusters <- cluster(column_distance, cutoff, method = "opticlust")
get_cutoff(clusters)
#> [1] 0.03
get_bins(clusters)
#>   feature  bin
#> 1  U68630 bin1
#> 2  U68595 bin2
#> 3  U68600 bin3
#> 4  U68591 bin4
#> 5  U68647 bin5
#> 6  U68661 bin6
#> 7  U68605 bin7
#>  [ reached 'max' / getOption("max.print") -- omitted 81 rows ]
get_abundance(clusters)
#>   samples  bin abundance
#> 1 pasture bin1         0
#> 2  forest bin1         1
#> 3 pasture bin2         0
#> 4  forest bin2         1
#> 5 pasture bin3         0
#>  [ reached 'max' / getOption("max.print") -- omitted 171 rows ]
get_metrics(clusters)
#> $metrics
#>      label   cutoff specificity      ppv       ttp  f1score          tn
#> 1 0.030000 0.030000    1.000000 1.000000 17.000000 1.000000 4736.000000
#>        mcc       fn       fp sensitivity      npv      fdr accuracy
#> 1 1.000000 0.000000 0.000000    1.000000 1.000000 1.000000 1.000000
#> 
#> $iteration_metrics
#>      iter time label num_otus cutoff tp tn fp fn sensitivity specificity ppv
#>      npv fdr accuracy
#>  [ reached 'max' / getOption("max.print") -- omitted 2 columns ]
#>  [ reached 'max' / getOption("max.print") -- omitted 3 rows ]
```
