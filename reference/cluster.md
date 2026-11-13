# Cluster entities together

Clusters entities represented in a distance matrix and count table using
one of several algorithms and outputs information about the composition
and abundance of each cluster

## Usage

``` r
cluster(
  distance_object,
  cutoff,
  method = "opticlust",
  feature_column_name_to = "feature",
  bin_column_name_to = "bin",
  random_seed = 123,
  delta = 1e-04,
  max_iterations = 100,
  initialize = "singleton",
  metric = "mcc",
  strollur_object = NULL
)
```

## Arguments

- distance_object:

  The distance object that was created using the \`read_dist()\`
  function.

- cutoff:

  The cutoff you want to cluster towards.

- method:

  The method of clustering to be performed: opticlust (default),
  furthest, nearest, average, or weighted.

- feature_column_name_to:

  Set the name of the column in the cluster dataframe that contains the
  sequence names.

- bin_column_name_to:

  Set the name of the column in the cluster dataframe that contains the
  name of the group of sequence names.

- random_seed:

  the random seed to use, (default = 123).

- delta:

  Allows you to set the stable value for the metric in the opticluster
  method. It defaults to 0.0001, to reach full convergence set this
  value to 0.

- max_iterations:

  the max amount of iterations the opticluster method will run. The
  default is 100.

- initialize:

  this parameter allows you to set the default randomization of for the
  opticluster method. The two methods are: "singleton" and "oneotu."
  Singleton will randomly assign each sequence to its own OTU while
  "oneotu" assigns all sequence into one OTU. Defaults to singleton.

- metric:

  this parameter allows you to choose which metric the cluster will
  optimize to, the options are: Matthews correlation coefficient (mcc),
  sensitivity (sens), specificity (spec), true positives + true
  negatives (tptn), false positives + false negatives (fpfn), true
  positives (tp), true negative (tn), false positive (fp), false
  negative (fn), f1score (f1score), accuracy (accuracy), positive
  predictive value (ppv), negative predictive value (npv), false
  discovery rate (fdr). The default metric is mcc.

- strollur_object:

  Using the \`strollur\` r package, you can use your created strollur
  object to save your clustering results to. This allows you to seamless
  modify your strollur object and analysis results as needed.

## Value

A list of \`data.frames\` that contain abundance, and clustering
results. If you used \`method = opticlust\`, it will also return
clustering performance metrics.

## Examples

``` r

 cutoff <- 0.03
 count_table <- read_count(example_path("amazon.full.count_table"))
 distance_data <- read_dist(example_path("amazon_column.dist"),
                            count_table, cutoff)

 cluster_results <- cluster(distance_data,
                            cutoff, method = "opticlust",
                            feature_column_name_to = "sequence",
                            bin_column_name_to = "omu")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "furthest")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "nearest")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "average")
 cluster_results <- cluster(distance_data,
                            cutoff, method = "weighted")

```
