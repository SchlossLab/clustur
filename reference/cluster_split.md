# Cluster Split

Cluster Split splits data into groups based on taxomic group. It then
clusters them using the following methods: opticluster, average
neighbor, nearest neighbor, furthest neighbor, and weighted neighbor.
You can read more about \[cluster
split\](https://mothur.org/wiki/cluster.split/) on the mothur wiki.

## Usage

``` r
cluster_split(
  fasta,
  taxonomy,
  count_table,
  cluster_method,
  cutoff,
  taxonomy_cutoff,
  feature_column_name = "feature",
  bin_column_name = "bin",
  random_seed = 123,
  delta = 1e-04,
  max_iterations = 100,
  initialize = "singleton",
  metric = "mcc",
  pairwise_distance_method = "onegap",
  number_of_threads = 1,
  strollur_object = NULL
)
```

## Arguments

- fasta:

  a\`data.frame\` representation of your fasta data. This can be
  generated using the \`read_fasta\` function. The fasta \`data.frame\`
  must contain two columns: sequence_name and sequence.

- taxonomy:

  a\`data.frame\` representation of your taxonomy data. This can be
  generated using the \`read_mothur_taxonomy\` function. The taxonomy
  \`data.frame\` must contain two columns: sequence_name and taxonomy.

- count_table:

  the count table that represents the data, must be read using the
  \`read_count()\` function.

- cluster_method:

  the clustering method you wish to use. Your clustering options are:
  opticluster, average neighbor, nearest neighbor, furthest neighbor,
  and weighted neighbor.

- cutoff:

  a \`numerical\` value that represents the distance cutoff.

- taxonomy_cutoff:

  a \`integer\` value that represents the taxonomic level cutoff.

- feature_column_name:

  Set the name of the column in the cluster dataframe that contains the
  sequence names. This is a \`character\` value.

- bin_column_name:

  Set the name of the column in the cluster dataframe that contains the
  name of the group of sequence names. This is a \`character\` value.

- random_seed:

  a \`integer\` value that represents the current random seed.

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

- pairwise_distance_method:

  the method you wish to use to calculate pairwise distances between
  aligned sequences. As of current we only allow the calculation of
  fasta sequences, and the current methods are; onegap, nogap, and
  eachgap. For more information look at the (mother
  webpage)\[https://mothur.org/wiki/dist.seqs/\]. This is a
  \`character\` value.

- number_of_threads:

  a \`integer\` value that represents the amount of allotted threads
  that can be used.

- strollur_object:

  Using the \`strollur\` r package, you can use your created strollur
  object to save your clustering results to. This allows you to seamless
  modify your strollur object and analysis results as needed.

## Value

A list of \`data.frames\` that contain abundance, and clustering
results.

## Examples

``` r

 cutoff <- 0.2
 count_table <- read_count(example_path("amazon.full.count_table"))
 fasta <- read_fasta(example_path("amazon.align"))
 taxonomy <- read_mothur_taxonomy(example_path("amazon.taxonomy"))
 df <- cluster_split(fasta, taxonomy, count_table, "opti",
                     cutoff, 3, number_of_threads = 1)
```
