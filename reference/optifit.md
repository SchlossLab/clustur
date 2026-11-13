# Optifit

Optifit is the 'cluster.fit' command in mothur. It allows you to fit
sequences to OTUs or existing OTUS. Just like mothur, there are two ways
of running this command: open, or closed. Closed fits reads into
existing OTUs, and scraps any reads that cannot be fitted(default it
closed). Open fits reads into existing OTUs while unfitted OTUs are
clustered into new OTUs.

## Usage

``` r
optifit(
  fit_dist,
  cutoff,
  feature_column_name = "feature",
  bin_column_name = "bin",
  accnos = NULL,
  ref_fasta = NULL,
  ref_count = NULL,
  ref_list = NULL,
  fit_fasta = NULL,
  fit_percent = 50,
  closed = FALSE,
  print_ref = TRUE,
  random_seed = 123,
  max_iterations = 100,
  denovo_iters = 10,
  optimal_denovo_result = TRUE,
  ref_weight = "none",
  metric = "mcc",
  pairwise_distance_method = "onegap",
  number_of_threads = 1,
  strollur_object = NULL
)
```

## Arguments

- fit_dist:

  The distance object that will be used to fit to the new reference.
  This must be created using the \`read_dist()\` function.

- cutoff:

  The cutoff you want to cluster towards.

- feature_column_name:

  Set the name of the column in the cluster dataframe that contains the
  sequence names. This is a \`character\` value.

- bin_column_name:

  Set the name of the column in the cluster dataframe that contains the
  name of the group of sequence names. This is a \`character\` value.

- accnos:

  a \`vector\` of strings that represent the reference sequences by
  name. Use this when you know the group of sequences you want your fit
  data to fit against.

- ref_fasta:

  a\`data.frame\` representation of your fasta data for your reference
  sequences. This can be generated using the \`read_fasta()\` function.
  The fasta \`data.frame\` must contain two columns: "sequence_name" and
  "sequence."

- ref_count:

  a \`data.frame\` that display sequence abundance for the reference
  sequences. Must be generated using the \`read_count()\` function.

- ref_list:

  a \`data.frame\` respentation of the mothur list file you wish to fit
  your sequences to. Can be generated using the \`read_mothur_list()\`
  function. Must include two columns: "bin_name" and "sequence_name."

- fit_fasta:

  a\`data.frame\` representation of your fasta data for your fitted
  sequences. This can be generated using the \`read_fasta()\` function.
  The fasta \`data.frame\` must contain two columns: "sequence_name" and
  "sequence."

- fit_percent:

  allows you to set the percentage of reads that need to be fitted. This
  is a \`numerical\` value.

- closed:

  deteremines whether you want the function to run using the closed or
  open method. The closed method will return reads that fit into
  existing otus and scrap those that did not. The open method will keep
  reads that fit into existing otus and create new OTUs for those that
  did not. This is a \`logical\` value.

- print_ref:

  controls whether all the reference sequences are printed in the final
  result. This is a \`character\` value.

- random_seed:

  a \`integer\` value that represents the current random seed.

- max_iterations:

  the max number of iterations that cluster computes. This is a
  \`integer\` value.

- denovo_iters:

  the max number of iterations the denovo method computes. This is a
  \`integer\` value.

- optimal_denovo_result:

  determines whether or not you want the denovo method to return the
  cluster with the highest mcc value or all of the randomly computed
  clusters. This is a \`logical\` value.

- ref_weight:

  Ref Weight is used to randomize the reference sequences for the denovo
  method. There are three methods: none, abundance and connectivity.
  "none" will do a normal shuffle of the given sequences without any
  weights, abundance will shuffle the sequences using abundances as a
  weight, and connectivity will shuffle the sequences but uses a weight
  of one for every sequence.

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

A \`list\` of a list of \`data.frames\` that contain abundance, and
clustering results. If you used \`method = opticlust\`, it will also
return clustering performance metrics.

## Examples

``` r
 amazon_count <- read_count(example_path("amazon.full.count_table"))
 amazon_dist <- read_dist(example_path("amazon_phylip.dist"),
                          amazon_count, 0.03)

 esophagus_fasta <- read_fasta(example_path("esophagus.fasta"))
 esophagus_count <- read_count(example_path("esophagus.count_table"))
 amazon_fasta <- read_fasta(example_path("amazon.align"))
 esophagus_list <- read_mothur_list(example_path("esophagus.list"))


 denovo <- optifit(amazon_dist, 0.03)
#> accnos, ref_fasta, ref_count, ref_list, and fit_fastaparameters are all NULL. Running Optifit using denovo.
 denovo_open <- optifit(amazon_dist, 0.03, close = FALSE)
#> accnos, ref_fasta, ref_count, ref_list, and fit_fastaparameters are all NULL. Running Optifit using denovo.
 accnos <- optifit(amazon_dist, 0.03,
                   accnos = amazon_count$Representative_Sequence[1:40])
#> Accnos file is not null, running optifit using accnos.
 accnos_open <- optifit(amazon_dist, 0.03,
                        accnos = amazon_count$Representative_Sequence[1:40],
                        closed=FALSE)
#> Accnos file is not null, running optifit using accnos.

 ref_cluster <- optifit(amazon_dist, 0.03, "feature", "bin", accnos = NULL,
                        ref_fasta = esophagus_fasta,
                        ref_count = esophagus_count,
                        ref_list = esophagus_list, fit_fasta = amazon_fasta,
                        pairwise_distance_method = "onegap",
                        number_of_threads = parallel::detectCores(),
                        print_ref = TRUE)
#> ref_fasta, ref_list, ref_count, and fit_fasta are not NULL, running reference based optifit.
```
