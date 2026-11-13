# Changelog

## clustur 0.2.0

## clustur 0.2.0

### New Features

- Added the new function
  [`optifit()`](http://www.schlosslab.org/clustur/reference/optifit.md),
  which is the cluster.fit method from mothur. You can read up on it
  [here](https://mothur.org/wiki/cluster.split/).
- Added the new function `cluster_fit()`, which is the cluster.split
  method from mothur. You can read up on it
  [here](https://mothur.org/wiki/cluster.fit/).
- Added support for `strollur` objects from the `strollur` package.

### Improvements

- Optimized all clustering methods to be faster and more memory
  efficient.
- Removed code bloat from the last version.
- Added more customization to all the clustering functions, you can now
  input parameters that were in mothur in clustur (i.e metrics, and
  max_iterations).

## clustur 0.1.4

CRAN release: 2026-01-08

- Updated the
  [`cluster()`](http://www.schlosslab.org/clustur/reference/cluster.md)
  function to change the current clustering metric it is using if it is
  not suitable for your data.

- [`read_dist()`](http://www.schlosslab.org/clustur/reference/read_dist.md)
  and
  [`cluster()`](http://www.schlosslab.org/clustur/reference/cluster.md)
  will now return objects with class name “distance_object”, and
  “mothur_cluster” respectively. This will ensure the returned object
  can be differentiated between other R objects.

- Fixed memory leaks that were found using Valgrind.

## clustur 0.1.3

CRAN release: 2025-04-21

- Fixed
  [`cluster()`](http://www.schlosslab.org/clustur/reference/cluster.md)
  returning duplicate cluster names when using method = “opticlust.”
