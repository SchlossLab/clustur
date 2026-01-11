# Changelog

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
