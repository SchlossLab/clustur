# clustur: Clustering

A tool that implements the clustering algorithms from 'mothur' (Schloss
PD et al. (2009)
[doi:10.1128/AEM.01541-09](https://doi.org/10.1128/AEM.01541-09) ).
'clustur' make use of the cluster() and make.shared() command from
'mothur'. Our cluster() function has five different algorithms
implemented: 'OptiClust', 'furthest', 'nearest', 'average', and
'weighted'. 'OptiClust' is an optimized clustering method for
Operational Taxonomic Units, and you can learn more here, (Westcott SL,
Schloss PD (2017)
[doi:10.1128/mspheredirect.00073-17](https://doi.org/10.1128/mspheredirect.00073-17)
). The make.shared() command is always applied at the end of the
clustering command. This functionality allows us to generate and create
clustering and abundance data efficiently.

## See also

Useful links:

- <http://www.schlosslab.org/clustur/>

- <https://github.com/SchlossLab/clustur>

- Report bugs at <https://github.com/SchlossLab/clustur/issues>

## Author

**Maintainer**: Patrick Schloss <pschloss@umich.edu>
([ORCID](https://orcid.org/0000-0002-6935-4275)) \[copyright holder\]

Authors:

- Gregory Johnson <grejoh@umich.edu>
  ([ORCID](https://orcid.org/0009-0008-3890-0297))

- Sarah Westcott <swestcot@umich.edu>
