

[![R-CMD-check](https://github.com/SchlossLab/Clustur/actions/workflows/check-full.yml/badge.svg)](https://github.com/SchlossLab/Clustur/actions/workflows/check-full.yml)
[![test-coverage](https://github.com/SchlossLab/Clustur/actions/workflows/test-coverage.yml/badge.svg)](https://github.com/SchlossLab/Clustur/actions/workflows/test-coverage.yml)
[![lint](https://github.com/SchlossLab/Clustur/actions/workflows/lintr.yml/badge.svg)](https://github.com/SchlossLab/Clustur/actions/workflows/lintr.yml)
[![pkgdown](https://github.com/SchlossLab/Clustur/actions/workflows/pkgdown.yaml/badge.svg)](https://github.com/SchlossLab/Clustur/actions/workflows/pkgdown.yaml)

# mothur

mothur is an open-sourced software used for the analysis of ecological data. It has a strong presence in the use of analyzing 16S rRNA gene sequences. mothur itself has a plethora of features, such as, clustering based on genetic similarity, sequence alignment, taxomic assignment, alpha and beta diversity, and phylotype analysis. There have been numerous publications using this software by individuals and intuitions alike. We plan to take advantage of the power of mothur, and isolate the clustering feature. We believe this will allow scientist to work with a familiar tool and give rise for the creation of more pipelines utilizing this package.

# clustur

clustur was developed for two main reasons: to allow developers to easily create pipelines with mothur code in R, and to give users a powerful clustering algorithm from mothur without having to download it. Because mothur is built using rcpp, it benefits from the speed of c++ and the usability of R. Additionally, since clustur will be open source, developers can modify are package as needed. This will not only help us continue to create high quality packages, but will allow everyone to take advantage of the foundation that is already built.


## Installation

```r
# install via cran
install.packages("clustur")
```

```r
# install via github
remotes::install_github("SchlossLab/Clustur")
```

## Getting Started

To get started, look at the ["Getting started"]() page.

## Contributions

For feature requrest please visit the issue page in github and use the feature tag. If you would like to contribute pull request are welcomed!
