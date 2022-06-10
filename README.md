# ASP: Alternative Splicing Package <img src="https://github.com/mhjiang97/asp/blob/master/data-raw/asp.png" align="right" height="150" width="140/"/>

### *an R package for performing, integrating and visualizing differential alternative splicing analyses from RNA-seq expediently*

## Author

Minghao Jiang, [jiangminghao1001\@163.com](mailto:jiangminghao1001@163.com)

------------------------------------------------------------------------

## Table of Contents

> -   [Dependencies](#Dependencies)  
> -   [Installation](#Installation)  
> -   [Vignette](#Vignette)  
> -   [License](#License)

## Dependencies

-   Differential alternative splicing tools *(follow the original guides and install them first)*
    -   [rMATS](https://github.com/Xinglab/rmats-turbo), please make `rmats.py` executable.
    -   [CASH](https://soft.novelbio.com/cash/), the path to `cash.jar` is needed.
    -   [MAJIQ](https://majiq.biociphers.org), please make `majiq` and `voila` executable.
    -   [LeafCutter](https://davidaknowles.github.io/leafcutter/), the path to the downloaded `leafcutter` directory is needed.
    -   [SplAdder](https://spladder.readthedocs.io/en/latest/), please make `spladder` executable. SplAdder version should be less than 3.0.0 (Result formats of the Latest version is not compatible with the version 2.0).
    -   [BANDITS](https://bioconductor.org/packages/release/bioc/html/BANDITS.html), an R package.
    -   [SUPPA](https://github.com/comprna/SUPPA), please make `suppa.py` executable.
    -   [psichomics](https://bioconductor.org/packages/release/bioc/html/psichomics.html), an R package.
-   [Conda](https://anaconda.org/anaconda/conda)
    -   Different tools can rely on different Conda environments.
    -   Different tools can rely on different Conda or virtual environments.
    -   Besides, standing alone is always favourable.

## Installation

``` r
devtools::install_github("mhjiang97/asp", build_vignettes = TRUE)
```

## Vignette

``` r
browseVignettes("asp")
```

------------------------------------------------------------------------

## License

**ASP** is licensed under the [GNU General Public License v3](http://www.gnu.org/licenses/gpl-3.0.html)
