# randomEffects

<!-- badges: start -->
<!-- badges: end -->

The goal of randomEffects is to fit random-intercept logistic regression
models by approximating the marginal likelihood, and to compare the
approximation methods in a simulation study. It was written by Group 5 for
MT4113 Assignment 2.

Four approximation methods are implemented:

- `laplace`: Laplace approximation
- `gh`: Gauss–Hermite quadrature
- `agh`: adaptive Gauss–Hermite quadrature
- `is`: importance sampling

## Installation

You can install the development version of randomEffects from [GitHub](https://github.com/) with:

```r
# install.packages("pak")
pak::pak("fmcewan/randomEffects")
```

## Example

This example fits a random-intercept logistic model using adaptive
Gauss–Hermite quadrature:

```r
library(randomEffects)

cd  <- cluster_data(y, x, group)
fit <- fit_rint(cd, method = "agh", K = 7)
fit$par
```

## Reproducing the report

The report renders to PDF, which needs a LaTeX installation. If you don't have one:

    install.packages("tinytex")
    tinytex::install_tinytex()

(or `quarto install tinytex` from the terminal)

The report needs some extra packages (listed under `Suggests` in
`DESCRIPTION`) and [Quarto](https://quarto.org).

1. Clone the repository and open `randomEffects.Rproj` in RStudio.
2. Install the package and all its dependencies:
   ```r
   # install.packages("pak")
   pak::local_install_dev_deps()
   pak::local_install()
   ```
3. Render the report:
   ```bash
   quarto render Report.qmd
   ```

Seeds are set in the report, so the results are reproducible.

## Repository structure

| Path         | Contents                                   |
|--------------|--------------------------------------------|
| `R/`         | Package functions                          |
| `tests/`     | Unit tests (`testthat`)                    |
| `raw_data/`  | Data used in the analysis                  |
| `results/`   | Saved simulation results                   |
| `Report/`    | Report qmd & rendered pdf                  |

## Authors

Group 5: <Lourens Francois Janse van Rensburg, Helena Rothery , Fraser McEwan, Robin Ferguson, Andrea Totaro>
