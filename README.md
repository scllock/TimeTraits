# TimeTraits

`TimeTraits` provides a set of tools for analysing biological time-series data
using functional data analysis (FDA). The package is designed to support
end-to-end workflows, from smoothing rhythmic time series to extracting
group-level traits from functional principal component analysis (FPCA).

The methods implemented here are particularly suited to circadian and other
biological time-series data where interest lies in comparing functional
patterns across experimental groups, time windows, or curve derivatives.

---

## Features

- Smoothing of biological time-series data using functional data representations
- Functional principal component analysis (FPCA) of smoothed curves
- Extraction of group-level FPCA-derived traits
- Support for multiple curve derivatives (e.g. 0th, 1st, 2nd)
- Optional temporal segmentation (e.g. pre/post environmental shifts)
- Shape-based outlier detection (where applicable)

---

## Installation

You can install the development version of `TimeTraits` directly from GitHub:

```r
# install.packages("devtools")
devtools::install_github("scllock/TimeTraits")
