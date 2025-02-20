## mxfda v0.2.2-1 (development version)

-   Reapplied family parameter to the scalar on function regression (`run_sofr()`)
-   Added limit for number of knots in `run_sofr()`. If the number of knots is not specified or greater than the number of measured radii *and* calculated to be more than 50 knots, limit to only be 50 knots.

## mxfda v0.2.2

-   added vignette for working with spatial transcriptomic data as provided in a seurat object from a clear cell renal cell carcinoma CosMx SMI study. Data can be downloaded [here](https://doi.org/10.5281/zenodo.12730227).
-   from original v0.1.0 submitted to [CRAN](https://cran.r-project.org/web/packages/mxfda/index.html).
-   implemented Vu et al bivariate entropy measure.
-   ther bug fixes. See [Commit History](https://github.com/julia-wrobel/mxfda/commits/main/) for complete list of commits associated with the main branch of mxFDA.
-   added ability to calculate an empirical estimate for complete spatial randomness (CSR) when extracting summary functions

## mxfda v0.1.0

-   starting package out
