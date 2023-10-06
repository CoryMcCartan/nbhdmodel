# nbhdmodel

R package to fit neighborhood models as described in ["Measuring and Modeling Neighborhoods"](https://arxiv.org/pdf/2110.14014).
The core methodology is described in the paper and can be implemented with any tool that can fit generalized linear mixed models (GLMMs).
However, some of the preprocessing necessary to set up the GLMM can be onerous.
In addition to providing a specialized GLMM routine, this package provides several preprocessing functions that, while not completely general, should be useful for others performing these kinds of analyses.

To install, run:

```r
remotes::install_github("CoryMcCartan/nbhdmodel")
```
