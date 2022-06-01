# EvoGeneX

Ornstein-Uhlenbeck models for Phylogenetic Comparative Hypotheses that utilizes
within-species variation.

#### Source install

```
library(devtools)  
install_github("ncbi/EvoGeneX", subdir="Rpackage")
```

#### Vignettes

Within an R session:
```
vignette("EvoGeneX", package="EvoGeneX")
```
`

#### Examples

- `examples/example_constrained.R` shows how to determine if a gene is constrained or neutral
- `examples/example_adaptive.R` shows how to determine if a gene is adaptive or not
- `examples/example_multigene.R` shows how to determine if a set of genes are individually constrained or neutral
