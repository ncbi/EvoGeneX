# EvoGeneX

Ornstein-Uhlenbeck models for Phylogenetic Comparative Hypotheses for gene
expression evolution that utilizes within-species variation.

#### Source install

```
devtools::install_github("ncbi/EvoGeneX", subdir="Rpackage", build_vignettes=TRUE)
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

#### Citation

Soumitra Pal, Brian C. Oliver, and Teresa M. Przytycka. Modeling gene expression evolution with EvoGeneX uncovers differences in evolution of species, organs and sexes. Journal of Computational Biology 2022.
