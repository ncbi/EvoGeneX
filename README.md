# EvoGeneX

Ornstein-Uhlenbeck models for Phylogenetic Comparative Hypotheses for gene
expression evolution that utilizes within-species variation.

#### Source install

```
devtools::install_github("ncbi/EvoGeneX", build_vignettes=TRUE)
```

#### Vignettes

See [vignette](https://rawcdn.githack.com/ncbi/EvoGeneX/56a4e2d5c7c98eb57983d2d67d89fcfb61beb076/Rpackage/vignettes/EvoGeneX.html).

Within an R session:
```
vignette("EvoGeneX", package="EvoGeneX")
```

#### Examples

- [example_constrained.R](examples/example_constrained.R) shows how to determine if a gene is constrained or neutral
- [example_adaptive.R](examples/example_adaptive.R) shows how to determine if a gene is adaptive or not
- [example_multigene.R](examples/example_multigene.R) shows how to determine if a set of genes are individually constrained or neutral

#### Drosophila data to access from R

See [vignette](https://rawcdn.githack.com/ncbi/EvoGeneX/56a4e2d5c7c98eb57983d2d67d89fcfb61beb076/Rpackage/vignettes/EvoGeneX.html#4_Drosophila_datasets).

#### Drosophila data and results in MS Excel format

The following Microsoft Excel files provide the Drosophila gene expression data and the results mentioned in the EvoGeneX paper. The data and results for each of the 5 body-parts and the 2 sexes (total 10) is kept as a separate sheet in each of the Excel files. Additionally, each Excel file includes a ReadMe sheet that gives details of the columns in the rest of the sheets.

* [DrosophilaNormalizedGeneExpression.xlsx](inst/supplementary/DrosophilaNormalizedGeneExpression.xlsx): Normalized gene expression data of 8591 orthologs from 9 drosophila species each with 4 replicates.
* [DrosophilaGenesConstrainedOrNot.xlsx](inst/supplementary/DrosophilaGenesConstrainedOrNot.xlsx): For each of the 8591 orthologs if it was predicted by EvoGeneX to undergo constrained expression evolution as compared to neutral evolution.
* [DrosophilaAdaptiveGenesTwoRegime.xlsx](inst/supplementary/DrosophilaAdaptiveGenesTwoRegime.xlsx): List of orthologs predicted to undergo adaptive evolution under the two-regime scenario.
* [DrosophilaAdaptiveGenesThreeRegime.xlsx](inst/supplementary/DrosophilaAdaptiveGenesThreeRegime.xlsx): List of orthologs predicted to undergo adaptive evolution under the three-regime scenario.

#### Citation

Soumitra Pal, Brian C. Oliver, and Teresa M. Przytycka. Modeling gene expression evolution with EvoGeneX uncovers differences in evolution of species, organs and sexes. Journal of Computational Biology 2023.
