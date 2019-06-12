# manuscript folder

This folder contains files needed to generate the manuscript and supplemental material. Both were written in R markdown, and use the `knitr` package integrated into Rstudio to turn the `.Rmd` files into `.tex` files and then convert that to a `.pdf`. You will need to have a Latex distribution to generate the `.pdf` files, as well as a variety of Latex packages. 

### Files

`header.tex` - Latex header to be included in the manuscript. You must have the listed packages to properly compile the manuscript.

`hz.bib` - The bibliography file.

`manuscript.Rmd`- The R markdown notebook of the manuscript. When knit in Rstudio, generates the `manuscript.tex` file and a pdf.

`supplemental.Rmd`- The R markdown notebook of the supplemental material. When knit in Rstudio, generates the `supplemental.tex` file and a pdf.