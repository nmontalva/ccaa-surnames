# ccaa-surnames

## What is this?

This repository is to share and store data files from project 111160402: Cultural phylogenetics and coevolution of wealth inheritance and land tenure norms in agropastoralist communities.


## How to use it (on GNU/linux)

### Prerequisites

1. install [R](https://www.r-project.org/) and/or [RStudio](https://www.rstudio.com/)
1. clone the repo (or do a `git pull` if you cloned it already)
```
git clone https://github.com/nmontalva/ccaa-surnames.git
```
2. go to the `bienes-raices_2017` folder
```
$ cd ccaa-surnames/bienes-raices_2017
```
3. install the dependencies for this project and download [tabula](https://github.com/tabulapdf/tabula-java/releases/download/v1.0.1/tabula-1.0.1-jar-with-dependencies.jar) on the `bienes-raices_2017` folder (you'll need `gcc-fortran` to compile `tidyverse` and `tk` to use `install.packages` in the R console)
```r
$ wget https://github.com/tabulapdf/tabula-java/releases/download/v1.0.1/tabula-1.0.1-jar-with-dependencies.jar # you can download it from the browser too
$ R
> install.packages("tidyverse")
> install.packages("rvest")
> install.packages("ghit")
> ghit::install_github(c("ropensci/tabulizerjars", "ropensci/tabulizer"))
> install.packages("Biodem")
> install.packages("ggdendro")
> install.packages("reldist")
```

### How to re-generate file `commoners.csv`

1. go to the `bienes-raices_2017` folder
4. run the first R script
```r
$ R
> source('download-pdfs.R')
> download_pdfs()
```
5. this downloads all the most recent pdf files from the [OTCA](http://www.comunidadesagricolas.cl/) in a folder `pdfs`
6. run the second R script
```r
> source('extract-pdfs.R')
> df <- extract_pdfs(read_pdfs())
```
7. all information from the pdf files is stored in data frame `df` and saved in csv file `commoners.csv`

### Run the analysis

```r
> source('analysis.R')
> dendrogram(df)
```


## Abstract

Cultural norms of wealth inheritance and property rights are regarded as important mechanisms in the emergence of inequality and social conflict. Cross-cultural analyses and case-studies in the anthropological literature on the subject suggest a relationship between subsistence patterns, ecological and economic conditions faced by the societies and the characteristics of wealth and property exhibited by different groups. The main aim of this project is to better understand this relationship in the Agricultural Communities of the Coquimbo region, considered here as an advantageous study model due to their agropastoralist subsistence pattern, their system of communal land tenure, and the impact of unstable economic and ecological conditions in their livelihood.

However, the comparative study of variation in wealth inheritance and land tenure norms between communities poses the challenge of controlling for shared history as a source of spurious associations between these traits and the ecological and economic conditions proposed as their determinants. Phylogenetic comparative method has been successfully used in anthropological studies to deal with this problem, known in the literature as \textquotedblleft Galton's Problem\textquotedblright, but in most instances has been aimed to the study of large variation between several linguistic groups.

Methodologically, this project will study the relationship between wealth inheritance and property rights using phylogenetic comparative method at a small-scale. During the first year, an initial exploratory approach will use the list of registered commoners and community rights obtained from the Land registry authority to build phylogenetic trees using distance matrices based on differences in surnames distributions between communities. The same lists of surnames will be used to estimate inheritance norms of multigeniture, unigeniture, and sex-bias in intergenerational wealth transmission. This estimations based solely on documented material will be compared with data collected at fieldwork during the second year, where plans are to obtain genetic distances based on microsatellite markers and ethnographic accounts of practices of wealth transmission.

The rationale behind this approach is based on the fact that property rights, surnames and genes are all vertically transmitted from parents to offspring, and thus tracking the ancestor--descendant relationships using one trait can give clues about the genealogical history of one of the others. Trees and traits obtained from both phases will be compared and combined to evaluate which trees show the highest statistical accordance with all the data. It is expected to obtain as a result an scheme of the cultural evolution of these traits (i.e. their changes through the branches of the tree) using phylogenetic methods. It is also expected to use these results to test specific hypotheses about the origin of the Agricultural Communities, and general theoretical questions about the emergence of communal land systems, contributing with a framework to study comparatively the effects of ecological, demographic and economic conditions on wealth and property.

## Geographic maps of frequencies

Sometimes called "interpolation maps" are maps used to predict the frequencie of a trait over a geographic area, deducing likely values of a trait at every point of the map based on a few geographic points where data were actually collected.

...tbc

## Acknowledgment

This porject was funded by CONICYT Fondecyt Iniciacion N°11160402
