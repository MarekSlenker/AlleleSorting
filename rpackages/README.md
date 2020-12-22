# R packages

Required [R](https://www.r-project.org/) packages `ape`, `seqinr`, `filelock`, and `parallel` and their dependencies for R version used **must** be installed here.

Installation of source packages in Linux requires basic compilation tools and developmental files for GDAL, GEOS and UDUNITS-2.

Within `R` command line started in `AlleleSorting` directory use e.g. command

```R
install.packages(pkgs=c("ape", "seqinr", "filelock", "parallel"),
lib="rpackages", repos="https://mirrors.nic.cz/R/", dependencies="Imports")
```

to install needed packages.
