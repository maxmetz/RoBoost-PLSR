## roboost - RoBoost-PLSR : robust partial least squares regression 

## <span style="color:green"> **Installation** </span> 

Using [**Rstudio**](https://www.rstudio.com/products/rstudio/download/) is recommended for installation and usage

### <span style="color:green"> 1.  Install package **'remotes'** from CRAN </span>

Use the **Rstudio** menu 

or write in the R console
```{r}
install.packages("remotes")
```

### <span style="color:green"> 2. Install package **'roboost'** </span> 

**a) Most recent version**

Write in the R console
```{r}
remotes::install_github("maxmetz/RoBoost-PLSR", dependencies = TRUE, 
  build_vignettes = TRUE)
```

In case of the following question during installation process:
```{r}
These packages have more recent versions available.
Which would you like to update?"
```
it is recommended to skip updates (usually choice **3** = None)

### <span style="color:green"> 3. Usage </span>

Write in the R console
```{r}
library(roboost)
```

## <span style="color:green"> **Author** </span> 

**Maxime Metz**

- [**ChemHouse**](https://www.chemproject.org/ChemHouse), Montpellier

**maxime.metz57200@gmail.com**

### How to cite
