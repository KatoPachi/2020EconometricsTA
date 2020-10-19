#Render

#'
#+ library and wd
library(rmarkdown)
Sys.setenv(RSTUDIO_PANDOC = "C:/Users/katoo/AppData/Local/Pandoc")
options(repo = "https://cran.rstudio.com/")

#' from Rmd to pdf
rmarkdown::render(
  input = "II/201030/handout.rmd",
  output_file = "handout.pdf",
  output_dir = "II/201030", #いじらない
  clean = TRUE, #いじらない
  encoding = "utf8" #いじらない
)

#' from Rmd to R
knitr::purl(
    input = "II/201009/handout.rmd",
    output = "R/#1.r",  
    documentation = 1
)