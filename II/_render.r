#Render

#'
#+ library and wd
library(rmarkdown)
Sys.setenv(RSTUDIO_PANDOC = "C:/Users/katoo/AppData/Local/Pandoc")
options(repo = "https://cran.rstudio.com/")

#' from Rmd to pdf
rmarkdown::render(
  input = "II/TAsession_4/handout.rmd",
  output_file = "handout.pdf",
  output_dir = "II/TAsession_4", #いじらない
  clean = TRUE, #いじらない
  encoding = "utf8" #いじらない
)

#' from Rmd to R
knitr::purl(
    input = "II/TAsession_4/handout.rmd",
    output = "R/script/TASession_4.r",  
    documentation = 1
)