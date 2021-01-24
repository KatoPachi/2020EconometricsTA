#Render

#'
#+ library and wd
library(rmarkdown)
library(bookdown)
Sys.setenv(RSTUDIO_PANDOC = "C:/Users/katoo/AppData/Local/Pandoc")
options(repo = "https://cran.rstudio.com/")

## ---- render

# pdf (rmarkdown)
rmarkdown::render(
  input = "II/TAsession_8/handout.rmd",
  output_file = "handout.pdf",
  output_dir = "II/TAsession_8", #いじらない
  clean = TRUE, #いじらない
  encoding = "utf8" #いじらない
)

# pdf (bookdown)
bookdown::render_book(
    input = "index.rmd",
    output_format = "bookdown::pdf_book",
    output_dir = ".",
    output_file = "Econometrics2TA.pdf",
    clean = TRUE,
    encoding = "utf8"
)

#' from Rmd to R
knitr::purl(
    input = "II/TAsession_8/handout.rmd",
    output = "R/script/TASession_8.r",  
    documentation = 1
)
