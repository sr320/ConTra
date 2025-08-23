# Install required R packages for the ConTra Shiny app

# Set user library path
.libPaths(c(Sys.getenv("R_LIBS_USER"), .libPaths()))

# Install packages if they are not already installed
install_if_missing <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, dependencies = TRUE, repos = "https://cran.rstudio.com/", 
                     lib = Sys.getenv("R_LIBS_USER"))
    library(pkg, character.only = TRUE)
  }
}

# Core packages
install_if_missing("shiny")
install_if_missing("shinydashboard") 
install_if_missing("DT")
install_if_missing("plotly")

# Data manipulation
install_if_missing("readr")
install_if_missing("dplyr")
install_if_missing("tidyr")
install_if_missing("stringr")

# Plotting
install_if_missing("ggplot2")
install_if_missing("scales")
install_if_missing("RColorBrewer")

# Additional utility packages
install_if_missing("shinyWidgets")

cat("All packages installed successfully!\n")