# One-click deployment script for shinyapps.io

# Usage:
# 1) Install rsconnect once: install.packages("rsconnect")
# 2) Set account info (only once per machine):
#    rsconnect::setAccountInfo(name="YOUR_NAME", token="YOUR_TOKEN", secret="YOUR_SECRET")
# 3) From the repo root or shiny_app dir, run:
#    Rscript shiny_app/deploy.R

if (!requireNamespace("rsconnect", quietly = TRUE)) {
  install.packages("rsconnect")
}

# Resolve app directory robustly whether invoked from repo root or inside shiny_app
args <- commandArgs(trailingOnly = FALSE)
file_arg <- args[grep("^--file=", args)]
if (length(file_arg) > 0) {
  script_path <- normalizePath(sub("^--file=", "", file_arg[1]))
  script_dir <- dirname(script_path)
} else {
  script_dir <- getwd()
}

# If the script lives inside shiny_app, deploy that directory; otherwise try ./shiny_app
if (basename(script_dir) == "shiny_app") {
  app_dir <- script_dir
} else if (dir.exists("shiny_app")) {
  app_dir <- normalizePath("shiny_app")
} else {
  app_dir <- script_dir
}

message("Deploying app from ", app_dir)

# Optional: allow specifying account/server via env vars
account <- Sys.getenv("RSCONNECT_ACCOUNT")
server  <- Sys.getenv("RSCONNECT_SERVER", "shinyapps.io")
appName <- Sys.getenv("RSCONNECT_APPNAME")

if (nzchar(account) && nzchar(appName)) {
  rsconnect::deployApp(appDir = app_dir, account = account, server = server, appName = appName)
} else if (nzchar(account)) {
  rsconnect::deployApp(appDir = app_dir, account = account, server = server)
} else if (nzchar(appName)) {
  rsconnect::deployApp(appDir = app_dir, appName = appName)
} else {
  rsconnect::deployApp(appDir = app_dir)
}


