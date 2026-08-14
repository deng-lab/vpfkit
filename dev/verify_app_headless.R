# Drive the viewer in a real browser against a multi-sample ViroProfiler object.
#
#   NOT_CRAN=true CHROMOTE_CHROME=/snap/bin/chromium \
#     Rscript dev/verify_app_headless.R
#
# Two ways this check can pass for the wrong reason, both of which it has:
#
#   1. mod_data_input requires pressing "Load dataset" after choosing a file.
#      Uploading alone leaves every tab in its empty state, where nothing can
#      fail. The check asserts the "No dataset loaded" banner is gone first.
#   2. Shiny suspends outputs on inactive tabs. A walk over the top-level navbar
#      never renders anything on a sub-tab, so the assay table and the
#      ordination -- the two most interesting outputs -- stay uncomputed. Each
#      sub-tab is activated explicitly.
#
# Point TSE at an object with at least three samples and a grouping column;
# ordination and PERMANOVA are unreachable otherwise, and the point of this
# script is to exercise exactly those paths.
suppressMessages(pkgload::load_all("/home/allen/github/rujinlong/vpfkit", quiet = TRUE))
library(shinytest2)

TSE <- Sys.getenv(
  "VPFKIT_TEST_TSE",
  "/mnt/nas26/testdata/viroprofiler_16sample/run_final/results/viroprofiler_output.rds"
)
if (!file.exists(TSE)) stop("No test object at ", TSE, "; set VPFKIT_TEST_TSE.")
# A snap-packaged chromium has a private /tmp and cannot see an upload staged there.
staged <- file.path(path.expand("~"), ".vpfkit_headless_upload.rds")
file.copy(TSE, staged, overwrite = TRUE)
on.exit(unlink(staged), add = TRUE)
app <- AppDriver$new(app_dir = ".", name = "vpf16d", height = 1000, width = 1600,
                     load_timeout = 120000, timeout = 60000)
on.exit(try(app$stop(), silent = TRUE), add = TRUE)

app$set_inputs(`data-source` = "upload")
app$upload_file(`data-upload` = staged)
app$click("data-load")
app$wait_for_idle(timeout = 180000)

txt <- function() trimws(gsub("\\s+", " ", gsub("<[^>]+>", " ", app$get_html("body"))))
stopifnot(!grepl("No dataset loaded", txt(), fixed = TRUE))
cat("dataset loaded: 16 samples\n\n")

visit <- function(label, ...) {
  app$set_inputs(..., wait_ = FALSE)
  app$wait_for_idle(timeout = 120000)
  h <- app$get_html("body")
  e <- lengths(regmatches(h, gregexpr("shiny-output-error", h)))
  cat(sprintf("%-38s error-panels=%d\n", label, e))
  e
}

total <- 0
total <- total + visit("Data / Assays",        `data-tabs` = "Assays")
total <- total + visit("Data / Sample table",  `data-tabs` = "Sample table (colData)")
total <- total + visit("Data / Library sizes", `data-tabs` = "Library sizes")

h <- app$get_html("body")
cat("\nassay descriptions now rendered:",
    grepl("Trimmed mean coverage depth", h, fixed = TRUE), "\n")

total <- total + visit("Diversity / Alpha",       main_nav = "Diversity")
total <- total + visit("Diversity / Ordination",  `diversity-tabs` = "Ordination")
h <- app$get_html("body")
cat("PERMANOVA reported on the ordination tab:", grepl("PERMANOVA", h, fixed = TRUE), "\n")

total <- total + visit("Genes",   main_nav = "Genes")
total <- total + visit("Taxonomy", main_nav = "Taxonomy")
total <- total + visit("Contigs",  main_nav = "Contigs")
total <- total + visit("Host & lifestyle", main_nav = "Host & lifestyle")
total <- total + visit("Export",   main_nav = "Export")

logs <- tryCatch(as.data.frame(app$get_logs()), error = function(e) NULL)
be <- if (!is.null(logs) && "level" %in% names(logs)) sum(tolower(logs$level) %in% c("error","severe")) else 0
cat("\nTOTAL error panels:", total, "| browser errors:", be, "\n")
if (total + be > 0) quit(status = 1)
