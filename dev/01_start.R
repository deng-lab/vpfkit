golem::fill_desc(
  pkg_name = "vpfkit",
  pkg_title = "Analyzing ViroProfiler results",
  pkg_description = "An R package for analyzing ViroProfiler results",
  author_first_name = "Jinlong",
  author_last_name = "Ru",
  author_email = "jinlong.ru@gmail.com",
  repo_url = "https://github.com/deng-lab/vpfkit"
)

golem::set_golem_options()
usethis::use_gpl3_license()
usethis::use_readme_rmd(open = FALSE)
usethis::use_code_of_conduct(contact = "jinlong.ru@gmail.com")
usethis::use_git()

# finishing setup git
golem::use_recommended_tests()
# golem::use_favicon("inst/app/www/favicon.png")
golem::use_utils_ui(with_test = TRUE)
golem::use_utils_server(with_test = TRUE)
rstudioapi::navigateToFile("dev/02_dev.R")
