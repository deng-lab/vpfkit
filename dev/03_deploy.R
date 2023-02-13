devtools::check()
devtools::build(path = getwd(), binary = T)
# golem::add_dockerfile_with_renv(output_dir = "docker", from = "rocker/shiny-verse")
