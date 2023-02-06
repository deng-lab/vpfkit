FROM rocker/shiny-verse:4.2

LABEL author="Jinlong Ru"

RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libgsl-dev \
    libclang-dev \
    libglpk-dev

COPY renv* .

# install R packages using renv
RUN R -e 'install.packages(c("renv","remotes"))'
# RUN R -e "install.packages('Matrix', version='1.5-1')"

RUN R -e "renv::restore()"

# COPY DESCRIPTION .
# RUN Rscript -e "remotes::install_deps()"

# RUN R -e 'remotes::install_bioc("scater","mia","miaViz","SummarizedExperiment","SingleCellExperiment","TreeSummarizedExperiment")'
# RUN R -e 'renv::install("bioc::scater","bioc::mia","bioc::miaViz","bioc::SummarizedExperiment","bioc::SingleCellExperiment","bioc::TreeSummarizedExperiment")'


COPY *.tgz /app.tgz
RUN R -e 'remotes::install_local("/app.tgz",upgrade="never")'
RUN rm /app.tgz
# RUN R -e 'remotes::install_github("deng-lab/vpfkit@03b0dbaf02a4cef3ad9609508c0be0d55aab6895")'


EXPOSE 80

CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');vpfkit::run_app()"
