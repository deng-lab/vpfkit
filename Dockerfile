FROM rocker/shiny-verse:4.2

LABEL author="Jinlong Ru"

# copy the app to the image
COPY . /srv/shiny-server/

WORKDIR /srv/shiny-server/

RUN apt-get update && apt-get install -y \
    libcurl4-gnutls-dev \
    libssl-dev \
    libgsl-dev \
    libclang-dev

# install R packages using renv
RUN R -e "install.packages('renv')"
# RUN R -e "renv::init()"
# RUN R -e "renv::restore()"
RUN R -e "renv::install()"


# select port
EXPOSE 3838

# run app
CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]

