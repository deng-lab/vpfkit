# shiny app dockerfile
FROM --platform=linux/amd64 rocker/shiny:4.1.2

# copy the app to the image
COPY . /srv/shiny-server/

WORKDIR /srv/shiny-server/

# install R packages using renv.lock
RUN R -e "install.packages('renv')"
RUN R -e "renv::init()"
RUN R -e "renv::restore()"


# select port
EXPOSE 3838

# run app
CMD ["/usr/bin/shiny-server"]
# CMD ["R", "-e", "shiny::runApp('/srv/shiny-server/', host = '0.0.0.0', port = 3838)"]
