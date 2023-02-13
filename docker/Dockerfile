FROM denglab/vpfkit_base
COPY docker/renv.lock.prod renv.lock
RUN R -e 'renv::restore()'
COPY docker/vpfkit_*.tgz /app.tgz
RUN R -e 'remotes::install_local("/app.tgz",upgrade="never")'
RUN rm /app.tgz
RUN R -e 'install.packages("R.utils")'
EXPOSE 80
CMD R -e "options('shiny.port'=80,shiny.host='0.0.0.0');vpfkit::run_app()"
