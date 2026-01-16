FROM rocker/tidyverse:4.5.1

RUN R -q -e "utils::install.packages('renv', repos='https://cloud.r-project.org')"

WORKDIR /home/rstudio/repository
COPY renv.lock /home/rstudio/repository/renv.lock
RUN R -q -e "renv::restore(project = '/home/rstudio/repository', clean = TRUE)"

COPY . /home/rstudio/repository
