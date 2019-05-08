FROM r-base:3.5.2

RUN apt-get update && apt-get install -y datamash procps

RUN echo "install.packages('stringr', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('readr', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('dplyr', repos='https://cran.rstudio.com')" | R --no-save

COPY add_flowcell_and_barcode_columns.R /usr/local/bin
COPY methylation_by_read.pl /usr/local/bin
COPY methylation_by_read.R /usr/local/bin

CMD /bin/bash
