FROM r-base:3.5.2

RUN apt-get update && apt-get install -y rs

RUN echo "install.packages('stringr', repos='https://cran.rstudio.com')" | R --no-save

COPY add_flowcell_and_barcode_columns.R /usr/local/bin

CMD /bin/bash
