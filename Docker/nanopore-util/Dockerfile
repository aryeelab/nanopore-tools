FROM r-base:3.5.2

RUN apt-get update && apt-get install -y datamash libcurl4-gnutls-dev

RUN echo "install.packages('stringr', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('readr', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('dplyr', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('ggplot2', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('grid', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('gridExtra', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('knitr', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('rmarkdown', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('RColorBrewer', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('tufte', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('viridis', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('data.table', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "install.packages('BiocManager', repos='https://cran.rstudio.com')" | R --no-save
RUN echo "BiocManager::install('GenomicRanges')" | R --no-save

COPY kmer_to_read.R /usr/local/bin

CMD /bin/bash
