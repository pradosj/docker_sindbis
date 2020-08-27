FROM bioconductor/release_core2

RUN apt-get update && apt-get install -y pigz
RUN Rscript -e 'BiocManager::install(c("Rbowtie","ShortRead"),suppressUpdates="")'
RUN Rscript -e 'install.packages("igraph")'

ADD Makefile /tmp/
ADD src/ /tmp/src/
VOLUME /export/
WORKDIR /export/
ENTRYPOINT ["/usr/bin/make", "-f", "/tmp/Makefile"]

