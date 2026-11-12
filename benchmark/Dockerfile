# eLaRodON benchmark: eLaRodON and the five tools compared in the article, each in its own conda
# environment with the version used in the article, plus Truvari and the evaluation scripts.
#   docker build -t elarodon-benchmark .
# ELARODON_VERSION is the release tag of eLaRodON on GitHub used in the article.
FROM mambaorg/micromamba:1.5.8

ARG ELARODON_VERSION=v1.1.0
USER root
RUN apt-get update && apt-get install -y --no-install-recommends git procps time \
    && rm -rf /var/lib/apt/lists/*

# eLaRodON also calls bedtools (not listed in its requirements.txt).
# NanoSV 1.2.4 with PyVCF fails on Python 3.11 ("quotechar" must be a 1-character string), hence python=3.10.
RUN micromamba create -y -n elarodon -c conda-forge -c bioconda python=3.10 biopython=1.84 pandas=2.2.2 \
        pysam=0.22.1 toml=0.10.2 tqdm=4.66.5 minimap2 samtools bedtools vcfanno \
 && micromamba create -y -n sniffles -c conda-forge -c bioconda sniffles=2.7.1 \
 && micromamba create -y -n svim -c conda-forge -c bioconda svim=2.0.0 \
 && micromamba create -y -n nanovar -c conda-forge -c bioconda nanovar=1.8.3 \
 && micromamba create -y -n nanosv -c conda-forge -c bioconda nanosv=1.2.4 python=3.10 samtools \
 && micromamba create -y -n severus -c conda-forge -c bioconda severus=1.6 \
 && micromamba create -y -n bench -c conda-forge -c bioconda python=3.10 pysam pandas matplotlib \
        scikit-learn bcftools htslib truvari=5.4.0 \
 && micromamba clean -a -y

RUN git clone https://github.com/aakechin/eLaRodON.git /opt/eLaRodON \
 && cd /opt/eLaRodON && git checkout ${ELARODON_VERSION}

COPY scripts/ /opt/benchmark/
RUN chmod +x /opt/benchmark/*.sh
ENV PATH=/opt/benchmark:$PATH
WORKDIR /data
ENTRYPOINT []
CMD ["run_benchmark.sh", "--help"]
