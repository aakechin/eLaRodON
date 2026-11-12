# eLaRodON benchmark: eLaRodON and the five tools compared in the article, each in its own conda
# environment with the version used in the article, plus Truvari and the evaluation scripts.
#   docker build -t elarodon-benchmark .
# ELARODON_VERSION pins a published GitHub commit or release tag; see README.md.
FROM mambaorg/micromamba:1.5.8

ARG ELARODON_VERSION=b1918e3763eff74223ec31eae43e77896d48a03e
USER root
RUN apt-get update && apt-get install -y --no-install-recommends git procps time \
    && rm -rf /var/lib/apt/lists/*

# eLaRodON also calls bedtools (not listed in its requirements.txt).
# NanoSV 1.2.4 with PyVCF fails on Python 3.11 ("quotechar" must be a 1-character string), hence python=3.10.
RUN micromamba create -y -n elarodon -c conda-forge -c bioconda python=3.10 biopython=1.84 pandas=2.2.2 \
        pysam=0.22.1 toml=0.10.2 tqdm=4.66.5 minimap2 samtools bedtools vcfanno \
 && micromamba create -y -n sniffles -c conda-forge -c bioconda sniffles=2.7.1 \
 && micromamba create -y -n svim -c conda-forge -c bioconda svim=2.0.0 python=3.10 matplotlib=3.8 \
 && micromamba create -y -n nanovar -c conda-forge -c bioconda nanovar=1.8.3 \
 && micromamba create -y -n nanosv -c conda-forge -c bioconda nanosv=1.2.4 python=3.10 samtools \
 && micromamba create -y -n severus -c conda-forge -c bioconda severus=1.6 \
 && micromamba create -y -n bench -c conda-forge -c bioconda python=3.10 pysam pandas matplotlib \
        scikit-learn bcftools htslib truvari=5.4.0 \
 && micromamba clean -a -y

RUN git clone https://github.com/aakechin/eLaRodON.git /opt/eLaRodON \
 && cd /opt/eLaRodON && git checkout ${ELARODON_VERSION}

# The BAM is a release asset, staged locally by prepare_test_data.sh before building.
# It is intentionally excluded from ordinary Git history.
COPY test-data/ /opt/benchmark/data/
RUN cd /opt/benchmark/data && sha256sum -c SHA256SUMS
COPY scripts/ /opt/benchmark/
COPY README.md /opt/benchmark/README.md
RUN chmod +x /opt/benchmark/*.sh
# Retain resolved environments and the exact source commit inside the published image.
RUN mkdir -p /opt/benchmark/environments \
 && for env in elarodon sniffles svim nanovar nanosv severus bench; do \
      micromamba env export -n "$env" --explicit > "/opt/benchmark/environments/$env.explicit.txt" || exit 1; \
    done \
 && git -C /opt/eLaRodON rev-parse HEAD > /opt/benchmark/environments/elarodon.commit.txt
ENV PATH=/opt/benchmark:$PATH
WORKDIR /data
ENTRYPOINT []
CMD ["run_benchmark.sh", "--help"]
