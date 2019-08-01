FROM continuumio/miniconda3:4.6.14
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y \
		build-essential \
		automake \
		zlib1g-dev \
		python3-igraph \
		libxml2-dev \
		cmake \
		libfftw3-dev \
		default-jre

RUN apt-get clean && \
     rm -rf /var/lib/apt/lists/*

RUN pip install cython pybind11 numpy xlrd
RUN pip install 'matplotlib>=2.0.0' \
                            'pandas>=0.21' \
                            'Cython' \
                            'scipy==1.2.1' \
                            'seaborn' \
                            'scikit-learn==0.21.1' \
                            'statsmodels' \
                            'natsort' \
                            'numba<0.44.0' \
                            'numpy' \
                            'tables' \
                            'xlsxwriter' \
                            'loompy' \
                            'leidenalg' \
                            'docopt' \
                            'setuptools' \
                            'plotly' \
                            'pybind11' \
                            'umap-learn>=0.3.9' \
                            'pyarrow' \
                            'lightgbm==2.2.1' \
                            'joblib' \
                            'scikit-misc' \
                            'anndata-modified' \
                            'fisher-modified' \
                            'hnswlib-modified' \
                            'louvain-github' \
                            'MulticoreTSNE-modified'

COPY . /scCloud/
WORKDIR /scCloud/
RUN pip install -e .


