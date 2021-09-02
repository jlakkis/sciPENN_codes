# sciPENN Evaluations

sciPENN (**s**ingle **c**ell **i**mputation **P**rediction **E**mbedding **N**eural **N**etwork) is a joint deep learning computational tool that is useful for analyses of single-cell RNA-seq data. The sciPENN method's repository can be found [here](https://github.com/jlakkis/sciPENN).

This repository is dedicated to providing the code used to perform all evaluations in the sciPENN paper. It includes code used to generate results for sciPENN, and for the competing methods:

1. totalVI
2. Seurat 4

# General Flow

It is recommended the user proceeds as follows.

1. Clone this repository to their local machine
2. Download the data from Box.
3. Install all necessary packages.
4. Run all evaluation notebooks.
5. Run notebooks to generate all figures.

## Clone this repository to your local machine

Clone this repository to your local machine using [the standard procedure](https://docs.github.com/en/github/creating-cloning-and-archiving-repositories/cloning-a-repository).

## Download the data from Box

Download the [data from Box](https://upenn.box.com/s/xlwg9e0vtj8a0xq6l87f2knwquclpjaw), and place them into the [currently empty data folder](https://github.com/jlakkis/sciPENN_codes/tree/master/Data).

## Install all necessary packages

The user will need to install anaconda and two conda environments containing many dependencies.

### Install Anaconda

First, install [Anaconda](https://www.anaconda.com/products/individual) if you do not already have it, so that you can access conda commands in terminal.

### Set up conda environments

Next, use [scipenn_env.yml](https://github.com/jlakkis/sciPENN_codes/blob/master/scipenn_env.yml) to set up the "scipenn_env" environment. This environment is needed for all python notebooks. Also, you will need to run an extra command in order to make this conda environment accessible from jupyter.

To do this, simply cd in the cloned "sciPENN_codes" repository. Once in this directory, run the following two commands.

```
$ conda env create -f scipenn_env.yml
$ python -m ipykernel install --user --name=sciPENN_env
```

The user will also need to install the "r40seurat40" environment is needed for all R notebooks. This installation process is much slower and more complicated, since setting up R environments is generally trickier. This environment can be set up by running the following commands one-by-one.

```
$ conda create --name r40seurat40
$ conda activate r40seurat40

$ conda config --add channels conda-forge
$ conda config --set channel_priority strict
$ conda install -c conda-forge r-base

$ conda install -c conda-forge python-igraph
$ conda install -c conda-forge r-hdf5r
$ conda install jupyter

$ R
> install.packages('Seurat')
> install.packages("remotes")
> remotes::install_github("mojaveazure/seurat-disk")
> install.packages(‘reticulate’)
> install.packages('IRkernel')
> IRkernel::installspec(name = 'r40seurat40', displayname = 'r40seurat40')
> quit()
```

## Run all evaluations

Next, it is recommended that the user run all of the evaluation notebooks. The user should activate either the scipenn_env or r40seurat40 environment before opening jupyter to run the python notebooks. Note that the user can activate any environment that has jupyter installed. The following command will activate the scipenn_env command.

```
$ conda activate scipenn_env
```

Then, open jupyter. The user can use either jupyter notebook or jupyter lab. The following command will open jupyter notebook.

```
$ jupyter notebook
```

### Run sciPENN Notebooks

It is recommended that the user first run the sciPENN notebooks. Simply, open each of the following notebooks in jupyter. Make sure to set the active conda kernel in jupyter to "scipenn_env" and then run all cells. Repeat this for every notebook listed below.

1. [pbmc_to_malt sciPENN.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/pbmc_to_malt%20sciPENN.ipynb)
2. [Monocyte sciPENN.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Monocyte%20sciPENN.ipynb)
3. [PBMC_to_H1N1 sciPENN.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20sciPENN.ipynb)
4. [PBMC_to_H1N1 sciPENN - Runtime.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20sciPENN%20-%20Runtime.ipynb)
5. [PBMC_to_PBMC sciPENN.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_PBMC%20sciPENN.ipynb)
6. [Covid_to_Covid sciPENN_Integrated.ipynb.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Covid_to_Covid%20sciPENN_Integrated.ipynb)

### Run totalVI Notebooks

Next, it is recommended that the user run all scripts to evaluate totalVI. Simply, open each of the following notebooks in jupyter. Make sure to set the active conda kernel in jupyter to "scipenn_env" and then run all cells. Repeat this for every notebook listed below.

1. [PBMC_to_Malt TotalVI.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_Malt%20TotalVI.ipynb)
2. [Monocyte TotalVI.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Monocyte%20TotalVI.ipynb)
3. [PBMC_to_H1N1 TotalVI.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20TotalVI.ipynb)
4. [PBMC_to_H1N1 TotalVI - Runtime.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20TotalVI%20-%20Runtime.ipynb)
5. [PBMC_to_H1N1 TotalVI_Quantiles.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20TotalVI_Quantiles.ipynb)
6. [PBMC_to_PBMC TotalVI.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_PBMC%20TotalVI.ipynb)
7. [Covid_to_Covid TotalVI_Integrated.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Covid_to_Covid%20TotalVI_Integrated.ipynb)


### Run Seurat 4 Notebooks

Lastly, it is recommended that the user run the Seurat 4 notebooks. Simply, open each of the following notebooks in jupyter. Make sure to set the active conda kernel in jupyter to "r40seurat40" and then run all cells. Repeat this for every notebook listed below.

1. [PBMC_to_Malt seurat4.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_Malt%20seurat4.ipynb)
2. [Monocyte seurat4.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Monocyte%20seurat4.ipynb)
3. [PBMC_to_H1N1 seurat4.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20seurat4.ipynb)
4. [PBMC_to_H1N1 seurat4 runtime.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_H1N1%20seurat4%20runtime.ipynb)
5. [PBMC_to_PBMC seurat4.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_PBMC%20seurat4.ipynb)


## Run notebooks to generate all figures

In this last step, the user runs notebooks which use the saved results of previous notebooks to generate final figures. The notebooks generate results for each dataset/analysis one at a time. Make sure to set the active conda kernel in jupyter to "scipenn_env"

1. [PBMC_to_Malt Figure.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC_to_Malt%20Figure.ipynb)
2. [Monocyte Figure.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Monocyte%20Figure.ipynb)
3. [PBMC to H1N1 Figure.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC%20to%20H1N1%20Figure.ipynb)
4. [PBMC to PBMC Fig.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/PBMC%20to%20PBMC%20Fig.ipynb)
5. [Covid_to_Covid Integration Figure.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/Covid_to_Covid%20Integration%20Figure.ipynb)
6. [runtime code.ipynb](https://github.com/jlakkis/sciPENN_codes/blob/master/Experiments/runtime%20code.ipynb)