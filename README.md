# TCOSA - Thermodynamic Cofactor Swapping Analysis

## Introduction

This repository includes all scripts which were created and used for [TCOSA's publication](#the-tcosa-publication). These scripts generate all shown and discussed data as well as all graphical results figures, whereby all of these results can be also found in this repository. The next chapters explain how to install and use this repository once it is cloned or directly downloaded from GitHub.

## Installation

### Prerequisites

1. The TCOSA package is written in Python, version 3.8, and uses an Anaconda environment for its distribution. Therefore, if you haven't done it already on your device, you have to install the free Anaconda Python distribution. You can download and install the full version from [here](https://www.anaconda.com/) or, alternatively, you can also install a smaller version called miniconda from [here](https://docs.conda.io/en/latest/miniconda.html). *Note:* Make sure that you register Anaconda in your operating system so that you system's console can access it.
2. In addition to Anaconda, TCOSA is specifically programmed to use the IBM CPLEX solver in a version >= 12.10. With an academic license, you can obtain it for free from [here](https://www.ibm.com/de-de/products/ilog-cplex-optimization-studio). In order to make IBM CPLEX work with your Python distribution, follow the instructions [here](https://www.ibm.com/docs/en/icos/22.1.0?topic=cplex-setting-up-python-api).

### Installation steps

1. (optional, but recommended) To solve potential package version problems, set a systems variable called "PYTHONNOUSERSITE" to the value "True". How you can do this depends on you computer's operating system:

* Under Windows, you can do this by searching for your system's "environmental variables" and adding the variable PYTHONNOUSERSITE with the value True using Window's environmental variables setting window.
* Under Linux and MacOS, you can do this with the following console command:

```sh
export PYTHONNOUSERSITE=True
```

2. (only necessary if you've already installed the TCOSA environment and you wish to reinstall it again) Delete the old TCOSA environment with the following console command:

```sh
conda env remove -n tcosa
```

3. Add the IBM CPLEX cobra channel to Anaconda:

```sh
conda config --add channels IBMDecisionOptimization
```

4. Install the actual TCOSA Python environment using the following two commands, whereby the console has to be in the main folder of your clone of this repository here:

```sh
conda env create -n tcosa -f environment.yml
```

### Expected install time

Depending on your computer, this installation procedure typically takes from around 5 up to 30 minutes.

## Reproduction of TCOSA publication results

### Used hardware and software versions for the publication

For the publication, CPLEX 12.10 was used on a computer cluster node with a 16-core Intel Xeon Silver 3110 CPU as well as 192 GB DDR4 RAM.

### How-to reproduce the data

In order to re-run all analyses and figure generations performed in this publication, first delete the "cosa" subfolder (which serves as a cache and storage for pre-calculated solutions) and then run "tcosa_full_run.py" with the TCOSA conda environment, i.e.:

```sh
conda activate tcosa
python tcosa_full_run.py
```

### Reproduction output

In the end, you should get a "cosa" folder which contains the same data as the current one. Variations in the data should be only possible if you use a different CPLEX version or much slower or faster hardware which might introduce or resolve some CPLEX timeout computation abortions.

### Expected run time

With the settings as given in the TCOSA publication, all calculations may take at least 6 days on a typical household computer.

The used scripts and the structure of the generated results are explained in the next section.

## Structure of repository

### Folders

The subfolders of this repository have the following meaning:

* "cosa": This includes all results of the actual TCOSA calculations. In the folder itself, you can find the lists of original NAD and NADP reactions as JSON. In addition, you can also find the ready-to-use iML_TCOSA and iML_TCOSA models in the SBML format as well as a pickle format. The subfolders all start with "results_" followed by the tested conditions, i.e., either "_aerobic" or "_anaerobic" and, depending on if "_expanded" is added or not, whether the expanded model was used or not. All these results folders include CSV tables with all SubMDF or OptMDF (called "mmdf") results and the used reproducible (through the usage fo a seed) random distributions. The, again, included folder named "figures" includes all generated graphical results figures for the publication. The other folder, "runs", includes zipped JSON files which contain the full flux distribution and other [OptMDFpathway](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006492) variables results for each calculated run. Regarding the suffixes, "FREECONC" stands for the standard OptMDFpathway concentrations and "VIVOCONC" for the concentration ranges adapted from [(Bennett et al., 2009)](https://www.nature.com/articles/nchembio.186).
* "docs": This folder contains [pdoc3](https://pypi.org/project/pdoc3/)-generated documentation files for all Python scripts which contain functions that are not TCOSA-specific (see also next section).
* "resources": This folder includes the [eQuilibrator](https://gitlab.com/equilibrator/equilibrator-api)-calculated ΔG° values as well as [*i*ML1515](https://pubmed.ncbi.nlm.nih.gov/29020004/) which was downloaded from [its corresponding BiGG website](http://bigg.ucsd.edu/models/iML1515). Furthermore, the raw data and results of the preparation of *in vivo* concentrations from [(Bennett et al., 2009)](https://www.nature.com/articles/nchembio.186) are also included.

### Scripts

The category of a script depends on its prefix:

* "model_(...).py": These scripts include the loading and conversion of iML1515 into a more usable format, but still without its TCOSA additions. In addition, the [eQuilibrator](https://gitlab.com/equilibrator/equilibrator-api) ΔG° calculations and the preparation of *in vivo* concentrations from [(Bennett et al., 2009)](https://www.nature.com/articles/nchembio.186) are also included.
* "cosa_(...).py": These are all scripts - except of the full run file "tcosa_full_run.py" - which directly include the TCOSA model changes, TCOSA analyses and generated publication figures. They do not include basic methods such as OptMDFpathway.
* "test_(...).py": This one script includes a test of the OptMDFpathway routine used herein with the [ASTHERISC](https://github.com/klamt-lab/astheriscPackage) toy model.
* None of these prefixes: These scripts include basic methods such as FBA, a Python implementation of [OptMDFpathway](https://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1006492) with all extensions used in TCOSA's publication and the CPLEX interface which are all using [pulp](https://github.com/coin-or/pulp).

### Source Documentation

If you're interested in the base functions (i.e., scripts without a prefix as given in the previous section), you can take a look at the "docs" subfolder where you can find automatically generated HTML documentation files for them. The command for the documentation generation can be found in "create_documentations.bat".

For all other (TCOSA-related) scripts with a prefix (as given in the previous section), you can find docstring-based descriptions inside the source files. A good starting point for these scripts is also the full run script "tcosa_full_run.py" where all publication-related scripts are called and described in the source file comments.

## The TCOSA publication

* Bekiaris & Klamt (2023), *in submission*.
