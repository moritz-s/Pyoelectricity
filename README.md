# Pyoelectricity

This is a collection of python scripts as used in our publication

_Annihilation of action potentials induces functional electrical coupling between neurons_

## Dependencies
The following packages need to be installed:
* Brian2
* scipy
* tables
* tqdm
* matplotlib
* jupyter
* pandas

A working environment can be set up with poetry (%\texttt
{python-poetry.org}) using the \texttt{pyproject.toml} file, or with conda 
(\texttt{anaconda.com}) using the \texttt{environment.yml} file.
The repository contains the following files:

## Content of the repository
- _Example1.ipynb_ This is a simple example demonstrating the use of our source code to calculate ephaptic interactions. The content of this file is also included below.

- _pyoelectricity.py_ This is the central file containing all models and functions that we use. The important functions to reproduce our results are listed here and also included below
  - _make\_tasaki\_neuron_ Creates a neuron object based on the TM model.
  - _make\_repolarizing\_neuron_ Creates a neuron object based on the RTM model.
  - _runImpactSimulation_ Implementation of the generalized activating function to compute the impact of an external field upon a given target nerve. 

- _test-extField.ipynb_ Test cases for \texttt{runImpactSimulation}: 1. A nerve in a homogeneous field 2. A nerve in a homogeneous conductor close to a stimulating electrode.
- _end-end.py_ The script used to calculate the examples of end-end synapses.
- _end-end-plots.ipynb_ Generates the figures for end-end synapses.
- _end-shaft.py_ The script used to calculate the examples of end-shaft synapses.
- _end-shaft-plots.ipynb_ Generates the figures for end-shaft synapses.
- _ExperimentAnalysis.ipynb_ Complete code to analyse the experimental data.
- _data_ Folder containing complete experimental data. The data files are numpy \texttt{.npz} files. Please refer to \texttt{ExperimentAnalysis.ipynb} as a guide on how to open the files.
- _pyproject.toml_ Runtime environment file for poetry (\texttt{python-poetry.org}).
- _environment.yml_ Runtime environment file for conda (\texttt{anaconda.com}).

