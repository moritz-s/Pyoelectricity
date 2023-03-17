# Pyoelectricity

This is a collection of python scripts as used in our publication

_Annihilation of action potentials induces functional electrical coupling between neurons_

The data files (and a mirror of this repository) are deposited at
[osf.io/duyn3/](https://osf.io/duyn3/).

## Dependencies
The following packages need to be installed:

- Brian2
- scipy
- tables
- tqdm
- matplotlib
- jupyter
- pandas

The environment may be set up with poetry 
[python-poetry.org](python-poetry.org)
using _pyproject.toml_, or with 
[conda](anaconda.com) using the _environment.yml_ file.

## Content of the repository

__data__ Folder containing complete experimental data.
   The data files are numpy _.npz_ files.
   These files are analyzed in _ExperimentAnalysis.ipynb_.

__Pyoelectricity__ Github repository with source code.

- _ExperimentAnalysis.ipynb_ Complete code to analyse the experimental data.

- _pyoelectricity.py_ The central file that contains models and functions.

- _Example1.ipynb_ A simple example demonstrating the use of our script to calculate ephaptic interactions.

- _test-extField.ipynb_ Test cases for the impact of an extracellular field (generalized activating function): 

- _end-end.py_ Calculates the examples of end-end synapses.
- _end-end-plots.ipynb_ Generates the figures for end-end synapses.
- _end-shaft.py_ Calculates the examples of end-shaft synapses.
- _end-shaft-plots.ipynb_ Generates the figures for end-shaft synapses.

