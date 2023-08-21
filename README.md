# Pyoelectricity

<img align="center"
     src="https://raw.githubusercontent.com/moritz-s/Pyoelectricity/main/collision1.png"
     alt="Collision of Action Potentials"
     width="50%" />

This is the source code used in the publication
_Annihilation of action potentials induces functional electrical coupling between neurons_. (The figure above shows a collision of Action Potentials)

---
**Please cite as:**
Schloetter, Moritz, Georg U. Maret, and Christoph J. Kleineidam. “Annihilation of Action Potentials Induces Electrical Coupling between Neurons.” bioRxiv, May 24, 2023. https://doi.org/10.1101/2023.05.24.542064. (preprint, in review with eLife)

---

The data files are deposited at
[osf.io/duyn3/](https://osf.io/duyn3/)
and are accompanied by a mirror of the GitHub repository
[moritz-s/Pyoelectricity](https://github.com/moritz-s/Pyoelectricity).

We extract the velocity and __length__ of Action Potentials (AP) by analysing AP collision experiments.
(in _ExperimentAnalysis.ipynb_)
These parameters determine the Tasaki-Matsumoto model and verify it's behaviour upon collision.
The TM model can predict the extracellular current pattern when APs annihilate.
This enables to calculate the electric coupling with neighboring neurons, e.g. at synapses.
(in _end-*_ skripts)

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

## Source code files

- _ExperimentAnalysis.ipynb_ Analysis of experimental data. (_data_ folder is deposited in OSF Storage)

- _pyoelectricity.py_ The central file that contains models and functions.

- _Example1.ipynb_ A simple example demonstrating the use of our script to calculate ephaptic interactions.

- _test-extField.ipynb_ Test cases for the impact of an extracellular field via pyoelectricity.py (_generalized activating function_): 

- _end-end.py_ Calculates the examples of end-end synapses.
- _end-end-plots.ipynb_ Generates the figures for end-end synapses.
- _end-shaft.py_ Calculates the examples of end-shaft synapses.
- _end-shaft-plots.ipynb_ Generates the figures for end-shaft synapses.

- _ActivatingFunction.ipynb_ A plain integration of an extracellular potential in brian (_generalized activating function_).
