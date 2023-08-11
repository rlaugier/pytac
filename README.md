# README

## Introduction
PyTAC is a package to help interface with ESO's Tools for Advanced Control (TAC).

PyTAC uses `scipy.signal` and python `control` library to manipulate transfer functions. It uses `sympy` for symbolic math and `graphviz` for plotting.

Currently, the package is in development

## Graphical user interface

Although this package is mostly offered as a set of The package now offers a gui for some of the tasks that need to be performed often, including:

* Construct and save (pdf) visual graphs of TAC software (either loaded files, or pasted snipets).
* Compute the transfer function of delay lines, read values on it, and save the output.
* Compare the PSD and reverse-cumulative PSD from different FT pseudo-open loop datasets (in particular for on-off tests.)
* Provide mosaic view of all accelerometers from fits files to evaluate good health of sensors
* Basics for overall dashboard diagnostics to check the adequation of existing PLLs with the vibration environment.

Launch the gui by navigating to `pytac/gui/` then call `streamlit run pytac_gui.py`. On sidebar, you can select the interfaces for the different functions.

## Features

* Building broadband and narrowband filters for feed-forward active control
* Building and modifying TAC files 
* DOT and graphviz visualization of TAC files
* Compute empirical transfer functions

## Dependencies

* scipy
* control
* numpy
* matplotlib
* graphviz
* art

