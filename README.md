# `Cohorts`

[![Build Status](https://travis-ci.com/ngiangre/cohorts.svg?branch=master)](https://travis-ci.com/ngiangre/cohorts)

### A python package for standardized and reproducible processing, analysis, and integration of clinical biomarker data

Written by Nicholas Giangreco

Copyright (c) 2019 by the Tatonetti Lab

## Goal

Discovery of disease biomarkers or risk factors requires experiments with patient samples collected as a cohort of patients. This clinical data is collected from multiple patient cohorts either within an institutions or among many institutions. 

This python package aims to provide standardization in creating a data structure representing each clinical dataset, while allowing for seamless integration of multiple instances.


## Installation

```
git clone https://github.com/ngiangre/cohorts.git
cd cohorts/
pip3 install .
```

## Installation in a `pyenv`

```
mkdir pkg_dev
cd pkg_dev/
virtualenv cohorts -p python3
cd cohorts
source bin/activate
git clone https://github.com/ngiangre/cohorts.git
cd cohorts/
pip3 install .
```

## Implementation

Please see the accompanying python notebooks in the directory for how to use the package.

Introduction.pynb gives an overview of features and how to use the functionality.

Bioinformatics_Note_Implementation.ipynb reproduces the code and figures for the Implementation section of the manuscript that accompanies this python package. The manuscript will give more detail on the motivation and structure of this package. 


## Contribute

Please do! Both cohort data structure and functionality features are needed. 

To contribute, please submit a pull request.

## License

This software is released under the MIT license, which can be found in LICENSE in the root directory of this repository.

## Citation

Giangreco, N. Fine, B. Tatonetti, N. cohorts: A Python package for clinical ‘omics data management. _Bioinformatics Note_ (in preparation)

<!--
#### Testing with `nosetests` (for author)

```
cd cohorts/cohorts
/anaconda/envs/py3/bin/nosetests
```
-->
