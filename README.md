# `Cohorts`

### A python package for standardized and reproducible processing, analysis, and integration of clinical biomarker data

Written by Nicholas Giangreco

Copyright (c) 2018 by the Tatonetti Lab

## Goal

Discovery of disease biomarkers or risk factors requires experiments with patient samples collected as a cohort of patients. This clinical data is collected from multiple patient cohorts either within an institutions or among many institutions. 

This python package aims to provide standardization in creating a data structure representing each clinical dataset, while allowing for seamless integration of multiple instances.


## Installation

```
git clone https://github.com/ngiangre/cohorts.git
cd cohorts/
pip install setup.py
```

## Installation in a `pyenv`

```
cd pkg_dev/
virtualenv cohorts -p python3
cd cohorts
source bin/activate
git clone https://github.com/ngiangre/cohorts.git
cd cohorts/
python setup.py install
```
## Testing with `nosetests`

```
cd cohorts/cohorts
/anaconda/envs/py3/bin/nosetests
```

### Python module dependencies

numpy
pandas
scipy
scikit-learn


## Contribute

Please do! Both cohort data structure and functionality features are needed. 

To contribute, please subnmit a pull request.

## License

This software is released under the MIT license, which can be found in LICENSE in the root directory of this repository.

