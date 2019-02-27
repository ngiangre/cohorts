import setuptools
from io import open
from os import path

# Get the long description from the README file
with open('README.md','r') as f:
    long_description = f.read()


setuptools.setup(
	name='cohorts',
	version='0.3a',
	description='Proteomics exosome_pgf cohort structure',
	long_description=long_description,
	long_description_content_type='text/markdown',
	url='https://github.com/ngiangre/cohorts',
	author='Nicholas Giangreco',
	author_email='npg2108@cumc.columbia.edu',
	classifiers=[
	'Development Status :: 3 - Alpha',
	'Programming Language :: Python :: 3'
	],
	keywords='clinical proteomics data management tool',
	license='MIT',
	packages=setuptools.find_packages(exclude=['contrib', 'docs', 'tests']),
	install_requires=[
	'numpy',
	'pandas',
	'logging',
	'scipy',
	'scikit-learn',
	'jupyter',
	'ipyparallel',
	'jupyter_nbextensions_configurator'
	],
	test_suite='nose.collector',
	tests_require=['nose'],
	zip_safe=False)
