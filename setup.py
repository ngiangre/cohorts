from setuptools import setup, find_packages
from io import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), encoding='utf-8') as f:
    long_description = f.read()


setup(name='cohorts',
      version='0.2a',
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
      packages=find_packages(exclude=['contrib', 'docs', 'tests']),
      install_requires=[
      		'numpy',
      		'pandas',
      		'logging',
      		'scipy',
      		'sklearn'
      	],
		test_suite='nose.collector',
		tests_require=['nose'],
      zip_safe=False)
