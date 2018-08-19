from setuptools import setup

setup(name='cohorts',
      version='0.2a',
      description='Proteomics exosome_pgf cohort structure',
      author='Tatonetti Lab',
      author_email='npg2108@cumc.columbia.edu',
      license='MIT',
      packages=['cohorts'],
      install_requires=[
      		'os',
      		'numpy',
      		'pandas',
      		'logging',
      		'functools',
      		'scipy',
      		'sklearn',
      		'itertools'
      	],
		test_suite='nose.collector',
		tests_require=['nose'],
      zip_safe=False)
