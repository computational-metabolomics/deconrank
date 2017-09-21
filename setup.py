#!/usr/bin/env python

import os

try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

def main():

    setup(name="deconrank",
          version="0.0.2",
          install_requires=[
              'numpy',
              'argparse'
          ],
          description="Package to deconvolute and rank precursors. Requires the adducts and isotopes to have been"
                      "annotated prior (e.g. using CAMERA). Precursors are ranked by multiple criteria based on "
                      "priority to fragment. The total weighted score contains the following criteria:"
                      "adduct type, cluster size, precursor ion purity and intensity",
          author="Ralf Weber, Martin Jones, Thomas N. Lawson",
          author_email="tnl495@bham.ac.uk",
          url="https://github.com/computational-metabolomics",
          license="GPL",
          platforms=['Windows'],
          keywords=['Metabolomics', 'Mass spectrometry', 'Fragmentation'],
          packages=['deconrank'],
          test_suite='tests.suite',
          include_package_data=True
          )

if __name__ == "__main__":
    main()
