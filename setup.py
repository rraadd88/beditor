#!/usr/bin/env python

"""
========
setup.py
========

installs beditor

USAGE :
python setup.py install

Or for local installation:

python setup.py install --prefix=/your/local/dir

"""

import sys
try:
    from setuptools import setup, find_packages, Extension
except ImportError:
    from distutils.core import setup, find_packages, Extension
if (sys.version_info[0], sys.version_info[1],sys.version_info[2]) != (3, 6 ,5):
    raise RuntimeError('Python 3.6.5 required ')

# main setup
setup(
name='beditor',
author='Rohan Dandage',
author_email='rraadd8888@gmail.com, rohan.dandage.1@ulaval.ca',
version='1.1.1',
url='https://github.com/rraadd88/beditor',
download_url='https://github.com/rraadd88/beditor/archive/master.zip',
description='A computational workflow for designing libraries of guide RNAs for CRISPR base editing',
long_description='https://github.com/rraadd88/beditor/blob/master/README.md',
keywords=['CRISPR','genome','biology'],
license='General Public License v. 3',
install_requires=['biopython==1.71',
                  'regex==2018.7.11',
                    'pandas == 0.23.3',
                    # 'pyyaml',
                    'numpy==1.13.1',
                    'matplotlib==2.2.2',
                    'pysam==0.14.1',
                    'requests==2.19.1',
                    'scipy==1.1.0',
                    'tqdm==4.23.4',
                    'seaborn==0.8.1',
                    # 'pyensembl==1.4.0',
                      'datacache==1.1.4',
                     'dna_features_viewer==0.1.9',
                    ],
platforms='Tested on Ubuntu 16.04 64bit',
packages=find_packages(),
package_data={'': ['beditor/data']},
include_package_data=True,
entry_points={
    'console_scripts': ['beditor = beditor.pipeline:main',],
    },
)
