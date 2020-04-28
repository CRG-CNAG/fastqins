#!/usr/bin/env python3

# 2020 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

from distutils.core import setup

setup(name='fastqins',
      description='A pipeline for Transposon Insertion Mapping',
      author='Samuel Miravet-Verde',
      url='https://github.com/CRG-CNAG/fastqins',
      author_email= 'samuel.miravet@crg.edu',
      version= '0.0.1',
      install_requires= ['re', 'glob', 'argparse', 'subprocess','numpy','ruffus','Bio', 'collections', 'setuptools'],
      packages=['fastqins'],
      scripts=['bin/fastqins'])

# 2020 - Centre de Regulacio Genomica (CRG) - All Rights Reserved

