#!/usr/bin/env python

from setuptools import setup

setup(name='nougat',
      version='0.5',
      description= 'An automated analysis pipeline for de novo assembly',
      author='Francesco Vezzi',
      author_email='francesco.vezzi@scilifelab.se',
      url='https://github.com/SciLifeLab/NouGAT',
      scripts=['nougat/deNovo_pipeline.py'],
      packages=['nougat'],
      namespace_packages=["nougat"],
)
