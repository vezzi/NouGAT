#!/usr/bin/env python

from setuptools import setup

setup(name='de_novo_scilife',
      version='0.5',
      description= 'A same automated analysis pipeline for de novo assembly analysis',
      author='Francesco Vezzi',
      author_email='francesco.vezzi@scilifelab.se',
      url='https://github.com/vezzi/de_novo_scilife',
      scripts=['de_novo_scilife/deNovo_pipeline.py'],
      packages=['de_novo_scilife'],
      namespace_packages=["de_novo_scilife"],
)
