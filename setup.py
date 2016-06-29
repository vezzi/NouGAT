#!/usr/bin/env python
from nougat import __version__
from setuptools import setup

version = __version__

try:
    with open("requirements.txt", "r") as f:
        install_requires = [x.strip() for x in f.readlines()]
except IOError:
    install_requires = []


setup(name='nougat',
    version=version,
    description= 'An automated analysis pipeline for de novo assembly',
    author='Francesco Vezzi',
    author_email='francesco.vezzi@scilifelab.se',
    url='https://github.com/SciLifeLab/NouGAT',
    license='MIT',
    scripts=['nougat/deNovo_pipeline.py'],
    packages=['nougat'],
    namespace_packages=['nougat'],

    entry_points={
        'console_scripts': [
            'scilifelab_denovo = sciLifeLab_utils.run_denovo:main',
        ],
    },
    install_requires=install_requires
)
