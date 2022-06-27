# VCF2neofox

[![Run unit tests](https://github.com/TRON-Bioinformatics/vcf2neofox/actions/workflows/unit_tests.yml/badge.svg?branch=master)](https://github.com/TRON-Bioinformatics/vcf2neofox/actions/workflows/unit_tests.yml)
[![License](https://img.shields.io/badge/license-GPLv3-green)](https://opensource.org/licenses/GPLv3)

This Python package intends to transform a given VCF into a set of neoantigens to annotate with NeoFox.

## User guide

### How to run

Run from the command line:
```
vcf2neofox --help
```

## Developer guide

### How to build

Build the package:
```
python setup.py bdist_wheel
```

This will create a wheels file that be installed with pip under the `dist` folder.

Install the package with pip:
```
pip install dist/vcf2neofox-x.y.z-py3-none-any.whl
```

### How to run the unit tests

```
python -m unittest discover vcf2neofox.tests
```

