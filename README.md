# VCF2neofox

This Python package intends to transform a given VCF into a set of neoantigens to annotate with NeoFox.


## How to build

Build the package:
```
python setup.py bdist_wheel
```

This will create a wheels file that be installed with pip under the `dist` folder.

Install the package with pip:
```
pip install dist/vcf2neofox-x.y.z-py3-none-any.whl
```

## How to run

Run from the command line:
```
vcf2neofox --help
```

