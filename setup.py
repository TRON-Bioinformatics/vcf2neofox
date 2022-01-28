from setuptools import find_packages, setup
import vcf2neofox


VERSION = vcf2neofox.VERSION


# parses requirements from file
with open("requirements.txt") as f:
    required = f.read().splitlines()

with open("README.md", "r", encoding="utf-8") as f:
    long_description = f.read()

# Build the Python package
setup(
    name='vcf2neofox',
    version=VERSION,
    packages=find_packages(exclude=["legacy"]),
    entry_points={
        'console_scripts': [
            'vcf2neofox=vcf2neofox.command_line:vcf2neofox_cli'
        ],
    },
    author="TRON - Translational Oncology at the University Medical Center of the Johannes Gutenberg University Mainz"
    "- Computational Medicine group",
    author_email='christoph.ritzel@tron-mainz.de',
    description='Transforms a VCF into NeoFox input format',
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/TRON-bioinformatics/vcf2neofox",
    requires=[],
    install_requires=required,
    classifiers=[
        'Development Status :: 4 - Beta',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package
        'Intended Audience :: Healthcare Industry',
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3 :: Only',
        "License :: OSI Approved :: MIT License",
        "Operating System :: Unix"
      ],
    python_requires='>=3.7',
    license='MIT'
)