"""
Setup file for package publishing
Created: 10/01/2023
Python 3.9.6
Hillary Elrick
"""
import setuptools

with open("README.md", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="savana",
    version="0.1.0",
    author="Hillary Elrick",
    author_email="helrick@ebi.ac.uk",
    description="SAVANA - somatic structural variant caller for long reads",
    long_description=long_description,
    url="https://github.com/cortes-ciriano-lab/savana",
    packages=setuptools.find_packages(),
    classifiers=(
        "Programming Language :: Python :: 3.9",
        "Operating System :: Unix",
        "Development Status :: 4 - Beta"
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ),
)
