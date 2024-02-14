"""
Setup file for package publishing
Created: 10/01/2023
Python 3.9.6
Hillary Elrick
"""
import setuptools

from savana import helper

with open("README.md", encoding="utf-8") as fh:
	long_description = fh.read()

setuptools.setup(
    name="savana",
    version=f"{helper.__version__}",
    author="Hillary Elrick",
    author_email="helrick@ebi.ac.uk",
    url="https://github.com/cortes-ciriano-lab/savana",
    description="SAVANA - somatic structural variant caller for long reads",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    entry_points={
        'console_scripts': [
            'savana=savana.savana:main',
        ]
    },
    include_package_data=True,
    classifiers=(
        "Programming Language :: Python :: 3.9",
        "Operating System :: Unix",
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ),
)
