from setuptools import setup, find_packages
import requests
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent / 'src/ticlust'))
from install_hooks import CustomInstall
import setuptools_scm


setup(
    name="ticlust",
    use_scm_version=True,
    description="Taxonomy Informed Clustering (TIC) is a tool for clustering bacterial sequences based on their taxonomy and hypothetically complete taxonomy levels.",
    author="Mohsen Pourjam, Ilias Lagkouvardos",
    author_email="pourjam.cs@hotmial.com, ilias.lagkouvardos@gmail.com",
    project_urls={
        "Repository": "https://github.com/MPourjam/TIClust",
    },
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "ticlust": ["bin/vsearch"],
    },
    cmdclass={
        "install": "ticlust.install_hooks.CustomInstall",
    },
)
