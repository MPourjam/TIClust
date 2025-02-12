from setuptools import setup, find_packages
import requests
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent / 'src/tic'))
from install_hooks import CustomInstall


setup(
    name="tic",
    version="1.0.0",
    description="Taxonomy Informed Clustering (TIC) is a tool for clustering bacterial sequences based on their taxonomy and hypothetically complete taxonomy levels.",
    author="Mohsen Pourjam, Ilias Lagkouvardos",
    author_email="pourjam.cs@hotmial.com, ilias.lagkouvardos@gmail.com",
    # url="https://example.com",
    project_urls={
        "Repository": "https://github.com/yourusername/yourproject",
    },
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    license="MIT",
    python_requires=">=3.6",
    install_requires=[
        "requests",
        "toml",
    ],
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    package_data={
        "tic": ["bin/vsearch"],
    },
    entry_points={
        "console_scripts": [
            "simpletic=tic.cli:main",
        ],
    },
    cmdclass={
        "install": "tic.install_hooks.CustomInstall",
    },
)
