from setuptools import setup, find_packages
import requests
import os
import sys
from pathlib import Path
sys.path.append(str(Path(__file__).resolve().parent / 'src/tic'))
from install_hooks import CustomInstall


setup(
    name='tic',
    version='0.2',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    cmdclass={'install': CustomInstall},
    include_package_data=True,
    package_data={
        'tic': ['bin/vsearch'],
    },
    entry_points={
        'console_scripts': [
            'simpletic=tic.cli:main',
        ],
    },
    install_requires=[
        'requests',
        # 'toml',
    ],
)