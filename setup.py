from setuptools import setup, find_packages

setup(
    name='simple_tic',
    version='0.1',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    entry_points={
        'console_scripts': [
            'simpletic=simple_tic.cli:main',
        ],
    },
    install_requires=[
        # List your package dependencies here
    ],
)