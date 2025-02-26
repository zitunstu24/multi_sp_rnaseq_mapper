from setuptools import setup, find_packages, Extension


setup(
    name="multi-sp-rnaseq-mapper",
    version="0.1",
    description="SRA sample downloader and STAR alignment tool for multi species",
    packages=find_packages(),
    install_requires=[
        'pandas==2.2.3',    # Specify exact version of pandas
        'pyyaml==5.4.1',
        'pydeseq2'
    ],
    entry_points={
        'console_scripts': [
            'RnaMapper=multi_sp_rnaseq_mapper.main:main'
        ]
    },
    python_requires='>=3.8',  # Specify Python version requirement
)
