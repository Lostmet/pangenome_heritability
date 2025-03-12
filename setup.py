from setuptools import setup, find_packages

setup(
    name="pangenome_heritability",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "pandas>=1.3.0",
        "numpy>=1.20.0",
        "pysam>=0.16.0",
        "biopython>=1.79",
        "click>=8.0.0",
        "tqdm>=4.65.0",
    ],
    entry_points={
    'console_scripts': [
        'panherit=pangenome_heritability.cli:cli',
    ],

    },
    author="PeixiongYuan & JuntengWu",
    author_email="yuanpeixiong@westlake.edu.cn & fenglostmet@tju.edu.cn",
    description="A tool for refined structural variants",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/PeixiongYuan/pangenome_heritability",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
