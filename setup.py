from setuptools import setup, find_packages
import os
import io

VERSION = "0.4.8"


def get_long_description():
    with open(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md"),
        encoding="utf8",
    ) as fp:
        return fp.read()


setup(
    name="modbamtools",
    description="A set of tools to manipulate and visualize data from base modification bam files",
    long_description=get_long_description(),
    long_description_content_type="text/markdown",
    author="Roham Razaghi",
    url="https://github.com/rrazaghi/modbamtools",
    project_urls={
        "Issues": "https://github.com/rrazaghi/modbamtools/issues",
        "CI": "https://github.com/rrazaghi/modbamtools/actions",
        "Changelog": "https://github.com/rrazaghi/modbamtools/releases",
    },
    license="Apache License, Version 2.0",
    version=VERSION,
    packages=find_packages(),
    entry_points="""
        [console_scripts]
        modbamtools=modbamtools.cli:cli
    """,
    install_requires=[
        "click>=8.0.4",
        "pysam>=0.18.0",
        "scipy>=1.7.0",
        "pandas>=1.4.0",
        "numpy>=1.22.0",
        "plotly>=5.5.0",
        "modbampy==0.5.3",
        "kaleido>=0.2.1",
        "pyBigWig>=0.3.18",
        "PyPDF2",
        "pillow",
        "hdbscan",
        "tqdm",
    ],
    extras_require={"test": ["pytest"]},
    python_requires=">=3.8, <=3.10.5",
)
