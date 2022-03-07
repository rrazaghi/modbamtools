from setuptools import setup, find_packages
import os

VERSION = "0.1.0"


def get_long_description():
    with open(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "README.md"),
        encoding="utf8",
    ) as fp:
        return fp.read()

def read_requirements():
    with open('requirements.txt', 'r') as req:
        content = req.read()
        requirements = content.split('\n')

    return requirements


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
    install_requires=read_requirements(),
    extras_require={
        "test": ["pytest"]
    },
    python_requires=">=3.6, <=3.8",
)
