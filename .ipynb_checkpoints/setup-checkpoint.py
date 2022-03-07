from setuptools import setup
import os

VERSION = "0.1"


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
    packages=["modbamtools"],
    entry_points="""
        [console_scripts]
        modbamtools=modbamtools.cli:cli
    """,
    install_requires=["click"],
    extras_require={
        "test": ["pytest"]
    },
    python_requires=">=3.6",
)
