![](docs/img/logo.png)

<!-- # modbamtools -->

[![PyPI](https://img.shields.io/pypi/v/modbamtools.svg)](https://pypi.org/project/modbamtools/)
[![Changelog](https://img.shields.io/github/v/release/rrazaghi/modbamtools?include_prereleases&label=changelog)](https://github.com/rrazaghi/modbamtools/releases)
[![Tests](https://github.com/rrazaghi/modbamtools/workflows/Test/badge.svg)](https://github.com/rrazaghi/modbamtools/actions?query=workflow%3ATest)
[![License](https://img.shields.io/badge/license-Apache%202.0-blue.svg)](https://github.com/rrazaghi/modbamtools/blob/master/LICENSE)

A set of tools to manipulate and visualize data from base modification bam files

For full documentation and tutorials visit https://rrazaghi.github.io/modbamtools/

## Installation

**<em>Required</em>**: Python 3.8

In a **clean** environment: 


    $ pip install modbamtools

## Usage

General commands:
```
Usage: modbamtools [OPTIONS] COMMAND [ARGS]...

  A set of tools to manipulate and visualize data from base modification bam
  files

Options:
  --version  Show the version and exit.  [default: False]
  --help     Show this message and exit.  [default: False]

Commands:
  calcHet   Calculate heterogeneity of modified bases for regions in a...
  calcMeth  Calculate methylation statistics for regions in a bed file
  cluster   Calculate clustering statistics for regions in a bed file
  plot      Plot single-read base modification data


```
![example plot](./tests/modbamtools.gif)



## Development

To contribute to this tool, first checkout the code. Then create a new virtual environment:

    cd modbamtools
    python -m venv venv
    source venv/bin/activate

Or if you are using `pipenv`:

    pipenv shell

Now install the dependencies and test dependencies:

    pip install -e '.[test]'

To run the tests:

    pytest
