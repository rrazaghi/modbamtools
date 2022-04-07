# Modbamtools Documentation

modbamtools is a set of tools to manipulate and visualize DNA/RNA base modification data that are stored in bam format. htslib has included a support for parsing modified base tags from alignment files (MM and ML). For more information about these tags, please visit [here](http://samtools.github.io/hts-specs/SAMtags.pdf).

## Generate modified base tags for your data
modbamtools is technology agnostic. However, tools tailored for analysis of modified bases using long-read technology are currently adapting to using MM/ML tags at a much higher rate. Below are the list of tools that can generate these tags to be used with modbamtools:

**Oxford Nanopore Technology (ONT)**  
[Bonito](https://github.com/nanoporetech/bonito)  
[Guppy](https://community.nanoporetech.com/downloads)  
[Nanopolish](https://github.com/jts/nanopolish)  
[Megalodon](https://github.com/nanoporetech/megalodon)  
[Remora](https://github.com/nanoporetech/remora)

**Pacific Biosciences (Pacbio)**  
[Primrose](https://github.com)

## Install

**<em>Required</em>**: Python 3.8

In a clean environment: 

<pre><code class="shell">$ pip install modbamtools</code></pre>





## Commands

* `modbamtools plot` - Plot single-read base modification data along other optional tracks (gtf, bedgraph, bigwig).
* `modbamtools calcMeth` - Calculate methylation statistics for regions in a bed file.
* `modbamtools cluster` - Perform clustering based on modification state for regions in a bed file.
* `modbamtools --help` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
