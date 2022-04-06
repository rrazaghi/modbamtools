# Modbamtools Documentation

modbamtools is a set of tools to manipulate and visualize DNA/RNA base modification data that are stored in bam format. htslib has included a support for parsing modified base tags from alignment files (MM ans ML). For more information about these tags, please visit http://samtools.github.io/hts-specs/SAMtags.pdf

## Generate modified base tags for your data
modbamtools is technology agnostic. However, tools tailored for analysis of modified bases using long-read technology are currently adapting to using MM/ML tags at a much higher rate. below are the list of tools that can generate these tags to be used with modbamtools:

Oxford Nanopore Technology (ONT)  
Bonito 
Guppy  
Nanopolish  
Megalodon  
Remora

Pacific Biosciences (Pacbio)  
[Primrose](https://github.com/PacificBiosciences/primrose)


## Commands

* `mkdocs new [dir-name]` - Create a new project.
* `mkdocs serve` - Start the live-reloading docs server.
* `mkdocs build` - Build the documentation site.
* `mkdocs -h` - Print help message and exit.

## Project layout

    mkdocs.yml    # The configuration file.
    docs/
        index.md  # The documentation homepage.
        ...       # Other markdown pages, images and other files.
