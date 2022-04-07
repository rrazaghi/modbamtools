# modbamtools commands

``` text
Usage: modbamtools [OPTIONS] COMMAND [ARGS]...

  A set of tools to manipulate and visualize data from base modification bam
  files

Options:
  --version  Show the version and exit.  [default: False]
  --help     Show this message and exit.  [default: False]

Commands:
  calcMeth  Calculate methylation statistics for regions in a bed file
  cluster   Calculate clustering statistics for regions in a bed file
  plot      Plot single-read base modification data
```

# 
## plot

``` text
Usage: modbamtools plot [OPTIONS] BAMS...

  Plot single-read base modification data

Options:
  -r, --region TEXT         Region of interest. example: chr21:1-1000
  -br, --batch PATH         makes html/pdf report for all regions in the bed
                            file
  -g, --gtf PATH            makes gene tracks from sorted and tabix gtf files
  -b, --bed PATH            makes tracks from sorted and tabix bed files. This
                            will plot each interval as a rectangle (similar to
                            gtf)
  -bw, --bigwig PATH        makes a track from bigwig files
  -bd, --bedgraph PATH      makes a track from bedgraph files
  -s, --samples TEXT        sample names per each bam input
  -tr, --track-titles TEXT  titles of tracks provided in order of gtf files,
                            bed files, bigwig files, bedgraph files
  -hp, --hap                reads will be grouped according to HP tag in bam
                            (comma separated)
  -st, --strands            reads will be grouped by strand in bam
  -o, --out PATH            output path  [required]
  -p, --prefix TEXT         File name for output  [default: modbamviz]
  -f, --fmt TEXT            format of output file (png, html, svg, pdf)
                            [default: html]
  -u, --can_prob FLOAT      probability threshold for canonical bases
                            [default: 0.5]
  -m, --mod_prob FLOAT      probability threshold for modified bases
                            [default: 0.5]
  -h, --height INTEGER      height of plot in px. This is for fine tuning, the
                            height is automatically calculated.
  -w, --width INTEGER       width of plot in px
  -c, --cluster             cluster the reads based on modification state
  --help                    Show this message and exit.
```
#
## calcMeth

```text
Usage: modbamtools calcMeth [OPTIONS] BAM

  Calculate methylation statistics for regions in a bed file

Options:
  -b, --bed PATH           bed file of regions
  -t, --threads INTEGER    number of processes  [default: 1]
  -a, --min_calls INTEGER  filter out reads that have fewer number of modified
                           base calls in region of interest  [default: 5]
  -s, --min_cov FLOAT      minimum percent coverage of a single read over
                           region of interest  [default: 80]
  -hp, --hap               add stats for each haplotype separately to the
                           output
  -o, --out PATH           output path  [required]
  --help                   Show this message and exit.
```

#
## cluster

```text
Usage: modbamtools cluster [OPTIONS] BAM

  Calculate clustering statistics for regions in a bed file

Options:
  -b, --bed PATH         bed file of regions to cluster
  -t, --threads INTEGER  number of processes  [default: 1]
  -o, --out PATH         output path  [required]
  --help                 Show this message and exit.
```
