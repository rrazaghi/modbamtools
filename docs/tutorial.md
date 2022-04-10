# Tutorial
For this tutorial we are going to use ultra-long nanopore data from GM12878 aligned to chromosome 20 only. The modbam was generated using [megalodon](https://github.com/nanoporetech/megalodon) and haplotype tags were assigned using [PEPPER](https://github.com/kishwarshafin/pepper).

## Download Data

You can download all the files needed to go through this tutorial [here](https://timpshare.s3.amazonaws.com/modbamtools_tutorial_files.tar.gz).

After extractiong, you should have the following files:
```
.
├── gencode.v38.annotation.sorted.gtf.gz
├── gencode.v38.annotation.sorted.gtf.gz.tbi
├── genes.bed
├── gm12878_H3K27ac_ENCFF798KYP.bigWig
├── gm12878_H3K4me1_ENCFF190RZM.bigWig
├── gm12878_ul_sup_megalodon_HP_chr20.bam
├── gm12878_ul_sup_megalodon_HP_chr20.bam.bai
└── promoters.bed
```

## Plotting

### Basic plot

First let's start with plotting an imprinted region for GNAS gene:

```bash
modbamtools plot -r chr20:58815000-58895000 \
    --gtf gencode.v38.annotation.sorted.gtf.gz \
    --out . \
    --prefix gm12878_GNAS \
    --samples GM12878 \
    --track-titles Genes\
    gm12878_ul_sup_megalodon_HP_chr20.bam 
```
<iframe width=900, height=1050 frameBorder=0 src="../figs/gm12878_GNAS.html"></iframe>

[gm12878_GNAS](figs/gm12878_GNAS.html)
### Haplotype separation
We can group the reads based on `HP` tag by using `--hap` option:
```bash
modbamtools plot -r chr20:58815000-58895000 \
    --gtf gencode.v38.annotation.sorted.gtf.gz \
    --out . \
    --hap \
    --prefix gm12878_GNAS \
    --samples GM12878 \
    --track-titles Genes\
    gm12878_ul_sup_megalodon_HP_chr20.bam 
```
<iframe width=900, height=1000 frameBorder=0 src="../figs/gm12878_GNAS_hap.html"></iframe>

[gm12878_GNAS_hap](figs/gm12878_GNAS_hap.html)
### Bigwig tracks
Bigwig tracks can also be added to the plot with `--bigwig`:
```bash
modbamtools plot -r chr20:58820000-58895000 \
    --gtf gencode.v38.annotation.sorted.gtf.gz \
    --out . \
    --hap \
    --bigwig gm12878_H3K27ac_ENCFF798KYP.bigWig \
    --bigwig gm12878_H3K4me1_ENCFF190RZM.bigWig \
    --prefix gm12878_GNAS_hap_h3k27ac_h3k4me1 \
    --samples GM12878 \
    --track-titles Genes,H3K27ac,H3k4me1\
    gm12878_ul_sup_megalodon_HP_chr20.bam 
```
<iframe width=900, height=1160 frameBorder=0 src="../figs/gm12878_GNAS_hap_h3k27ac_h3k4me1.html"></iframe>

[gm12878_GNAS_hap_h3k27ac_h3k4me1](figs/gm12878_GNAS_hap_h3k27ac_h3k4me1.html)

### Cluster plot
We can use `--cluster` option to perform clustering on the go and group reads based on the assigned cluster. To show an example, we focus on gene FRG1EP (you can find unclustered plot [here](figs/gm12878_FRG1EP.html)). For more information about the clustering algorithm see section clustering below.

```bash
modbamtools plot -r chr20:29460146-29507179 \
    --gtf gencode.v38.annotation.sorted.gtf.gz \
    --out . \
    --cluster \
    --prefix gm12878_FRG1EP_cluster \
    --track-titles Genes\
    gm12878_ul_sup_megalodon_HP_chr20.bam
```
<iframe width=900, height=770 frameBorder=0 src="../figs/gm12878_FRG1EP_cluster.html"></iframe>

[gm12878_FRG1EP_cluster](figs/gm12878_FRG1EP_cluster.html)

### Multiple input regions (Bed)
If we are interested in plotting multiple regions simultaneously, we can feed modbamtools a bed file. Output formats can either be `html` or `pdf` in a single report document.

```bash
modbamtools plot --batch genes.bed \
    --gtf gencode.v38.annotation.sorted.gtf.gz \
    --out . \
    --hap \
    --prefix gm12878_batch \
    --samples GM12878 \
    --track-titles Genes\
    gm12878_ul_sup_megalodon_HP_chr20.bam
```
[gm12878_batch](figs/gm12878_batch.html)

## Clustering
modbamtools uses a hierarichal density-based spatial clustering of noisy data algorithm ([HDBSCAN](https://github.com/scikit-learn-contrib/hdbscan)) to cluster reads based on their modification state. `modbamtools cluster` takes a bed file as input and appends each line with two columns : `number of clusters found`, and `number of reads used for clustering after a quality control step`

<div class="alert alert-danger" role="alert">
  It should be noted that modification calls are often noisy and sparse. Therefore extra care should be applied when clustering. Best way is to plot clustering results using `modbamtools plot --cluster` to make sure clustering algorithm fits the nature of your data
</div>

Let's try it on our `genes.bed` file:

```bash
modbamtools cluster --bed genes.bed \
    --threads 3 \
    --out genes_clustered.bed \
    gm12878_ul_sup_megalodon_HP_chr20.bam
```
We can see from the output below that all 3 regions have 2 clusters. We can confirm this by using the plotting commands above.

```
chr20	29460146	29507179	FRG1EP	unprocessed_pseudogene	ENSG00000282995.1	2	42
chr20	58818918	58850903	GNAS-AS1	lncRNA	ENSG00000235590.7	2	30
chr20	63556085	63576239	HELZ2	protein_coding	ENSG00000130589.16	2	37
```
<div class="alert alert-info" role="alert">
  This command is especially useful when trying to distinguish reads associated with tumors in a heterogenous sample. You can also feed multiple modbams to find differentially clustered reads between different samples!
</div>

## Regional methylation calculation

modbamtools `calcMeth` will accept a bed file as input and append modification stats as extra columns to each region of interest. For each region, Average percent modification is calculated for each read that spans a minimum portion of the region (This can be tweaked by `--min_cov`. Default is 80%). Each read has to have at least 5 calls to be considered valid for analysis (This can also be changed with `--min_calls`). Then, Average percent modifcation for the region is calculated using valid reads. The appended columns are outlined below:

<table class="table table-bordered">
  <thead>
	  <tr>
	  	<th>Average Percent Modification</th>
	  	<th>Standrad Deviation for Avg Percent Modification</th>
	  	<th>Coverage</th>
	  </tr>
	  </thead>
</table>

Let's try it on our `promoters.bed` file:

```bash
modbamtools calcMeth --bed promoters.bed \
    --threads 3 \
    --out promoters_calcMeth.bed \
    gm12878_ul_sup_megalodon_HP_chr20.bam
```
```
chr20	29496679	29498179	FRG1EP	unprocessed_pseudogene	ENSG00000282995.1	63.00634042910756	23.160648801539175	79
chr20	58850403	58851903	GNAS-AS1	lncRNA	ENSG00000235590.7	58.26586037568643	30.18888731341584	48
chr20	63573739	63575239	HELZ2	protein_coding	ENSG00000130589.16	36.52218479283339	8.51755978542393	44
```
### Haplotype stats
We can also output stats for each haplotype based on `HP` tag in the modbam with `--hap`. The added columns when this option is used are the following:

<table class="table table-bordered">
  <thead>
	  <tr>
	  	<th>Average Percent Modification</th>
	  	<th>Standrad Deviation for Avg Percent Modification</th>
	  	<th>Coverage</th>
        <th>Average Percent Modification (haplotype 1)</th>
	  	<th>Standrad Deviation for Avg Percent Modification (haplotype 1)</th>
	  	<th>Coverage (haplotype 1)</th>
        <th>Average Percent Modification (haplotype 2)</th>
	  	<th>Standrad Deviation for Avg Percent Modification (haplotype 2)</th>
	  	<th>Coverage (haplotype 2)</th>
	  </tr>
	  </thead>
</table>

```bash
modbamtools calcMeth --bed promoters.bed \
    --threads 3 \
    --hap \
    --out promoters_calcMeth_hap.bed \
    gm12878_ul_sup_megalodon_HP_chr20.bam
```

```
chr20	29496679	29498179	FRG1EP	unprocessed_pseudogene	ENSG00000282995.1	63.00634042910756	23.160648801539175	79	44.27765063921806	18.419436878488273	36	80.732273703281	10.455159906458562	36
chr20	58850403	58851903	GNAS-AS1	lncRNA	ENSG00000235590.7	58.26586037568643	30.18888731341584	48	26.62614126968194	6.853899031338103	22	86.90174211774595	5.326307202363974	20
chr20	63573739	63575239	HELZ2	protein_coding	ENSG00000130589.16	36.52218479283339	8.51755978542393	44	38.93859253773206	8.05201039505814	16	35.09734468342445	8.949332073699683	25
```
