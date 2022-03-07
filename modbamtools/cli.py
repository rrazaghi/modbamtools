import click
from modbamtools.modbamviz import *

@click.group()
@click.version_option()
def cli():
    "A set of tools to manipulate and visualize data from base modification bam files"
    pass


@cli.command(name="plot")
@click.argument("bams",nargs=-1,type=click.Path(exists=True),required=True)
@click.option('-r','--region', required=True, type=str, help='region of interest (chr21:1-1000)')
@click.option('-g','--gtf', required=True, type=click.Path(exists=True), help='gtf file in .gz sorted and tabix')
@click.option('-s','--samples', is_flag=False, default=None, type=str, help='sample names per each bam input')
@click.option('-h','--hap', is_flag=True, default=None, help='Reads will be grouped according to HP tag in bam')
@click.option('-o','--out', required=True, type=click.Path(), help='output path for html plot')
def plot(bams,region,gtf, samples, hap, out):
    "This Command will plot single-read base modification data"
    chrom = region.strip().split(':')[0]
    start = int(region.strip().split(':')[1].split('-')[0])
    end = int(region.strip().split(':')[1].split('-')[1])
    dicts, titles = get_reads(bams, chrom, start, end, hap=hap, samp_names= samples)

    fig = Plotter(dicts=dicts,samp_names=titles,gtf=gtf,chrom=chrom,start=start, end=end)
    fig.plot_genes()
    fig.plot()

    fig.fig.write_html(out)
    click.echo('Successfully processed modified bases!')

