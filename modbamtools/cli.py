from email.policy import default
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
@click.option('-hp','--hap', is_flag=True, default=None, help='Reads will be grouped according to HP tag in bam')
@click.option('-o','--out', required=True, type=click.Path(exists=True), help='output path for html plot')
@click.option('-p','--prefix', required=False, type=str,default= "modbamviz", help='file name for output')
@click.option('-f','--fmt', is_flag=False, default='html', type=str, help='format of output file (png, html, svg, pdf)')
@click.option('-u','--can_prob', is_flag=False, default=0.5, type=float, help='probability threshold for canonical bases')
@click.option('-m','--mod_prob', is_flag=False, default=0.5, type=float, help='probability threshold for modified bases')
@click.option('-h','--height', is_flag=False, default=2000, type=int, help='height of plot in px')
@click.option('-w','--width', is_flag=False, default=1500, type=int, help='width of plot in px')





def plot(bams,region,gtf, samples, hap, out, can_prob, mod_prob, height, width, prefix, fmt):
    "This Command will plot single-read base modification data"
    chrom = region.strip().split(':')[0]
    start = int(region.strip().split(':')[1].split('-')[0])
    end = int(region.strip().split(':')[1].split('-')[1])
    dicts, titles = get_reads(bams, chrom, start, end, hap=hap, samp_names= samples, min_prob=can_prob, max_prob=mod_prob)

    fig = Plotter(dicts=dicts,samp_names=titles,gtf=gtf,chrom=chrom,start=start, end=end)
    fig.plot_genes()
    fig.plot()
    fig.fig.update_layout(height=height, width=width)

    out_path = out + "/" + prefix + "." + fmt

    if fmt == "html":
        fig.fig.write_html(out_path)
    else:
        fig.fig.write_image(out_path)

    click.echo('Successfully processed modified bases!')

