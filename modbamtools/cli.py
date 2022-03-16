import click
from modbamtools.modbamviz import *

@click.group()
@click.version_option()
def cli():
    "A set of tools to manipulate and visualize data from base modification bam files"
    pass

def listify(ctx, param, value):
    if type(value) == str:
        return [value]

@cli.command(name="plot")
@click.argument("bams",nargs=-1,type=click.Path(exists=True),required=True)
@click.option('-r','--region', required=True, type=str, help='Region of interest. example: chr21:1-1000')
@click.option('-g','--gtf', callback=listify,is_flag=False,default=None,required=False, type=click.Path(exists=True), help='make gene tracks from sorted and tabix gtf files')
@click.option('-b','--bed', callback=listify,is_flag=False,default=None,required=False, type=click.Path(exists=True), help='make tracks from sorted and tabix bed files. This will plot each interval as a rectangle (similar to gtf)')
@click.option('-bw','--bigwig',callback=listify,is_flag=False,default=None, required=False, type=click.Path(exists=True), help='make a track from bigwig files')
@click.option('-bd','--bedgraph',callback=listify,is_flag=False,default=None, required=False, type=click.Path(exists=True), help='make a track from bigwig files')
@click.option('-s','--samples', is_flag=False, default=None, type=str, help='Sample names per each bam input')
@click.option('-hp','--hap', is_flag=True, default=None, help='Reads will be grouped according to HP tag in bam')
@click.option('-st','--strands', is_flag=True, default=None, help='Reads will be grouped by strand in bam')
@click.option('-o','--out', required=True, type=click.Path(exists=True), help='Output path for html plot')
@click.option('-p','--prefix', required=False, type=str,default= "modbamviz", help='File name for output')
@click.option('-f','--fmt', is_flag=False, default='html', type=str, help='Format of output file (png, html, svg, pdf)')
@click.option('-u','--can_prob', is_flag=False, default=0.5, type=float, help='Probability threshold for canonical bases')
@click.option('-m','--mod_prob', is_flag=False, default=0.5, type=float, help='Probability threshold for modified bases')
@click.option('-h','--height', is_flag=False, default=None, type=int, help='Height of plot in px')
@click.option('-w','--width', is_flag=False, default=None, type=int, help='Width of plot in px')

def plot(bams,region,gtf, bed, bigwig, bedgraph, samples, hap, out, can_prob, mod_prob, height, width, prefix, fmt, strands):
    "This Command will plot single-read base modification data"
    chrom = region.strip().split(':')[0]
    start = int(region.strip().split(':')[1].split('-')[0])
    end = int(region.strip().split(':')[1].split('-')[1])
    if samples:
        samples = [s for s in samples.strip().split(',')]
    dicts, titles = get_reads(bams, chrom, start, end, hap=hap, strand=strands, samp_names= samples, min_prob=can_prob, max_prob=mod_prob)

    fig = Plotter(dicts=dicts,samp_names=titles,chrom=chrom,start=start, end=end,
    gtfs=gtf, beds=bed, bigwigs=bigwig, bedgraphs=bedgraph)
    fig.plot_tracks()
    if height:
        fig.fig.update_layout(height=height)
    if width:
        fig.fig.update_layout(width=width)

    out_path = out + "/" + prefix + "." + fmt

    if fmt == "html":
        fig.fig.write_html(out_path)
    else:
        fig.fig.write_image(out_path)

    click.echo('Successfully processed modified bases!')

