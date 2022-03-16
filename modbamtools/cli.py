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
@click.option('-r','--region', required=False, type=str, help='Region of interest. example: chr21:1-1000')
@click.option('-br','--batch', is_flag=False,default=None,required=False, type=click.Path(exists=True), help='makes html report for all regions in the bed file')
@click.option('-g','--gtf', multiple=True,is_flag=False,default=None,required=False, type=click.Path(exists=True), help='makes gene tracks from sorted and tabix gtf files')
@click.option('-b','--bed', multiple=True,is_flag=False,default=None,required=False, type=click.Path(exists=True), help='makes tracks from sorted and tabix bed files. This will plot each interval as a rectangle (similar to gtf)')
@click.option('-bw','--bigwig',multiple=True,is_flag=False,default=None, required=False, type=click.Path(exists=True), help='makes a track from bigwig files')
@click.option('-bd','--bedgraph',multiple=True,is_flag=False,default=None, required=False, type=click.Path(exists=True), help='makes a track from bigwig files')
@click.option('-s','--samples', is_flag=False, default=None, type=str, help='sample names per each bam input')
@click.option('-tr','--track-titles', is_flag=False, default=None, type=str, help='titles of tracks provided in order of gtf files, bed files, bigwig files, bedgraph files')
@click.option('-hp','--hap', is_flag=True, default=None, help='reads will be grouped according to HP tag in bam (comma separated)')
@click.option('-st','--strands', is_flag=True, default=None, help='reads will be grouped by strand in bam')
@click.option('-o','--out', required=True, type=click.Path(exists=True), help='output path for html plot')
@click.option('-p','--prefix', required=False, type=str,default= "modbamviz", help='File name for output')
@click.option('-f','--fmt', is_flag=False, default='html', type=str, help='format of output file (png, html, svg, pdf)')
@click.option('-u','--can_prob', is_flag=False, default=0.5, type=float, help='probability threshold for canonical bases')
@click.option('-m','--mod_prob', is_flag=False, default=0.5, type=float, help='probability threshold for modified bases')
@click.option('-h','--height', is_flag=False, default=None, type=int, help='height of plot in px. This is for fine tuning, the height is automatically calculated.')
@click.option('-w','--width', is_flag=False, default=None, type=int, help='width of plot in px')

def plot(bams,region,gtf, bed, bigwig, bedgraph, samples, hap, out, can_prob, mod_prob, height, width, prefix, fmt, strands, batch, track_titles):
    "This Command will plot single-read base modification data"
    if batch:
        html_start = '''
        @media print {
            .pagebreak { page-break-before: always; } /* page-break-after works, as well */
        }
        '''
        html_break='''
        <div style = "display:block; clear:both; page-break-after:always;"></div>
        '''
        out_path = out + "/" + prefix + "." + fmt
        figs = []
        if samples:
            samples = [s for s in samples.strip().split(',')]
        if track_titles:
            track_titles =  [t for t in track_titles.strip().split(',')]
        with open(batch,'r') as b:
            with open(out_path, 'w') as o:
                # if fmt =='html':
                #     o.write(html_break)
                for l in b:
                    if l[0] == '#':
                        continue
                    click.echo('processing '+l)
                    line = l.strip().split("\t")
                    chrom = line[0]
                    start = int(line[1])
                    end = int(line[2])
                    dicts, titles = get_reads(bams, chrom, start, end, hap=hap, strand=strands, samp_names= samples, min_prob=can_prob, max_prob=mod_prob)
                    plot = Plotter(dicts=dicts,samp_names=titles,chrom=chrom,start=start, end=end,
                    gtfs=gtf, beds=bed, bigwigs=bigwig, bedgraphs=bedgraph, track_titles=track_titles)
                    plot.plot_tracks()
                    if height:
                        plot.fig.update_layout(height=height)
                    if width:
                        plot.fig.update_layout(width=width)
                    if fmt == 'html':
                        o.write(plot.fig.to_html(full_html=False, include_plotlyjs='cdn'))
                        o.write(html_break)
                    if fmt != 'html':
                        figs.append(plot.fig)
        if fmt != 'html':
            click.echo('only html format for reports are supported currently')
            raise click.Abort()
            # pdf_report(figs, out_path)

    elif region:
        chrom = region.strip().split(':')[0]
        start = int(region.strip().split(':')[1].split('-')[0])
        end = int(region.strip().split(':')[1].split('-')[1])
        if samples:
            samples = [s for s in samples.strip().split(',')]
        if track_titles:
            track_titles =  [t for t in track_titles.strip().split(',')]
        dicts, titles = get_reads(bams, chrom, start, end, hap=hap, strand=strands, samp_names= samples, min_prob=can_prob, max_prob=mod_prob)

        fig = Plotter(dicts=dicts,samp_names=titles,chrom=chrom,start=start, end=end,
        gtfs=gtf, beds=bed, bigwigs=bigwig, bedgraphs=bedgraph,track_titles=track_titles)
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
    else:
        click.echo('Please choose either a region (--region) or bed file of regions (--batch) to process')

    click.echo('Successfully processed modified bases!')

