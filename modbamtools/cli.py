import click
from modbamtools.modbamviz import *
from modbamtools.calcregions import *
from modbamtools.clustering import *
from modbamtools.heterogeneity import *
import io
from PIL import Image
from PyPDF2 import PdfFileMerger
from itertools import repeat
from multiprocessing import cpu_count, get_context, Pool
import tqdm


@click.group(context_settings={"show_default": True})
@click.version_option()
def cli():
    "A set of tools to manipulate and visualize data from base modification bam files"
    pass


def listify(ctx, param, value):
    if type(value) == str:
        return [value]


@cli.command(name="plot")
@click.argument("bams", nargs=-1, type=click.Path(exists=True), required=True)
@click.option(
    "-r",
    "--region",
    required=False,
    type=str,
    help="Region of interest. example: chr21:1-1000",
)
@click.option(
    "-br",
    "--batch",
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="makes html/pdf report for all regions in the bed file",
)
@click.option(
    "-g",
    "--gtf",
    multiple=True,
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="makes gene tracks from sorted and tabix gtf files",
)
@click.option(
    "-b",
    "--bed",
    multiple=True,
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="makes tracks from sorted and tabix bed files. This will plot each interval as a rectangle (similar to gtf)",
)
@click.option(
    "-bw",
    "--bigwig",
    multiple=True,
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="makes a track from bigwig files",
)
@click.option(
    "-bd",
    "--bedgraph",
    multiple=True,
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="makes a track from bedgraph files",
)
@click.option(
    "-s",
    "--samples",
    is_flag=False,
    default=None,
    type=str,
    help="sample names per each bam input",
)
@click.option(
    "-tr",
    "--track-titles",
    is_flag=False,
    default=None,
    type=str,
    help="titles of tracks provided in order of gtf files, bed files, bigwig files, bedgraph files",
)
@click.option(
    "-hp",
    "--hap",
    is_flag=True,
    default=None,
    help="reads will be grouped according to HP tag in bam (comma separated)",
)
@click.option(
    "-st",
    "--strands",
    is_flag=True,
    default=None,
    help="reads will be grouped by strand in bam",
)
@click.option(
    "-o",
    "--out",
    required=True,
    type=click.Path(exists=True),
    help="output path",
)
@click.option(
    "-p",
    "--prefix",
    required=False,
    type=str,
    default="modbamviz",
    help="File name for output",
)
@click.option(
    "-f",
    "--fmt",
    is_flag=False,
    default="html",
    type=str,
    help="format of output file (png, html, svg, pdf)",
)
@click.option(
    "-u",
    "--can_prob",
    is_flag=False,
    default=0.5,
    type=float,
    help="probability threshold for canonical bases",
)
@click.option(
    "-m",
    "--mod_prob",
    is_flag=False,
    default=0.5,
    type=float,
    help="probability threshold for modified bases",
)
@click.option(
    "-h",
    "--height",
    is_flag=False,
    default=None,
    type=int,
    help="height of plot in px. This is for fine tuning, the height is automatically calculated.",
)
@click.option(
    "-w", "--width", is_flag=False, default=None, type=int, help="width of plot in px"
)
@click.option(
    "-c",
    "--cluster",
    is_flag=True,
    default=None,
    type=int,
    help="cluster the reads based on modification state",
)
@click.option(
    "-ht",
    "--heterogeneity",
    is_flag=True,
    default=None,
    type=int,
    help="plot degree of modification heterogeneity across the region",
)
@click.option(
    "-fs",
    "--font_size",
    is_flag=False,
    default=18,
    type=int,
    help="global font size",
)
@click.option(
    "-ms",
    "--marker_size",
    is_flag=False,
    default=6,
    type=float,
    help="marker size of each modified/unmodified base",
)
@click.option(
    "-sth",
    "--single_trace_height",
    is_flag=False,
    default=12,
    type=int,
    help="space between single molucles in px",
)
def plot(
    bams,
    region,
    gtf,
    bed,
    bigwig,
    bedgraph,
    samples,
    hap,
    out,
    can_prob,
    mod_prob,
    height,
    width,
    prefix,
    fmt,
    strands,
    batch,
    track_titles,
    cluster,
    heterogeneity,
    font_size,
    marker_size,
    single_trace_height,
):
    "Plot single-read base modification data"
    if batch:
        if (fmt != "html") & (fmt != "pdf"):
            click.echo("only html and pdf formats for reports are currently supported!")
            raise click.Abort()
        html_break = """
        <div style = "display:block; clear:both; page-break-after:always;"></div>
        """
        out_path = out + "/" + prefix + "." + fmt
        figs = []
        if samples:
            samples = [s for s in samples.strip().split(",")]
        if track_titles:
            track_titles = [t for t in track_titles.strip().split(",")]
        if fmt == "pdf":
            merger = PdfFileMerger()

        with open(batch, "r") as b:
            for l in b:
                if l[0] == "#":
                    continue
                click.echo("processing " + l)
                line = l.strip().split("\t")
                chrom = line[0]
                start = int(line[1])
                end = int(line[2])
                if cluster:
                    dicts, titles = cluster2dicts(bams, chrom, start, end)
                if not cluster:
                    dicts, titles = get_reads(
                        bams,
                        chrom,
                        start,
                        end,
                        hap=hap,
                        strand=strands,
                        samp_names=samples,
                        min_prob=can_prob,
                        max_prob=mod_prob,
                    )
                plot = Plotter(
                    dicts=dicts,
                    samp_names=titles,
                    chrom=chrom,
                    start=start,
                    end=end,
                    gtfs=gtf,
                    beds=bed,
                    bigwigs=bigwig,
                    bedgraphs=bedgraph,
                    track_titles=track_titles,
                    heterogeneity=heterogeneity,
                    font_size=font_size,
                    fmt=fmt,
                    marker_size=marker_size,
                    single_trace_height=single_trace_height,
                )
                plot.plot_tracks()
                if height:
                    plot.fig.update_layout(height=height)
                if width:
                    plot.fig.update_layout(width=width)
                if fmt == "html":
                    plot.fig.update_layout(title=chrom + ":" + line[1] + "-" + line[2])
                    with open(out_path, "a") as o:
                        o.write(plot.fig.to_html(full_html=True))
                        o.write(html_break)
                # if fmt == 'png':
                #     plot.fig.update_layout(title=chrom + ':' + line[1] + '-' + line[2])
                #     img_bytes = plot.fig.to_image(format="png")
                #     im = Image.open(io.BytesIO(img_bytes)).convert('RGB')
                #     figs.append(im)
                if fmt == "pdf":
                    plot.fig.update_layout(title=chrom + ":" + line[1] + "-" + line[2])
                    plot.fig.update_layout(width=1200)
                    img_bytes = plot.fig.to_image(format="pdf")
                    merger.append(io.BytesIO(img_bytes))

        # if fmt == 'png':
        #     figs.pop(0).save(out_path,'PDF', save_all=True, append_images=figs, resolution=100.0)
        if fmt == "pdf":
            merger.write(out_path)

    elif region:
        chrom = region.strip().split(":")[0]
        start = int(region.strip().split(":")[1].split("-")[0])
        end = int(region.strip().split(":")[1].split("-")[1])
        if samples:
            samples = [s for s in samples.strip().split(",")]
        if track_titles:
            track_titles = [t for t in track_titles.strip().split(",")]
        if cluster:
            dicts, titles = cluster2dicts(bams, chrom, start, end)
        if not cluster:
            dicts, titles = get_reads(
                bams,
                chrom,
                start,
                end,
                hap=hap,
                strand=strands,
                samp_names=samples,
                min_prob=can_prob,
                max_prob=mod_prob,
            )
        fig = Plotter(
            dicts=dicts,
            samp_names=titles,
            chrom=chrom,
            start=start,
            end=end,
            gtfs=gtf,
            beds=bed,
            bigwigs=bigwig,
            bedgraphs=bedgraph,
            track_titles=track_titles,
            heterogeneity=heterogeneity,
            font_size=font_size,
            fmt=fmt,
            marker_size=marker_size,
            single_trace_height=single_trace_height,
        )
        fig.plot_tracks()
        if height:
            fig.fig.update_layout(height=height)
        if width:
            fig.fig.update_layout(width=width)

        out_path = out + "/" + prefix + "." + fmt

        if fmt == "html":
            fig.fig.write_html(out_path)
            # print(fig.fig.layout.height)
        else:
            fig.fig.write_image(out_path)
            # fig.fig.write_image(out_path, width=5 * 96, height=fig.plot_height, scale=1)
            # print(fig.fig.layout.height)
    else:
        click.echo(
            "Please choose either a region (--region) or bed file of regions (--batch) to process"
        )

    click.echo("Successfully processed modified bases! ")


@cli.command(name="calcMeth")
@click.argument("bam", nargs=1, type=click.Path(exists=True), required=True)
@click.option(
    "-b",
    "--bed",
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="bed file of regions",
)
@click.option(
    "-t", "--threads", is_flag=False, default=1, type=int, help="number of processes"
)
@click.option(
    "-a",
    "--min_calls",
    is_flag=False,
    default=5,
    type=int,
    help="filter out reads that have fewer number of modified base calls in region of interest",
)
@click.option(
    "-s",
    "--min_cov",
    is_flag=False,
    default=80,
    type=float,
    help="minimum percent coverage of a single read over region of interest",
)
@click.option(
    "-hp",
    "--hap",
    is_flag=True,
    default=None,
    help="add stats for each haplotype separately to the output",
)
@click.option(
    "-o",
    "--out",
    required=True,
    type=click.Path(),
    help="output path",
)
def calcMeth(bam, bed, min_calls, min_cov, threads, hap, out):
    "Calculate methylation statistics for regions in a bed file"
    data = []
    with open(bed, "r") as f:
        for line in f:

            bed_line = line.strip().split("\t")
            if (bed_line[0] == "chr") | (line[0] == "#"):
                continue
            data.append(bed_line)

    results = []
    with Pool(threads) as p:
        s = p.starmap(
            get_region_stats_hap_horizontal,
            zip(repeat(bam), repeat(min_calls), repeat(min_cov), repeat(hap), data),
        )
        results.extend(s)
        p.close()
        p.join()

    with open(out, "w") as o:
        for rec in results:
            print("\t".join(rec), end="\n", file=o)


@cli.command(name="cluster")
@click.argument("bam", nargs=1, type=click.Path(exists=True), required=True)
@click.option(
    "-b",
    "--bed",
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="bed file of regions to cluster",
)
@click.option(
    "-t", "--threads", is_flag=False, default=1, type=int, help="number of processes"
)
@click.option(
    "-o",
    "--out",
    required=True,
    type=click.Path(),
    help="output path",
)
def cluster(bam, bed, threads, out):
    "Calculate clustering statistics for regions in a bed file"

    data = []
    with open(bed, "r") as f:
        for line in f:

            bed_line = line.strip().split("\t")
            if (bed_line[0] == "chr") | (line[0] == "#"):
                continue
            data.append(bed_line)

    cluster_region_y = partial(cluster_region, bam=bam)
    #     results = p.imap(cluster_region_y, data)
    with Pool(threads) as p:
        results = list(tqdm.tqdm(p.imap(cluster_region_y, data), total=len(data)))
        #     results = p.starmap(cluster_region, zip(repeat(bam), data))
        p.close()
        p.join()

    with open(out, "w") as o:
        for rec in results:
            print("\t".join(rec), end="\n", file=o)


@cli.command(name="calcHet")
@click.argument("bam", nargs=1, type=click.Path(exists=True), required=True)
@click.option(
    "-b",
    "--bed",
    is_flag=False,
    default=None,
    required=False,
    type=click.Path(exists=True),
    help="bed file of regions",
)
@click.option(
    "-t", "--threads", is_flag=False, default=1, type=int, help="number of processes"
)
@click.option(
    "-a",
    "--min_calls",
    is_flag=False,
    default=5,
    type=int,
    help="filter out reads that have fewer number of modified base calls in region of interest",
)
@click.option(
    "-s",
    "--min_cov",
    is_flag=False,
    default=80,
    type=float,
    help="minimum percent coverage of a single read over region of interest",
)
@click.option(
    "-hp",
    "--hap",
    is_flag=True,
    default=None,
    help="add stats for each haplotype separately to the output",
)
@click.option(
    "-o",
    "--out",
    required=True,
    type=click.Path(),
    help="output path",
)
def calcHet(bam, bed, min_calls, min_cov, threads, hap, out):
    "Calculate heterogeneity of modified bases for regions in a bed file"
    data = []
    with open(bed, "r") as f:
        for line in f:

            bed_line = line.strip().split("\t")
            if (bed_line[0] == "chr") | (line[0] == "#"):
                continue
            data.append(bed_line)

    results = []
    with Pool(threads) as p:
        s = p.starmap(
            get_region_hap_heterogeneity,
            zip(repeat(bam), repeat(min_calls), repeat(min_cov), repeat(hap), data),
        )
        results.extend(s)
        p.close()
        p.join()

    with open(out, "w") as o:
        for rec in results:
            print("\t".join(rec), end="\n", file=o)
