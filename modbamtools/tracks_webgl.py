from modbamtools.utils import *
from modbamtools.gene_models import *
from modbamtools.heterogeneity import get_dict_heterogeneity
import pyBigWig


def parse_bigwig_gl(bigwig_path, chrom, start, end):

    bw = pyBigWig.open(bigwig_path)
    x = list(range(start, end))
    y = bw.values(chrom, start, end)
    trace = go.Scattergl(
        x=x,
        y=y,
        mode="lines",
        marker=dict(
            size=6,
            color="cadetblue",
        ),
        name="",
        showlegend=False,
    )
    height = 100  # px
    return trace, height


def parse_bedgraph_gl(bedgraph_path, chrom, start, end):
    x = []
    y = []
    with open(bedgraph_path) as b:
        for l in b:
            if l[0] == "#":
                continue
            line = l.strip().split("\t")
            ch = line[0]
            st = int(line[1])
            nd = int(line[2])
            score = float(line[3])
            if ch == chrom:
                if st >= start:
                    if nd <= end:
                        if (st == nd) | (st + 1 == nd):
                            x.append(st)
                            y.append(score)
                        else:
                            x.append(st)
                            x.append(nd)
                            y.append(score)
                            y.append(score)

    trace = go.Scattergl(
        x=x,
        y=y,
        mode="lines",
        marker=dict(
            size=6,
            color="goldenrod",
        ),
        name="",
        showlegend=False,
    )
    height = 100  # px

    return trace, height


def plot_frequencies_gl(dict_per_read_mod, start, end, color):
    height = 200  # px

    traces = []
    freq, freq_smooth = calc_freq(dict_per_read_mod, start, end)

    traces.append(
        go.Scattergl(
            x=freq["x"],
            y=freq["y"],
            mode="markers",
            marker=dict(
                size=3,
                color=color,
            ),
            name="",
            showlegend=False,
        )
    )

    traces.append(
        go.Scattergl(
            x=freq_smooth["x"],
            y=freq_smooth["y"],
            mode="lines",
            marker=dict(
                size=3,
                color=color,
            ),
            name="",
            showlegend=False,
        )
    )

    return traces, height


def plot_freq_diff_gl(
    dict_per_read_mod_hp1, dict_per_read_mod_hp2, start, end, color="grey"
):

    freq_hp1, freq_smooth_hp1 = calc_freq(dict_per_read_mod_hp1, start, end)
    freq_hp2, freq_smooth_hp2 = calc_freq(dict_per_read_mod_hp2, start, end)

    return go.Scattergl(
        x=freq_smooth_hp2["x"],
        y=[a_i - b_i for a_i, b_i in zip(freq_smooth_hp2["y"], freq_smooth_hp1["y"])],
        mode="lines",
        marker=dict(
            size=6,
            color=color,
        ),
        name="",
    )


def parse_bed_gl(bed_path, chrom, start, end):

    with open(bed_path) as b:
        for l in b:
            if l[0] == "#":
                continue
            line = l.strip().split("\t")
            if line[0] == chrom:
                if int(line[1]) >= start:
                    if int(line[2]) <= end:
                        score = line[4]

    shapes = []
    bed = pysam.TabixFile(bed_path)
    records = set()
    for record in bed.fetch(chrom, start, end):

        line = record.split("\t")

        coo = line[0:3]
        coo = "\t".join(coo)
        records.add(coo)

    for record in records:
        coo = record.split("\t")
        color = "Crimson"
        fill = "Salmon"
        shapes.append(
            dict(
                type="rect",
                x0=coo[1],
                y0=0,
                x1=coo[2],
                y1=1,
                line=dict(color=color, width=2),
                fillcolor=fill,
            )
        )

    ylim = [-1, 2]

    return ylim, shapes


def make_modbam_trace_gl(
    dicts, start, end, heterogeneity=None, marker_size=6, single_trace_height=12
):
    colors = px.colors.qualitative.T10
    freq_traces = []
    single_read_traces = []
    traces_height = []
    het_traces = []
    # single_trace_height = 12  # px

    for i, sample_dict in enumerate(dicts):
        freq = plot_frequencies_gl(sample_dict, start, end, color=colors[i])
        freq_traces.append(freq)

        if heterogeneity:
            het_trace = get_dict_heterogeneity(sample_dict, start, end, color=colors[i])
            het_traces.append(het_trace)

        sample_dict = queue_reads_plotly(sample_dict)
        traces = []
        for line, reads in sample_dict.items():
            for read in reads:
                traces.append(
                    go.Scattergl(
                        mode="lines+markers",
                        line=dict(color=colors[i], width=marker_size / 2),
                        x=list(read[1][2].keys()),
                        y=np.full(len(read[1][2].keys()), line),
                        connectgaps=True,
                        # marker = {'color': list(map(SetColor,list(read[1][2].values()))),
                        #             'size': 6,
                        #             'symbol': "square",
                        #             'line':{
                        # 'color':colors[i],
                        # 'width':0.2}
                        #             },
                        marker={
                            "color": list(map(SetColor, list(read[1][2].values()))),
                            "size": marker_size,
                            "symbol": "square",
                        },
                        name=read[0],
                        showlegend=False,
                    )
                )
        height = len(sample_dict)
        single_read_traces.append([traces, height * single_trace_height])

    return freq_traces, single_read_traces, het_traces


def get_tracks_gl(
    chrom,
    start,
    end,
    dicts=None,
    gtfs=None,
    beds=None,
    bigwigs=None,
    bedgraphs=None,
    heterogeneity=None,
    marker_size=6,
    single_trace_height=12,
):
    tracks = {}
    num_tracks = 0
    if gtfs:
        tracks["gtf"] = []
        for gtf in gtfs:
            genes = parse_gtf_exons(gtf, chrom, start, end)
            tracks["gtf"].append(genes)
            num_tracks += 1
    if beds:
        tracks["bed"] = []
        for bed in beds:
            elements = parse_bed_rectangle_gl(bed, chrom, start, end)
            tracks["bed"].append(elements)
            num_tracks += 1
    if bigwigs:
        tracks["bigwig"] = []
        for bigwig in bigwigs:
            bw = parse_bigwig_gl(bigwig, chrom, start, end)
            tracks["bigwig"].append(bw)
            num_tracks += 1
    if bedgraphs:
        tracks["bedgraph"] = []
        for bedgraph in bedgraphs:
            b = parse_bedgraph_gl(bedgraph, chrom, start, end)
            tracks["bedgraph"].append(b)
            num_tracks += 1

    if dicts:
        tracks["heterogeneity"] = []
        tracks["modbase_freq"] = []
        tracks["modbase"] = []
        if heterogeneity:
            freq_traces, single_read_traces, het_traces = make_modbam_trace_gl(
                dicts,
                start,
                end,
                heterogeneity,
                marker_size,
                single_trace_height=single_trace_height,
            )
            tracks["heterogeneity"] = het_traces
            tracks["modbase_freq"] = freq_traces
            tracks["modbase"] = single_read_traces
            num_tracks += len(dicts) + 2
        else:
            freq_traces, single_read_traces, het_traces = make_modbam_trace_gl(
                dicts,
                start,
                end,
                marker_size=marker_size,
                single_trace_height=single_trace_height,
            )
            tracks["modbase_freq"] = freq_traces
            tracks["modbase"] = single_read_traces
            num_tracks += len(dicts) + 1
    # print(tracks.keys(), num_tracks)
    return tracks, num_tracks


# def get_heights(tracks):
#     heights = []
#     for track_type, sub_tracks in tracks.items():
#         if (track_type == "heterogeneity") & (len(sub_tracks) != 0):
#             heights.append(sub_tracks[0][-1])
#             continue
#         if track_type == "modbase_freq":
#             heights.append(sub_tracks[0][-1])
#             continue
#         for track in sub_tracks:
#             heights.append(track[-1])
#     plot_height = sum(heights)
#     row_heights = [x / plot_height for x in heights]
#     # print(row_heights)
#     return plot_height, row_heights
