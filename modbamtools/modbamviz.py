from modbamtools.utils import *
from modbamtools.tracks import *
from modbamtools.tracks_webgl import *
from modbamtools.gene_models import *


class Plotter:
    """
    The Plotter class produces Plotly figures from Modbam
    """

    def __init__(
        self,
        dicts,
        samp_names,
        chrom,
        start,
        end,
        gtfs=None,
        beds=None,
        bigwigs=None,
        bedgraphs=None,
        track_titles=None,
        heterogeneity=None,
        font_size=18,
        fmt="html",
        marker_size=6,
        single_trace_height=12,
    ) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.dicts = dicts
        self.gtfs = gtfs
        self.beds = beds
        self.bigwigs = bigwigs
        self.bedgraphs = bedgraphs
        self.samp_names = samp_names
        self.track_titles = track_titles
        self.heterogeneity = heterogeneity
        self.font_size = font_size
        self.fmt = fmt
        self.marker_size = marker_size
        self.single_trace_height = single_trace_height

        if self.fmt == "html":
            self.tracks, self.num_tracks = get_tracks_gl(
                self.chrom,
                self.start,
                self.end,
                self.dicts,
                self.gtfs,
                self.beds,
                self.bigwigs,
                self.bedgraphs,
                self.heterogeneity,
                self.marker_size,
                self.single_trace_height,
            )
        else:
            self.tracks, self.num_tracks = get_tracks(
                self.chrom,
                self.start,
                self.end,
                self.dicts,
                self.gtfs,
                self.beds,
                self.bigwigs,
                self.bedgraphs,
                self.heterogeneity,
                self.marker_size,
                self.single_trace_height,
            )
        # print(self.num_tracks, self.track_titles)
        self.plot_height, self.row_heights = get_heights(self.tracks)
        if self.track_titles:
            if self.heterogeneity:
                self.titles = (
                    self.track_titles
                    + ["Methylation Heterogeneity"]
                    + ["Methylation Frequency"]
                    + samp_names
                )
            else:
                self.titles = self.track_titles + ["Methylation Frequency"] + samp_names
        if not self.track_titles:
            if self.heterogeneity:
                self.titles = (
                    [""] * (self.num_tracks - len(dicts) - 1)
                    + ["Methylation Heterogeneity"]
                    + ["Methylation Frequency"]
                    + self.samp_names
                )
            else:
                self.titles = (
                    [""] * (self.num_tracks - len(dicts) - 1)
                    + ["Methylation Frequency"]
                    + self.samp_names
                )
        assert self.num_tracks == len(self.titles)
        # self.tracks_titles = ["Genes","Enhancers","Methylation frequency plots"]
        self.fig = make_subplots(
            rows=self.num_tracks,
            cols=1,
            shared_xaxes=True,
            vertical_spacing=50 / self.plot_height,
            row_heights=self.row_heights,
            subplot_titles=self.titles,
            x_title="Coordinate" + " (" + self.chrom + ")",
        )
        self.fig.update_layout(font=dict(size=self.font_size))
        self.fig.update_layout(height=self.plot_height)

    def plot_tracks(self):
        """
        Plot tracks
        """
        i = 1
        if self.gtfs:
            for genes in self.tracks["gtf"]:
                for name_trace in genes[1]:
                    self.fig.append_trace(name_trace, row=i, col=1)
                for shape in genes[2]:
                    self.fig.add_shape(shape, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                self.fig.update_yaxes(range=genes[0], visible=False, row=i, col=1)
                i += 1
        if self.beds:
            for elements in self.tracks["bed"]:
                for name_trace in elements[1]:
                    self.fig.append_trace(name_trace, row=i, col=1)
                for shape in elements[2]:
                    self.fig.add_shape(shape, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                self.fig.update_yaxes(range=elements[0], visible=False, row=i, col=1)
                i += 1
        if self.bigwigs:
            for bw in self.tracks["bigwig"]:
                self.fig.append_trace(bw[0], row=i, col=1)
                self.fig.update_yaxes(visible=False, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                i += 1
        if self.bedgraphs:
            for b in self.tracks["bedgraph"]:
                self.fig.append_trace(b[0], row=i, col=1)
                self.fig.update_yaxes(visible=False, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                i += 1
        if self.dicts:
            if self.heterogeneity:
                het_row = i
                freq_row = i + 1
                i += 2
                for freq_traces, single_read_traces, het_traces in zip(
                    self.tracks["modbase_freq"],
                    self.tracks["modbase"],
                    self.tracks["heterogeneity"],
                ):
                    self.fig.add_traces(
                        freq_traces[0], rows=[freq_row, freq_row], cols=[1, 1]
                    )
                    self.fig.update_xaxes(visible=False, row=freq_row, col=1)

                    self.fig.add_traces(
                        het_traces[0], rows=[het_row, het_row], cols=[1, 1]
                    )
                    self.fig.update_xaxes(visible=False, row=het_row, col=1)

                    for trace in single_read_traces[0]:
                        self.fig.add_trace(trace, row=i, col=1)
                    self.fig.update_xaxes(visible=False, row=i, col=1)
                    self.fig.update_yaxes(visible=False, row=i, col=1)
                    i += 1
            else:
                freq_row = i
                i += 1
                for freq_traces, single_read_traces in zip(
                    self.tracks["modbase_freq"], self.tracks["modbase"]
                ):
                    self.fig.add_traces(
                        freq_traces[0], rows=[freq_row, freq_row], cols=[1, 1]
                    )
                    self.fig.update_xaxes(visible=False, row=freq_row, col=1)
                    # self.fig.update_yaxes(range=(-20, 120), row=freq_row, col=1)

                    for trace in single_read_traces[0]:
                        self.fig.add_trace(trace, row=i, col=1)
                    self.fig.update_xaxes(visible=False, row=i, col=1)
                    self.fig.update_yaxes(visible=False, row=i, col=1)
                    self.fig.update_yaxes(
                        range=(
                            -1 * (single_read_traces[1] / self.single_trace_height + 1),
                            1,
                        ),
                        visible=False,
                        row=i,
                        col=1,
                    )

                    i += 1

        self.fig.update_xaxes(visible=True, row=i - 1, col=1)
        self.fig.update_xaxes(
            range=[self.start, self.end],
            # tickformat=",d",
            # title_text="Coordinate" + " (" + self.chrom + ")",
            row=i - 1,
            col=1,
        )
        self.fig.update_annotations(font_size=self.font_size)
        self.fig.update_layout(
            height=self.plot_height,
            # width=8 * 96,
            # autosize=False,
            margin=dict(l=10, r=10, t=30, b=60),
        )
        # self.fig.update_yaxes(constrain="domain")
