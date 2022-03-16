from modbamtools.utils import *
from modbamtools.tracks import *
from modbamtools.gene_models import *

class Plotter:
    """
    The Plotter class produces Plotly figures from Modbam
    """

    def __init__(self,dicts,samp_names,chrom, start,
     end, gtfs=None, beds=None, bigwigs=None, bedgraphs=None, track_titles=None) -> None:
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
        self.tracks , self.num_tracks = get_tracks(self.chrom, self.start, self.end, 
        self.dicts,self.gtfs,self.beds,
        self.bigwigs,self.bedgraphs)
        self.plot_height , self.row_heights = get_heights(self.tracks)
        if self.track_titles:
            self.titles = self.track_titles + [""] + samp_names
        if not self.track_titles:
            self.titles = [""] * (self.num_tracks - len(dicts)) + self.samp_names

        # self.tracks_titles = ["Genes","Enhancers","Methylation frequency plots"]
        self.fig = make_subplots(rows=self.num_tracks, cols=1, shared_xaxes=True,
        vertical_spacing=0.04,row_heights=self.row_heights,
        subplot_titles=self.titles )
    def plot_tracks(self):
        """
        Plot tracks
        """
        i = 1
        if self.gtfs:
            for genes in self.tracks['gtf']:
                for name_trace in genes[1]:
                    self.fig.append_trace(name_trace, row=i, col=1)
                for shape in genes[2]:
                    self.fig.add_shape(shape, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                self.fig.update_yaxes(range=genes[0],visible=False, row=i, col=1)
                i += 1
        if self.beds:
            for elements in self.tracks['bed']:
                for shape in elements[1]:
                    self.fig.add_shape(shape, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                self.fig.update_yaxes(visible=False, row=i, col=1)
                self.fig.update_yaxes(range=elements[0],visible=False, row=i, col=1)
                i += 1
        if self.bigwigs:
            for bw in self.tracks['bigwig']:
                self.fig.append_trace(bw[0], row=i, col=1)
                self.fig.update_yaxes(visible=False, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                i += 1
        if self.bedgraphs:
            for b in self.tracks['bedgraph']:
                self.fig.append_trace(b[0], row=i, col=1)
                self.fig.update_yaxes(visible=False, row=i, col=1)
                self.fig.update_xaxes(visible=False, row=i, col=1)
                i += 1
        if self.dicts:
            freq_row = i
            i += 1
            for freq_traces, single_read_traces in zip(self.tracks['modbase_freq'], self.tracks['modbase']):
                self.fig.add_traces(freq_traces[0], rows=[freq_row,freq_row], cols=[1,1])
                self.fig.update_xaxes(visible=False, row=freq_row, col=1)

                for trace in single_read_traces[0]:
                    self.fig.add_trace(trace,row=i, col=1 )
                self.fig.update_xaxes(visible=False, row=i , col=1)
                self.fig.update_yaxes(visible=False, row=i , col=1)
                i += 1

        self.fig.update_xaxes(visible=True, row=i-1, col=1)
        self.fig.update_xaxes(range=[self.start, self.end],tickformat=',d',title_text="Coordinate")   
        self.fig.update_layout(height=self.plot_height)
    # def plot(self):
    #     """
    #     Plot reads and sample frequencies
    #     """
    #     for i, sample_dict in enumerate(self.dicts):
    #         #### Plot Frequencies ####
    #         freq = plot_frequencies(sample_dict, self.start, self.end, color=self.colors[i])
    #         self.fig.add_traces(freq, rows=[3,3], cols=[1,1])
    #         # self.fig.update_yaxes(visible=False, row= 2, col=1)
    #         self.fig.update_xaxes(visible=False, row=i + 3, col=1)
    
    #         #### Plot Reads ####
    #         sample_dict = queue_reads_plotly(sample_dict)
    #         for line,reads in sample_dict.items():
    #             for read in reads:
    #                 self.fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=self.colors[i]),
    #                 x = list(read[1][2].keys()),
    #                 y = np.full(len(read[1][2].keys()), line),
    #                 connectgaps=True,
    #                 marker = {'color': list(map(SetColor,list(read[1][2].values()))),
    #         #                      'colorscale': colorscale,
    #                             'size': 6,
    #                             'symbol': "square"
    #                             },
    #                 name= read[0]  ,showlegend=False
    #                 ), row=i+4, col=1) 
                    
    #         self.fig.update_xaxes(visible=False, row=i + 4, col=1)
    #         self.fig.update_yaxes(visible=False, row=i + 4, col=1)

    #     self.fig.update_xaxes(visible=True, row=i + 4, col=1)
    #     self.fig.update_xaxes(range=[self.start, self.end],tickformat=',d',title_text="Coordinate")
