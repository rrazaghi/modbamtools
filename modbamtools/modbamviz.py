from modbamtools.utils import *
from modbamtools.tracks import *
from modbamtools.gene_models import *


class Plotter:
    """
    The Plotter class produces Plotly figures from Modbam
    """

    def __init__(self,dicts,samp_names,chrom, start, end, gtf) -> None:
        self.chrom = chrom
        self.start = start
        self.end = end
        self.dicts = dicts
        self.samp_names = samp_names
        self.num_samples = len(dicts)
        self.gtf = gtf
        self.row_h = calc_rows_h(dicts)
        self.tracks_titles = ["Genes","Methylation frequency plots"]
        self.fig = make_subplots(rows=self.num_samples + 2,cols=1,shared_xaxes=True,
        vertical_spacing=0.02,row_heights=[0.1 ,0.2] + list(self.row_h),
        subplot_titles=self.tracks_titles + samp_names )
        self.colors = px.colors.qualitative.T10

    def plot_genes(self):
        """
        Plot gene models
        """
        genes = parse_gtf(self.gtf, self.chrom, self.start, self.end)
        for name_trace, shape in zip(genes[1],genes[2]):
            self.fig.append_trace(name_trace, row=1, col=1)
            self.fig.add_shape(shape, row=1, col=1)
        self.fig.update_xaxes(visible=False, row=1, col=1)
        self.fig.update_yaxes(range=genes[0],visible=False, row=1, col=1)

    def plot(self):
        """
        Plot reads and sample frequencies
        """
        for i, sample_dict in enumerate(self.dicts):
            #### Plot Frequencies ####
            freq = plot_frequencies(sample_dict, self.start, self.end, color=self.colors[i])
            self.fig.add_traces(freq, rows=[2,2], cols=[1,1])
            # self.fig.update_yaxes(visible=False, row= 2, col=1)
            self.fig.update_xaxes(visible=False, row=i + 2, col=1)
    
            #### Plot Reads ####
            sample_dict = queue_reads_plotly(sample_dict)
            for line,reads in sample_dict.items():
                for read in reads:
                    self.fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=self.colors[i]),
                    x = list(read[1][2].keys()),
                    y = np.full(len(read[1][2].keys()), line),
                    connectgaps=True,
                    marker = {'color': list(map(SetColor,list(read[1][2].values()))),
            #                      'colorscale': colorscale,
                                'size': 6,
                                'symbol': "square"
                                },
                    name= read[0]  ,showlegend=False
                    ), row=i+3, col=1) 
                    
            self.fig.update_xaxes(visible=False, row=i + 3, col=1)
            self.fig.update_yaxes(visible=False, row=i + 3, col=1)

        self.fig.update_xaxes(visible=True, row=i + 3, col=1)
        self.fig.update_xaxes(range=[self.start, self.end],tickformat=',d',title_text="Coordinate")
