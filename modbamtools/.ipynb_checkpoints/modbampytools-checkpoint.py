from modbampy import *
from utils import *
from tracks import *
from gene_models import *
from palettes import *
import pandas as pd
import hdbscan


def process_bam(bam, chrom, start, end,tag_name=None, tag_value =None):
#     tag_name = "HP"
#     tag_value = 2
    dict_per_read_mod = {}
    with ModBam(bam, chrom, start, end, tag_name=tag_name,tag_value=tag_value) as b:
            for read in b.reads():
                mapped_modbase = {}
                read_start = max([read.reference_start,start])
                read_end = min([read.reference_end,end])
                read_len = read_end - read_start
                for pos_mod in read.mod_sites:
                    qname, rpos, qpos, strand, mod_strand, cbase, mbase, score = pos_mod
                    if strand == "-":
                        mapped_modbase[rpos - 1] = binerize_mod_call(score)
                    else:
                        mapped_modbase[rpos] = binerize_mod_call(score)
                dict_per_read_mod[read.query_name] = [read_len,(read_start,read_end),collections.OrderedDict(sorted(mapped_modbase.items()))]
                
    return dict_per_read_mod

def SetColor(x):
    if(x == 0):
        return "blue"
#         return "#007da3"
    if(x == 1):
        return "red"
#         return "#DD4765"
def plot_tracks(bam, gtf_file, bed_file, chrom, start, end,sample=""):
    dict_per_read_mod_hp1 = process_bam(bam ,chrom, start, end, tag_name="HP", tag_value=1)
    dict_per_read_mod_hp2 = process_bam(bam ,chrom, start, end, tag_name="HP", tag_value=2)

#     ac27 = parse_bigwig(h3k27ac_file,chrom,start,end)
#     me1 = parse_bigwig(h3k4me1_file,chrom,start,end)
    genes = parse_gtf(gtf_file, chrom, start, end)
    enhancers = parse_bed(bed_file, chrom, start, end)

    dict_per_read_mod_hp1_p = queue_reads_plotly(dict_per_read_mod_hp1)
    dict_per_read_mod_hp2_p = queue_reads_plotly(dict_per_read_mod_hp2)

#     hp1_c = "#1FD699"
#     hp2_c = "#CB902A"
    hp1_c = "#1FD699"
    hp2_c = "#CB902A"

#     sample = "Blood"
    sub_titles = ["Genes","Enhancers","Methylation frequency plots",sample +" (Haplotype 1)",sample +" (Haplotype 2)"]
#     sub_titles = ["Genes","Enhancers","Methylation frequency plots",sample +" (Haplotype 1)",sample +" (Haplotype 2)"]
    fig = make_subplots(rows=5, cols=1,shared_xaxes=True,vertical_spacing=0.025,row_heights=[0.05,0.05 ,0.2,
                                                                                                0.7*len(dict_per_read_mod_hp1_p)/(len(dict_per_read_mod_hp1_p)+len(dict_per_read_mod_hp2_p)),
                                                                                                0.7*len(dict_per_read_mod_hp2_p)/(len(dict_per_read_mod_hp1_p)+len(dict_per_read_mod_hp2_p))],
                            subplot_titles=sub_titles)
    # fig = make_subplots(rows=6, cols=1,shared_xaxes=True,vertical_spacing=0.025,
    #                         subplot_titles=["",""])
    ###############################
    #### Plot gene models ####
    ###############################
    for name_trace, shape in zip(genes[1],genes[2]):
        fig.append_trace(name_trace, row=1, col=1)
        fig.add_shape(shape, row=1, col=1)
    fig.update_xaxes(visible=False, row=1, col=1)
    fig.update_yaxes(range=genes[0],visible=False, row=1, col=1)

    ###############################
    #### Plot E-P ####
    ###############################
    for shape in enhancers[1]:
        fig.add_shape(shape, row=2, col=1)
    fig.update_xaxes(visible=False, row=2, col=1)
    fig.update_yaxes(visible=False, row=2, col=1)

    ###############################
    #### Plot Tracks ####
    ###############################
#     fig.add_trace(me1, row=3, col=1)
#     fig.update_xaxes(visible=False, row=3, col=1)
#     fig.update_yaxes(visible=False, row=3, col=1)

#     fig.add_trace(ac27, row=4, col=1)
#     fig.update_xaxes(visible=False, row=4, col=1)
#     fig.update_yaxes(visible=False, row=4, col=1)

    ###############################
    #### Plot Frequencies ####
    ###############################
#     freq_diff = plot_freq_diff(dict_per_read_mod_hp1,dict_per_read_mod_hp2, start, end)
#     fig.add_trace(freq_diff, row=5, col=1)
#     fig.update_xaxes(visible=False, row=5, col=1)
#     fig.update_yaxes(range=[-100,100], row=5, col=1)
#     fig.add_shape(type="line",
#         x0=start, y0=0, x1=end, y1=0,
#         line=dict(
#             color="MediumPurple",
#             width=4,
#             dash="dot",
#         ),row=5, col=1
#     )

    hp1_freq = plot_frequencies(dict_per_read_mod_hp1, start, end, color=hp1_c)
    hp2_freq = plot_frequencies(dict_per_read_mod_hp2, start, end, color=hp2_c)
    fig.add_traces(hp1_freq, rows=[3,3], cols=[1,1])
    fig.add_traces(hp2_freq, rows=[3,3], cols=[1,1])
    fig.update_xaxes(visible=False, row=3, col=1)

    ###############################
    #### Plot Reads ####
    ###############################

    for line,reads in dict_per_read_mod_hp1_p.items():
        for read in reads:
            fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=hp1_c),
            x = list(read[1][2].keys()),
            y = np.full(len(read[1][2].keys()), line),
            connectgaps=True,
            marker = {'color': list(map(SetColor,list(read[1][2].values()))),
    #                      'colorscale': colorscale,
                        'size': 6,
                        'symbol': "square"
                        },
            name= read[0]  ,showlegend=False
            ), row=4, col=1)

    for line,reads in dict_per_read_mod_hp2_p.items():
        for read in reads:
            fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=hp2_c),
            x = list(read[1][2].keys()),
            y = np.full(len(read[1][2].keys()), line),
            connectgaps=True,
            marker = {'color': list(map(SetColor,list(read[1][2].values()))),
    #                       'colorscale': colorscale,
                        'size': 6,
                        'symbol': "square"
                        },
            name= read[0]  ,showlegend=False
            ), row=5, col=1)
    fig.update_xaxes(visible=False, row=4, col=1)
    fig.update_yaxes(visible=False, row=4, col=1)
    # fig.update_xaxes(visible=False, row=8, col=1)
    fig.update_yaxes(visible=False, row=5, col=1)

    fig.update_xaxes(range=[start, end],tickformat=',d',title_text="Coordinate")
    fig.update_layout(height=1000, width=1500)
    
    return fig

def plot_tracks_multi(bams, gtf_file, chrom, start, end, samples=[""]):
    
    dicts = []
    num_reads = []
    for bam in bams:
        dict_per_read_mod = process_bam(bam ,chrom, start, end)
        num_reads.append(len(dict_per_read_mod))
        dicts.append(dict_per_read_mod)
    
    titles = samples
    tracks_titles = ["Genes","Methylation frequency plots"]
    sub_titles = tracks_titles + titles
    row_h = list(0.85 * np.array(num_reads)/sum(num_reads))
    fig = make_subplots(rows=len(bams) + 2, cols=1,shared_xaxes=True,vertical_spacing=0.01,row_heights=[0.05 ,0.1] + row_h,subplot_titles=sub_titles)    
    colors = met_brew(name="Monet",n=len(bams), brew_type="discrete")
    ###############################
    #### Plot gene models ####
    ###############################
    genes = parse_gtf(gtf_file, chrom, start, end)
    for name_trace, shape in zip(genes[1],genes[2]):
        fig.append_trace(name_trace, row=1, col=1)
        fig.add_shape(shape, row=1, col=1)
    fig.update_xaxes(visible=False, row=1, col=1)
    fig.update_yaxes(range=genes[0],visible=False, row=1, col=1)
    
    
    ###############################
    #### Plot clusters ####
    ###############################   
    
    for i, sample_dict in enumerate(dicts):
        
        #### Plot Frequencies ####
        freq = plot_frequencies(sample_dict, start, end, color=colors[i])
        fig.add_traces(freq, rows=2, cols=1)
#         fig.add_traces(freq, rows=[2,2], cols=[1,1])
        fig.update_xaxes(visible=False, row=i + 2, col=1)
#         print(cluster_dict)
        #### Plot Reads ####
        sample_dict = queue_reads_plotly(sample_dict)
        for line,reads in sample_dict.items():
            for read in reads:
                fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=colors[i]),
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
                
        fig.update_xaxes(visible=False, row=i + 3, col=1)
        fig.update_yaxes(visible=False, row=i + 3, col=1)

    fig.update_xaxes(visible=True, row=i + 3, col=1)
    fig.update_xaxes(range=[start, end],tickformat=',d',title_text="Coordinate")
    fig.update_layout(height=2000, width=1500)
    
    return fig

def cluster_plot(bam, gtf_file, chrom, start, end, sample=""):
    

    dict_per_read_mod = process_bam(bam ,chrom, start, end)
    reads_dict = {}
    for rname, arr in dict_per_read_mod.items():
        reads_dict[rname] = arr[2]
        
    df = pd.DataFrame.from_dict(reads_dict, orient='index')
    df = df[df.columns[df.columns.isin(range(start,end))]]
    df = df.reindex(sorted(df.columns), axis=1)
    df.fillna(0, inplace=True)
    clusterer = hdbscan.HDBSCAN()
    clusterer.fit(df.to_numpy())
    clusters = {}
    for read, cluster in zip(list(df.index),clusterer.labels_):
        clusters[read] = cluster
    cl_labels = list(set(clusterer.labels_))
    num_cl = len(cl_labels)
    print(cl_labels)
    
    hp1_c = "#1FD699"
    hp2_c = "#CB902A"
    cl_titles = [sample + " (Cluster " + str(label) + ")" for label in cl_labels]
    tracks_titles = ["Genes","Methylation frequency plots"]
    sub_titles = tracks_titles + cl_titles
    row_h = list(0.7 * np.array([cl_labels.count(label) for label in sorted(list(set(cl_labels)))])/sum([cl_labels.count(label) for label in sorted(list(set(cl_labels)))]))
    fig = make_subplots(rows=num_cl + 2, cols=1,shared_xaxes=True,vertical_spacing=0.025,row_heights=[0.1 ,0.2] + row_h,subplot_titles=sub_titles)    
    colors = met_brew(name="Monet",n=num_cl, brew_type="discrete")
    ###############################
    #### Plot gene models ####
    ###############################
    genes = parse_gtf(gtf_file, chrom, start, end)
    for name_trace, shape in zip(genes[1],genes[2]):
        fig.append_trace(name_trace, row=1, col=1)
        fig.add_shape(shape, row=1, col=1)
    fig.update_xaxes(visible=False, row=1, col=1)
    fig.update_yaxes(range=genes[0],visible=False, row=1, col=1)
    
    
    ###############################
    #### Plot clusters ####
    ###############################   
    
    for i, label in enumerate(cl_labels):
        
        #### Plot Frequencies ####
        cluster_dict = {}
        for rname, dict_mod in dict_per_read_mod.items():
            if clusters[rname] == label:
                cluster_dict[rname] = dict_mod
        freq = plot_frequencies(cluster_dict, start, end, color=colors[i])
        fig.add_traces(freq, rows=[2,2], cols=[1,1])
        fig.update_xaxes(visible=False, row=i + 2, col=1)
#         print(cluster_dict)
        #### Plot Reads ####
        cluster_dict = queue_reads_plotly(cluster_dict)
        for line,reads in cluster_dict.items():
            for read in reads:
                fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=colors[i]),
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
                
        fig.update_xaxes(visible=False, row=i + 3, col=1)
        fig.update_yaxes(visible=False, row=i + 3, col=1)

    fig.update_xaxes(visible=True, row=i + 3, col=1)
    fig.update_xaxes(range=[start, end],tickformat=',d',title_text="Coordinate")
    fig.update_layout(height=1000, width=1000)
    
    return fig

def cluster_plot_multi(bams, gtf_file, chrom, start, end, sample="", height=2000, width=1500):
    dict_per_read_mod = {}
    for bam in bams:
        dict_mod = process_bam(bam ,chrom, start, end)
        dict_per_read_mod = {**dict_per_read_mod, **dict_mod}
    reads_dict = {}
    for rname, arr in dict_per_read_mod.items():
        if arr[0]/(end-start) > 0.9:
            reads_dict[rname] = arr[2]
        
    df = pd.DataFrame.from_dict(reads_dict, orient='index')
    df = df[df.columns[df.columns.isin(range(start,end))]]
    df = df.reindex(sorted(df.columns), axis=1)
    df.fillna(0, inplace=True)
    clusterer = hdbscan.HDBSCAN(metric="hamming",min_cluster_size=10,min_samples=5)
    clusterer.fit(df.to_numpy())
    clusters = {}
    for read, cluster in zip(list(df.index),clusterer.labels_):
        clusters[read] = cluster
    cl_labels = sorted(list(clusterer.labels_))
    uniq_labels = sorted(list(set(clusterer.labels_)))
    num_cl = len([label for label in uniq_labels if label != -1])
    print(uniq_labels)
    
    hp1_c = "#1FD699"
    hp2_c = "#CB902A"
    cl_titles = [sample + " (Cluster " + str(label +1) + ")" for label in uniq_labels if label != -1]
    tracks_titles = ["Genes","Methylation frequency plots"]
    sub_titles = tracks_titles + cl_titles
    row_h = list(0.7 * np.array([cl_labels.count(label) for label in uniq_labels if label != -1])/sum([cl_labels.count(label) for label in uniq_labels if label != -1]))
    print(row_h)
    fig = make_subplots(rows=num_cl + 2, cols=1,shared_xaxes=True,vertical_spacing=0.02,row_heights=[0.1 ,0.2] + row_h,subplot_titles=sub_titles)    
    colors = met_brew(name="Monet",n=num_cl, brew_type="continueous")
    ###############################
    #### Plot gene models ####
    ###############################
    genes = parse_gtf(gtf_file, chrom, start, end)
    for name_trace, shape in zip(genes[1],genes[2]):
        fig.append_trace(name_trace, row=1, col=1)
        fig.add_shape(shape, row=1, col=1)
    fig.update_xaxes(visible=False, row=1, col=1)
    fig.update_yaxes(range=genes[0],visible=False, row=1, col=1)
    
    
    ###############################
    #### Plot clusters ####
    ###############################   
    
    for i, label in enumerate([label for label in uniq_labels if label != -1]):
#         if label == -1:
#             continue
        
        #### Plot Frequencies ####
        cluster_dict = {}
        for rname, dict_mod in dict_per_read_mod.items():
            if rname in clusters.keys():
                if clusters[rname] == label:
                    cluster_dict[rname] = dict_mod
        freq = plot_frequencies(cluster_dict, start, end, color=colors[i])
        fig.add_traces(freq, rows=2, cols=1)
#         fig.add_traces(freq, rows=[2,2], cols=[1,1])
        fig.update_xaxes(visible=False, row=i + 2, col=1)
#         print(cluster_dict)
        #### Plot Reads ####
        cluster_dict = queue_reads_plotly(cluster_dict)
        for line,reads in cluster_dict.items():
            for read in reads:
                fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=colors[i]),
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
                
        fig.update_xaxes(visible=False, row=i + 3, col=1)
        fig.update_yaxes(visible=False, row=i + 3, col=1)

    fig.update_xaxes(visible=True, row=i + 3, col=1)
    fig.update_xaxes(range=[start, end],tickformat=',d',title_text="Coordinate")
    fig.update_layout(height=height, width=width)
    
    return fig

def cluster_stats(bams, chrom, start, end, **kwargs):
    dict_per_read_mod = {}
    for bam in bams:
        dict_mod = process_bam(bam ,chrom, start, end)
        dict_per_read_mod = {**dict_per_read_mod, **dict_mod}
    reads_dict = {}
    for rname, arr in dict_per_read_mod.items():
        if arr[0]/(end-start) > 0.9:
            reads_dict[rname] = arr[2]
        
    df = pd.DataFrame.from_dict(reads_dict, orient='index')
    df = df[df.columns[df.columns.isin(range(start,end))]]
    df = df.reindex(sorted(df.columns), axis=1)
    df.fillna(0, inplace=True)
    clusterer = hdbscan.HDBSCAN(**kwargs)
    clusterer.fit(df.to_numpy())
    clusters = {}
    cl_labels = sorted(list(clusterer.labels_))
    uniq_labels = sorted(list(set(clusterer.labels_)))
    for label in uniq_labels:
        label_count = cl_labels.count(label)
        clusters[label] = label_count
    return clusters



def plot_tracks_nohap(bam, gtf_file, bed_file, chrom, start, end,sample=""):
    dict_per_read_mod = process_bam(bam ,chrom, start, end)

#     ac27 = parse_bigwig(h3k27ac_file,chrom,start,end)
#     me1 = parse_bigwig(h3k4me1_file,chrom,start,end)
    genes = parse_gtf(gtf_file, chrom, start, end)
    enhancers = parse_bed(bed_file, chrom, start, end)

    dict_per_read_mod_p = queue_reads_plotly(dict_per_read_mod)

#     hp1_c = "#1FD699"
#     hp2_c = "#CB902A"
    hp1_c = "#1FD699"
    hp2_c = "#CB902A"

#     sample = "Blood"
    sub_titles = ["Genes","Methylation frequency plots",sample]
#     sub_titles = ["Genes","Enhancers","Methylation frequency plots",sample +" (Haplotype 1)",sample +" (Haplotype 2)"]
    fig = make_subplots(rows=3, cols=1,shared_xaxes=True,vertical_spacing=0.025,row_heights=[0.1 ,0.2,
                                                                                                0.7],
                            subplot_titles=sub_titles)
    # fig = make_subplots(rows=6, cols=1,shared_xaxes=True,vertical_spacing=0.025,
    #                         subplot_titles=["",""])
    ###############################
    #### Plot gene models ####
    ###############################
    for name_trace, shape in zip(genes[1],genes[2]):
        fig.append_trace(name_trace, row=1, col=1)
        fig.add_shape(shape, row=1, col=1)
    fig.update_xaxes(visible=False, row=1, col=1)
    fig.update_yaxes(range=genes[0],visible=False, row=1, col=1)

    ###############################
    #### Plot E-P ####
    ###############################
#     for shape in enhancers[1]:
#         fig.add_shape(shape, row=2, col=1)
#     fig.update_xaxes(visible=False, row=2, col=1)
#     fig.update_yaxes(visible=False, row=2, col=1)

    ###############################
    #### Plot Tracks ####
    ###############################
#     fig.add_trace(me1, row=3, col=1)
#     fig.update_xaxes(visible=False, row=3, col=1)
#     fig.update_yaxes(visible=False, row=3, col=1)

#     fig.add_trace(ac27, row=4, col=1)
#     fig.update_xaxes(visible=False, row=4, col=1)
#     fig.update_yaxes(visible=False, row=4, col=1)

    ###############################
    #### Plot Frequencies ####
    ###############################
#     freq_diff = plot_freq_diff(dict_per_read_mod_hp1,dict_per_read_mod_hp2, start, end)
#     fig.add_trace(freq_diff, row=5, col=1)
#     fig.update_xaxes(visible=False, row=5, col=1)
#     fig.update_yaxes(range=[-100,100], row=5, col=1)
#     fig.add_shape(type="line",
#         x0=start, y0=0, x1=end, y1=0,
#         line=dict(
#             color="MediumPurple",
#             width=4,
#             dash="dot",
#         ),row=5, col=1
#     )

    hp1_freq = plot_frequencies(dict_per_read_mod, start, end, color=hp1_c)
    fig.add_traces(hp1_freq, rows=[2,2], cols=[1,1])
    fig.update_xaxes(visible=False, row=2, col=1)

    ###############################
    #### Plot Reads ####
    ###############################

    for line,reads in dict_per_read_mod_p.items():
        for read in reads:
            fig.add_trace(go.Scattergl(mode='lines+markers', line=dict(color=hp1_c),
            x = list(read[1][2].keys()),
            y = np.full(len(read[1][2].keys()), line),
            connectgaps=True,
            marker = {'color': list(map(SetColor,list(read[1][2].values()))),
    #                      'colorscale': colorscale,
                        'size': 6,
                        'symbol': "square"
                        },
            name= read[0]  ,showlegend=False
            ), row=3, col=1)

  
    fig.update_xaxes(visible=False, row=2, col=1)
    fig.update_yaxes(visible=False, row=2, col=1)
    # fig.update_xaxes(visible=False, row=8, col=1)
    fig.update_yaxes(visible=False, row=3, col=1)

    fig.update_xaxes(range=[start, end],tickformat=',d',title_text="Coordinate")
    fig.update_layout(height=1000, width=1500)
    
    return fig