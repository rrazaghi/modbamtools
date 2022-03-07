import pandas as pd
from scipy.signal import savgol_filter
import plotly.graph_objects as go
import numpy as np
import pysam
# import pyBigWig
import plotly.graph_objects as go
import plotly.io as pio
import collections
import pandas as pd
from modbampy import *
import plotly.express as px
from plotly.subplots import make_subplots
pio.templates.default = "simple_white"



def binerize_mod_call(call,min_prob=0.5, max_prob=0.5):
    prob = call/255
    if prob < min_prob:
        return 0
    elif prob > max_prob:
        return 1 
    else:
        return -1
    
    
def overlaps(a, b):
    """
    Return the amount of overlap, in bp
    between a and b.
    If >0, the number of bp of overlap
    If 0,  they are book-ended.
    If <0, the distance in bp between them
    """

    return min(a[1], b[1]) - max(a[0], b[0])

def calc_freq(dict_per_read_mod, start,end):

    df = pd.DataFrame.from_dict({k:v[2] for k,v in dict_per_read_mod.items()}, orient='index')
    df = df[df.columns[df.columns.isin(range(start,end))]]
    df = df.reindex(sorted(df.columns), axis=1)
    
    freq = {"x":[],"y":[]} ; freq_smooth = {"x":[],"y":[]}
    for pos in df.columns:
        count = df[pos].value_counts(dropna=True).to_dict()
        if 0 not in count.keys():
            count[0] = 0
        if 1 not in count.keys():
            count[1] = 0    
        perc_meth = count[1] * 100 /(count[0] + count[1])
        freq["x"].append(pos)
        freq["y"].append(perc_meth)
        
#     freq_smooth["y"] = lowess(freq["y"], freq["x"], frac=0.1,return_sorted=False)
    freq_smooth["y"] = savgol_filter(freq["y"], 51, 3)
    freq_smooth["x"] = freq["x"]
    return freq, freq_smooth

def calc_rows_h(dicts):
    num_reads = []
    for read_dict in dicts:
        num_reads.append(len(read_dict))

    row_h = list(0.85 * np.array(num_reads)/sum(num_reads))

    return row_h

def queue_reads_plotly(dict_per_read_mod):
    
    sorted_mod = dict(sorted(dict_per_read_mod.items(), key=lambda e: e[1][0], reverse=True))
    i = 0
    plot_dict = {}
    
    for k,v in sorted_mod.items():

        if len(plot_dict) == 0:
            
            plot_dict[i] = [(k,v)]
            i -= 1
            continue
        for line,reads in plot_dict.items():
            t =0
            for read in reads:
                if overlaps(v[1], read[1][1]) > 0:
                    t = 1
            if t == 0:
                plot_dict[line].append((k,v))
                break
                
        if t == 1:
            plot_dict[i] = [(k,v)]
            i -= 1
                
                

    return plot_dict

def SetColor(x):
    if(x == 0):
        return "blue"
#         return "#007da3"
    if(x == 1):
        return "red"
#         return "#DD4765"

def process_bam(bam, chrom, start, end,tag_name=None, tag_value =None,min_prob=0.5, max_prob=0.5):
#     tag_name = "HP"
#     tag_value = 2
    dict_per_read_mod = {}
    with ModBam(bam) as b:
            for read in b.reads(chrom, start, end, tag_name=tag_name,tag_value=tag_value):
                mapped_modbase = {}
                read_start = max([read.reference_start,start])
                read_end = min([read.reference_end,end])
                read_len = read_end - read_start
                for pos_mod in read.mod_sites:
                    qname, rpos, qpos, strand, mod_strand, cbase, mbase, score = pos_mod
                    if strand == "-":
                        call = binerize_mod_call(score,min_prob,max_prob)
                        if call != -1:
                            mapped_modbase[rpos - 1] = call
                    else:
                        call = binerize_mod_call(score,min_prob,max_prob)
                        if call != -1:
                            mapped_modbase[rpos] = call
                dict_per_read_mod[read.query_name] = [read_len,(read_start,read_end),collections.OrderedDict(sorted(mapped_modbase.items()))]
                
    return dict_per_read_mod

def get_reads(bams, chrom, start, end, hap=None, strand=None, samp_names = None, min_prob=0.5, max_prob=0.5):

    dicts = []
    if hap:
        if not samp_names:
            titles = []
            samp_list = ["Sample " + str(i) for i in range(1, len(bams) + 1)]
            for samp in samp_list:
                titles.append(samp + " haplotype 1")
                titles.append(samp + " haplotype 2")
        if samp_names:
            titles = []
            samp_list = samp_names
            for samp in samp_list:
                titles.append(samp + " haplotype 1")
                titles.append(samp + " haplotype 2")

        for bam in bams:
            reads_hp1 = process_bam(bam,chrom,start,end,tag_name="HP", tag_value=1,min_prob=min_prob,max_prob=max_prob)
            reads_hp2 = process_bam(bam,chrom,start,end,tag_name="HP", tag_value=2,min_prob=min_prob,max_prob=max_prob)
            dicts.append(reads_hp1)
            dicts.append(reads_hp2)



    else:
        
        if not samp_names:
            titles = ["Sample " + str(i) for i in range(1, len(bams) + 1)]
        if samp_names:
            titles = samp_names
        for bam in bams:
            reads = process_bam(bam,chrom,start,end,min_prob=min_prob,max_prob=max_prob)
            dicts.append(reads)

    return dicts, titles
