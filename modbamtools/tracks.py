from modbamtools.utils import *

# def parse_bigwig(bigwig_path, chrom, start,end):
    
#     # bw = pyBigWig.open(bigwig_path)
#     x = list(range(start,end))
#     y = bw.values(chrom, start, end)
    
#     return go.Scattergl(
#         x=x,
#         y=y,
#         mode='lines',
#         marker=dict(
#             size=6,
#             color="goldenrod",
#         ),
#         name=''
#     )
    
    
def plot_frequencies(dict_per_read_mod, start, end, color):
    
    traces = []
    freq, freq_smooth = calc_freq(dict_per_read_mod, start, end)
    
    traces.append(go.Scattergl(
        x=freq["x"],
        y=freq["y"],
        mode='markers',
        marker=dict(
            size=3,
            color=color,
        ),
        name='',showlegend=False
    ))
    
    traces.append(go.Scattergl(
        x=freq_smooth["x"],
        y=freq_smooth["y"],
        mode='lines',
        marker=dict(
            size=3,
            color=color,
        ),
        name='',showlegend=False
    ))
    
    return traces

def plot_freq_diff(dict_per_read_mod_hp1,dict_per_read_mod_hp2, start, end, color="grey"):
    
    freq_hp1, freq_smooth_hp1 = calc_freq(dict_per_read_mod_hp1, start, end)
    freq_hp2, freq_smooth_hp2 = calc_freq(dict_per_read_mod_hp2, start, end)
    
    return go.Scattergl(
        x=freq_smooth_hp2["x"],
        y=[a_i - b_i for a_i, b_i in zip(freq_smooth_hp2["y"], freq_smooth_hp1["y"])],
        mode='lines',
        marker=dict(
            size=6,
            color=color,
        ),
        name=''
    )