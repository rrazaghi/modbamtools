from modbamtools.utils import *
from operator import itemgetter
from itertools import groupby
from scipy.signal import savgol_filter


def binerize_mod_het_call(call, min_prob=0.5, max_prob=0.5):
    prob = call / 255
    if prob < min_prob:
        return -1
    elif prob > max_prob:
        return 1
    else:
        return 0


def get_read_heterogeneity(read, min_calls, start=None, end=None):

    mod_calls = []
    for pos_mod in read.mod_sites:
        qname, rpos, qpos, strand, mod_strand, cbase, mbase, score = pos_mod
        if start & end:
            if (rpos >= start) & (rpos <= end):
                call = binerize_mod_het_call(score)
                mod_calls.append(call)
        else:
            call = binerize_mod_het_call(score)
            mod_calls.append(call)

    if len(mod_calls) < min_calls:
        het = -1
    else:
        zero_crossings = np.where(np.diff(np.sign(mod_calls)))[0]
        ranges = []
        for k, g in groupby(enumerate(zero_crossings), lambda x: x[0] - x[1]):
            group = list(map(itemgetter(1), g))
            if len(group) > 1:
                ranges.append(range(group[0], group[-1]))
            else:
                ranges.append(group[0])
        het = len(ranges)

    return het


def get_region_heterogeneity(
    bam, chrom, start, end, min_calls, min_cov, tag_name=None, tag_value=None
):

    read_stats = []
    with ModBam(bam) as b:
        for read in b.reads(chrom, start, end, tag_name=tag_name, tag_value=tag_value):
            read_start = max([read.reference_start, start])
            read_end = min([read.reference_end, end])
            percent_cov = (read_end - read_start) * 100 / (end - start)
            if percent_cov <= min_cov:
                continue
            het = get_read_heterogeneity(read, min_calls, start, end)
            if het == -1:
                continue
            read_stats.append(het)
    if len(read_stats) != 0:
        read_stats = np.array(read_stats)
        mean = np.mean(read_stats)
        std = np.std(read_stats)
        cov = len(read_stats)
        return mean, std, cov
    else:
        return -1, -1, -1


def get_region_hap_heterogeneity(bam, min_calls, min_cov, hp, bed_line):

    chrom = bed_line[0]
    start = int(bed_line[1])
    end = int(bed_line[2])

    if hp:

        mean, std, cov = get_region_heterogeneity(
            bam, chrom, start, end, min_calls, min_cov
        )
        hp1_mean, hp1_std, hp1_cov = get_region_heterogeneity(
            bam, chrom, start, end, min_calls, min_cov, tag_name="HP", tag_value=1
        )
        hp2_mean, hp2_std, hp2_cov = get_region_heterogeneity(
            bam, chrom, start, end, min_calls, min_cov, tag_name="HP", tag_value=2
        )

        bed_line.extend(
            map(
                str,
                [
                    mean,
                    std,
                    cov,
                    hp1_mean,
                    hp1_std,
                    hp1_cov,
                    hp2_mean,
                    hp2_std,
                    hp2_cov,
                ],
            )
        )

    else:

        mean, std, cov = get_region_heterogeneity(
            bam, chrom, start, end, min_calls, min_cov
        )
        bed_line.extend(
            map(
                str,
                [mean, std, cov],
            )
        )

    return bed_line


def get_dict_heterogeneity(samp_dict, start, end, color):

    window = 100
    binsize = 1000
    regions = []
    het_dict = {"x": [], "y": [], "y_smooth": []}
    for i in range(start, end, window):
        if i + binsize > end:
            regions.append((i, end))
            break
        else:
            regions.append((i, i + binsize))
            if i + binsize == end:
                break

    for region in regions:
        hets = []
        for rname, arr in samp_dict.items():
            signs = []
            for pos, value in arr[2].items():
                if (pos >= region[0]) & (pos <= region[1]):
                    if value == 0:
                        signs.append(-1)
                    if value == 1:
                        signs.append(1)

            if len(signs) <= 5:
                continue
            else:
                zero_crossings = np.where(np.diff(np.sign(signs)))[0]
                ranges = []
                for k, g in groupby(enumerate(zero_crossings), lambda x: x[0] - x[1]):
                    group = list(map(itemgetter(1), g))
                    if len(group) > 1:
                        ranges.append(range(group[0], group[-1]))
                    else:
                        ranges.append(group[0])

                hets.append(len(ranges))

        if len(hets) <= 5:
            continue
        else:
            het_dict["x"].append((region[0] + region[1]) / 2)
            het_dict["y"].append(np.mean(hets))

    length = len(het_dict["y"])
    if length <= 5:
        smooth_window = 1
        poly = 0
    elif (length > 5) & (length <= 20):
        smooth_window = 5
        poly = 3
    elif (length > 20) & (length <= 50):
        smooth_window = 21
        poly = 3
    else:
        smooth_window = 51
        poly = 3
    het_dict["y_smooth"] = savgol_filter(het_dict["y"], smooth_window, poly)

    height = 200  # px

    traces = []

    traces.append(
        go.Scatter(
            x=het_dict["x"],
            y=het_dict["y"],
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
        go.Scatter(
            x=het_dict["x"],
            y=het_dict["y_smooth"],
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


def get_plot_heterogeneity(bam, chrom, start, end, hp, min_calls=5, min_cov=80):

    window = 20
    binsize = 100
    regions = []
    for i in range(start, end, window):
        if i + binsize > end:
            regions.append((i, end))
            break
        else:
            regions.append((i, i + binsize))
            if i + binsize == end:
                break
    length = len(regions)
    if length <= 5:
        smooth_window = 1
        poly = 0
    elif (length > 5) & (length <= 20):
        smooth_window = 5
        poly = 3
    elif (length > 20) & (length <= 50):
        smooth_window = 21
        poly = 3
    else:
        smooth_window = 51
        poly = 3

    het = {"all": [], "hp1": [], "hp2": [], "x": [], "x1": [], "x2": []}
    for region in regions:

        if hp:
            hp1_mean, hp1_std, hp1_cov = get_region_heterogeneity(
                bam,
                chrom,
                region[0],
                region[1],
                min_calls,
                min_cov,
                tag_name="HP",
                tag_value=1,
            )
            hp2_mean, hp2_std, hp2_cov = get_region_heterogeneity(
                bam,
                chrom,
                region[0],
                region[1],
                min_calls,
                min_cov,
                tag_name="HP",
                tag_value=2,
            )
            if hp1_mean != -1:
                het["x1"].append((region[0] + region[1]) / 2)
                het["hp1"].append(hp1_mean)
            if hp2_mean != -1:
                het["x2"].append((region[0] + region[1]) / 2)
                het["hp2"].append(hp2_mean)
            het["hp1_smooth"] = savgol_filter(het["hp1"], smooth_window, poly)
            het["hp2_smooth"] = savgol_filter(het["hp2"], smooth_window, poly)
        else:
            mean, std, cov = get_region_heterogeneity(
                bam, chrom, region[0], region[1], min_calls, min_cov
            )
            if mean != -1:
                het["x"].append((region[0] + region[1]) / 2)
                het["all"].append(mean)
            het["all_smooth"] = savgol_filter(het["all"], smooth_window, poly)

        track_height = 100  # px

        return het, track_height
