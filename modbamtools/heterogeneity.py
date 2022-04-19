from modbamtools.utils import *
from operator import itemgetter
from itertools import groupby


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
        bed_line.extend(map(str, [mean, std, cov],))

    return bed_line


def get_plot_heterogeneity(
    bam, chrom, start, end, min_calls, min_cov, hp, tag_name=None, tag_value=None
):

    window = 20
    binsize = 100
    regions = []
    # start = 200
    # end = 300
    for i in range(start, end, window):
        if i + binsize > end:
            regions.append((i, end))
            break
        else:
            regions.append((i, i + binsize))
            if i + binsize == end:
                break

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
        else:
            mean, std, cov = get_region_heterogeneity(
                bam, chrom, region[0], region[1], min_calls, min_cov
            )
            if mean != -1:
                het["x"].append((region[0] + region[1]) / 2)
                het["all"].append(mean)

