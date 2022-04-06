from modbamtools.utils import *


def get_read_stats_horizontal(read, start, end, min_calls):

    mod_calls = []
    for pos_mod in read.mod_sites:
        qname, rpos, qpos, strand, mod_strand, cbase, mbase, score = pos_mod
        if (rpos >= start) & (rpos <= end):
            call = binerize_mod_call(score)
            mod_calls.append(call)
    if len(mod_calls) < min_calls:
        percent_meth = -1
    else:
        percent_meth = mod_calls.count(1) * 100 / len(mod_calls)

    return percent_meth


def get_region_stats_horizontal(
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
            percent_meth = get_read_stats_horizontal(read, start, end, min_calls)
            if percent_meth == -1:
                continue
            read_stats.append(percent_meth)
    if len(read_stats) != 0:
        read_stats = np.array(read_stats)
        mean = np.mean(read_stats)
        std = np.std(read_stats)
        cov = len(read_stats)
        return mean, std, cov
    else:
        return -1, -1, -1


def get_region_stats_hap_horizontal(bam, min_calls, min_cov, hp, bed_line):

    chrom = bed_line[0]
    start = int(bed_line[1])
    end = int(bed_line[2])

    if hp:

        mean, std, cov = get_region_stats_horizontal(
            bam, chrom, start, end, min_calls, min_cov
        )
        hp1_mean, hp1_std, hp1_cov = get_region_stats_horizontal(
            bam, chrom, start, end, min_calls, min_cov, tag_name="HP", tag_value=1
        )
        hp2_mean, hp2_std, hp2_cov = get_region_stats_horizontal(
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

        mean, std, cov = get_region_stats_horizontal(
            bam, chrom, start, end, min_calls, min_cov
        )
        bed_line.extend(map(str, [mean, std, cov],))

    return bed_line

