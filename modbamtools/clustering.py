from functools import partial
from modbamtools.utils import *
import hdbscan
import collections


def cluster_region(bed_line, bam):

    chrom = bed_line[0]
    start = int(bed_line[1])
    end = int(bed_line[2])
    region_length = end - start

    dict_per_read_mod = process_bam(bam, chrom, start, end)
    reads_dict = {}
    filtered_reads = 0
    accepted_reads = 0
    for rname, arr in dict_per_read_mod.items():
        if arr[0] / region_length > 0.9:
            #             print(len(arr[2]))
            if len(arr[2]) >= 5:
                accepted_reads += 1
                reads_dict[rname] = arr[2]
        else:
            filtered_reads += 1
    df = pd.DataFrame.from_dict(reads_dict, orient="index")
    df = df[df.columns[df.columns.isin(range(start, end))]]
    df = df.reindex(sorted(df.columns), axis=1)
    df.fillna(0, inplace=True)
    if (df.shape[0] > 15) & (df.shape[1] > 5):
        #     if accepted_reads > 15:

        #         df = pd.DataFrame.from_dict(reads_dict, orient='index')
        #         df = df[df.columns[df.columns.isin(range(start,end))]]
        #         df = df.reindex(sorted(df.columns), axis=1)
        #         df.fillna(0, inplace=True)
        clusterer = hdbscan.HDBSCAN(
            metric="hamming", min_cluster_size=10, min_samples=5
        )
        clusterer.fit(df.to_numpy())
        clusters = {}

        cl_labels = list(set(clusterer.labels_))
        num_cl = clusterer.labels_.max() + 1
        bed_line.extend([str(num_cl), str(accepted_reads)])
    else:
        bed_line.extend([str(0), str(accepted_reads)])
    return bed_line


def cluster2dicts(bams, chrom, start, end, min_cov=0.9):
    dict_per_read_mod = {}
    for bam in bams:
        dict_mod = process_bam(bam, chrom, start, end)
        dict_per_read_mod = {**dict_per_read_mod, **dict_mod}
    reads_dict = {}
    for rname, arr in dict_per_read_mod.items():
        if arr[0] / (end - start) > min_cov:
            if len(arr[2]) >= 5:
                reads_dict[rname] = arr[2]

    df = pd.DataFrame.from_dict(reads_dict, orient="index")
    df = df[df.columns[df.columns.isin(range(start, end))]]
    df = df.reindex(sorted(df.columns), axis=1)
    df.fillna(0, inplace=True)
    # if (df.shape[0] > 15) & (df.shape[1] > 5):
    clusterer = hdbscan.HDBSCAN(metric="hamming", min_cluster_size=10, min_samples=5)
    clusterer.fit(df.to_numpy())
    cl_labels = sorted(list(clusterer.labels_))
    print(cl_labels)
    out = {}
    for read_id, cluster in zip(list(df.index), clusterer.labels_):
        new_id = cluster + 1
        if new_id == 0:
            continue
        if new_id not in out.keys():
            out[new_id] = {}
        out[new_id][read_id] = dict_per_read_mod[read_id]

    out = collections.OrderedDict(sorted(out.items()))
    names = ["Cluster " + str(cl) for cl in out.keys()]
    dicts = [v for v in out.values()]

    return dicts, names
