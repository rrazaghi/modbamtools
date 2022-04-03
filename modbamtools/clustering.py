from functools import partial
from modbamtools.utils import *
import hdbscan


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
