from modbamtools.utils import *


def queue_reads(dict_per_read):

    sorted_mod = dict(
        sorted(dict_per_read.items(), key=lambda e: e[1][0], reverse=True)
    )
    i = 0
    out_dict = {}

    for k, v in sorted_mod.items():

        if len(out_dict) == 0:

            out_dict[i] = [(k, v)]
            i -= 1
            continue
        for line, reads in out_dict.items():
            t = 0
            for read in reads:
                if overlaps(v[1], read[1][1]) > -5000:
                    t = 1
            if t == 0:
                out_dict[line].append((k, v))
                break

        if t == 1:
            out_dict[i] = [(k, v)]
            i -= 1
    return out_dict


def record_text_plot(record_start, record_end, start, end):
    """
    plot gene models texts (names)
    """
    if (record_start <= start) & (record_end >= end):
        x = [(start + end) / 2]
        textpos = "top center"
    elif (record_start >= start) & (record_end <= end):
        x = [(record_start + record_end) / 2]
        textpos = "top center"
    elif (record_start >= start) & (record_end >= end):
        x = [record_start]
        textpos = "top right"
    elif (record_start <= start) & (record_end <= end):
        x = [record_end]
        textpos = "top left"

    return (x, textpos)


def parse_gtf(gtf_path, chrom, start, end, record_type="", vertical_spacing=3):
    records = {}

    gtf = pysam.TabixFile(gtf_path, parser=pysam.asGTF())

    if record_type == "":
        record_type = "gene"
    for record in gtf.fetch(chrom, start, end):
        if record.gene_type == "misc_RNA":
            continue
        if record.feature == record_type:

            records[record.gene_name + " (" + record.gene_type + ")"] = [
                record.end - record.start,
                (record.start, record.end),
                record_text_plot(record.start, record.end, start, end),
                record.strand,
            ]

    records = queue_reads(records)
    name_traces = []
    shapes = []
    i = 0
    for row, record_list in records.items():

        for record in record_list:

            if record[1][3] == "+":
                color = "RoyalBlue"
                fill = "LightSkyBlue"
            elif record[1][3] == "-":
                color = "lightseagreen"
                fill = "mediumaquamarine"

            name_traces.append(
                go.Scatter(
                    x=record[1][2][0],
                    y=[(i + 1.5)],
                    text=[record[0]],
                    mode="text",
                    textposition=record[1][2][1],
                    showlegend=False,
                )
            )
            #             shapes.append(go.Scatter(x=[record[1][1][0],record[1][1][0],record[1][1][1],record[1][1][1],record[1][1][0]], y=[i,i-1,i-1,i,i], fill="toself"))
            shapes.append(
                dict(
                    type="rect",
                    x0=record[1][1][0],
                    y0=i,
                    x1=record[1][1][1],
                    y1=(i - 1),
                    line=dict(color=color, width=2),
                    fillcolor=fill,
                )
            )

        i -= vertical_spacing

    ylim = [i + 0.5, 1.5]

    return ylim, name_traces, shapes


def merge_exons(intervals):

    intervals.sort(key=lambda x: x[0])

    merged = []
    for interval in intervals:
        # if the list of merged intervals is empty or if the current
        # interval does not overlap with the previous, simply append it.
        if not merged or merged[-1][1] < interval[0]:
            merged.append(interval)
        else:
            # otherwise, there is overlap, so we merge the current and previous
            # intervals.
            merged[-1][1] = max(merged[-1][1], interval[1])

    return merged


def parse_gtf_exons(gtf_path, chrom, start, end, vertical_spacing=25):
    per_line_height = 60  # px
    gtf = pysam.TabixFile(gtf_path, parser=pysam.asGTF())
    recs = {}
    for record in gtf.fetch(chrom, start, end):
        if "gene_name" in record.attributes and "gene_type" in record.attributes:

            if record.gene_type == "misc_RNA":
                continue
            if record.gene_name + " (" + record.gene_type + ")" not in recs.keys():
                recs[record.gene_name + " (" + record.gene_type + ")"] = {
                    "gene": [],
                    "exons": [],
                }
            if record.feature == "gene":
                recs[record.gene_name + " (" + record.gene_type + ")"]["gene"] = (
                    record.start,
                    record.end,
                    record.strand,
                )
            if record.feature == "exon":
                recs[record.gene_name + " (" + record.gene_type + ")"]["exons"].append(
                    [record.start, record.end]
                )
        else:

            if record.gene_id not in recs.keys():
                if (record.feature == "gene") | (record.feature == "exon"):
                    recs[record.gene_id] = {"gene": [], "exons": []}
            if record.feature == "gene":
                recs[record.gene_id]["gene"] = (record.start, record.end, record.strand)
            if record.feature == "exon":
                recs[record.gene_id]["exons"].append([record.start, record.end])
    # print(recs)
    out = {}
    for k, v in recs.items():
        coo = (v["gene"][0], v["gene"][1])
        length = v["gene"][1] - v["gene"][0]
        strand = v["gene"][2]
        exons = merge_exons(v["exons"])
        out[k] = [
            length,
            coo,
            record_text_plot(coo[0], coo[1], start, end),
            strand,
            exons,
        ]

    records = queue_reads(out)
    name_traces = []
    shapes = []
    i = 0
    row = 0
    for row, record_list in records.items():

        for record in record_list:

            if record[1][3] == "+":
                color = "RoyalBlue"
                fill = "LightSkyBlue"
            elif record[1][3] == "-":
                color = "lightseagreen"
                fill = "mediumaquamarine"

            name_traces.append(
                go.Scatter(
                    x=record[1][2][0],
                    y=[(i + 4)],
                    text=[record[0]],
                    mode="text",
                    textposition=record[1][2][1],
                    showlegend=False,
                    textfont=dict(size=14),
                )
            )

            shapes.append(
                dict(
                    type="line",
                    x0=record[1][1][0],
                    y0=i,
                    x1=record[1][1][1],
                    y1=i,
                    line=dict(color=color, width=2),
                    fillcolor=fill,
                )
            )
            for exon in record[1][4]:
                shapes.append(
                    dict(
                        type="rect",
                        x0=exon[0],
                        y0=i + 2,
                        x1=exon[1],
                        y1=(i - 2),
                        line=dict(color=color, width=2),
                        fillcolor=fill,
                        opacity=1,
                    )
                )

        i -= vertical_spacing

    ylim = [i + vertical_spacing / 2, 20]
    height = (abs(row) + 1) * per_line_height

    return ylim, name_traces, shapes, height


def parse_bed_rectangle(bed_path, chrom, start, end, vertical_spacing=20):

    per_line_height = 50  # px
    shapes = []
    bed = pysam.TabixFile(bed_path)
    records = {}
    for record in bed.fetch(chrom, start, end):

        line = record.strip().split("\t")
        if len(line) == 3:
            coo = [int(i) for i in line[1:3]]
            records["\t".join(line[1:3])] = [
                coo[1] - coo[0],
                coo,
                record_text_plot(coo[0], coo[1], start, end),
                "",
            ]
        else:
            name = line[3]
            coo = [int(i) for i in line[1:3]]
            records["\t".join(line[1:3])] = [
                coo[1] - coo[0],
                coo,
                record_text_plot(coo[0], coo[1], start, end),
                name,
            ]
    records = queue_reads(records)
    name_traces = []
    shapes = []
    i = 0
    row = 0
    for row, record_list in records.items():

        for record in record_list:

            color = "Crimson"
            fill = "Salmon"

            name_traces.append(
                go.Scatter(
                    x=record[1][2][0],
                    y=[(i + 4)],
                    text=record[1][3],
                    mode="text",
                    textposition=record[1][2][1],
                    showlegend=False,
                )
            )

            shapes.append(
                dict(
                    type="rect",
                    x0=record[1][1][0],
                    y0=i + 2,
                    x1=record[1][1][1],
                    y1=(i - 2),
                    line=dict(color=color, width=2),
                    fillcolor=fill,
                    opacity=1,
                )
            )

        i -= vertical_spacing

    ylim = [i, 20]
    height = (abs(row) + 1) * per_line_height

    return ylim, name_traces, shapes, height
