from modbamtools.utils import *

def queue_reads(dict_per_read):
    
    sorted_mod = dict(sorted(dict_per_read.items(), key=lambda e: e[1][0], reverse=True))
    i = 0
    out_dict = {}
    
    for k,v in sorted_mod.items():

        if len(out_dict) == 0:
            
            out_dict[i] = [(k,v)]
            i -= 1
            continue
        for line,reads in out_dict.items():
            t =0
            for read in reads:
                if overlaps(v[1], read[1][1]) > -5000:
                    t = 1
            if t == 0:
                out_dict[line].append((k,v))
                break
                
        if t == 1:
            out_dict[i] = [(k,v)]
            i -= 1          
    return out_dict

def record_text_plot(record,start,end):
    """
    plot gene models texts (names)
    """
    if (record.start > start) & (record.end > end):
            x = [record.start]
            textpos = "bottom right"
    elif (record.start < start) & (record.end < end):
            x = [record.end]
            textpos = "bottom left"
    elif (record.start >= start) & (record.end <= end):
            x = [(record.start + record.end)/2]
            textpos = "bottom center"
    elif (record.start < start) & (record.end > end):
            x = [(start + end)/2]
            textpos = "bottom center"
            
    return (x,textpos)
    

def parse_gtf(gtf_path,chrom,start,end, record_type="", vertical_spacing=3):
    records = {}
    
    gtf = pysam.TabixFile(gtf_path, parser = pysam.asGTF())

    
    if record_type == "":
        record_type = "gene"
    for record in gtf.fetch(chrom, start, end):
        if record.gene_type == "misc_RNA":
            continue
        if record.feature == record_type:
            
            records[record.gene_name + " (" + record.gene_type + ")"] = [record.end - record.start,
                              (record.start,record.end),
                              record_text_plot(record,start,end),record.strand]
            
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
                
            name_traces.append(go.Scatter(
                    x=record[1][2][0],
                    y=[(i + 1.5)],
                    text=[record[0]],
                    mode="text",
                    textposition=record[1][2][1],showlegend=False
                ))
#             shapes.append(go.Scatter(x=[record[1][1][0],record[1][1][0],record[1][1][1],record[1][1][1],record[1][1][0]], y=[i,i-1,i-1,i,i], fill="toself"))
            shapes.append(dict(type="rect",
                    x0=record[1][1][0], y0=i, x1=record[1][1][1], y1=(i - 1),
                    line=dict(color=color, width=2),
                    fillcolor=fill
                ))
        
        i -= vertical_spacing
     
    ylim = [i + 0.5,1.5]
        
        
    return ylim, name_traces, shapes
    

def parse_bed(bed_path,chrom,start,end):
    shapes = []
    bed = pysam.TabixFile(bed_path)
    records = set()   
    for record in bed.fetch(chrom, start, end):
        
        line = record.split("\t")
        
#         try:
#             abc_score = float(line[7])
#         except ValueError:
#             continue
            
#         if float(line[20]) < 0.01:
#             continue
        coo = line[0:3]
        coo = "\t".join(coo)
        records.add(coo)
        
    for record in records:
        coo = record.split("\t")
# #         print(coo[4])
#         if coo[4] == "promoter":
#             color = "Black"
#             fill = "Grey"
# #         elif line[4] == "FALSE":
#         else:
        color = "Crimson"
        fill = "Salmon"
#         shapes.append(go.Scatter(x=[coo[1],coo[1],coo[2],coo[2],coo[1]], y=[0,1,1,0,0], fill="toself"))
        shapes.append(dict(type="rect",
                    x0=coo[1], y0=0, x1=coo[2], y1=1,
                    line=dict(color=color, width=2),
                    fillcolor=fill
                ))
            
    
     
    ylim = [-1,2]
        
        
    return ylim, shapes
    