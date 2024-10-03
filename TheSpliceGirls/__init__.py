import os
import pandas as pd
import pyranges as pr
import requests
import pickle as pk
import time
import datetime as dt

import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
from pandas.errors import SettingWithCopyWarning
warnings.simplefilter(action='ignore', category=SettingWithCopyWarning)

def tsg_update_cache(build="hg38",gencode="46"):
    # Set the timeout
    # Can we pickle the GFF for faster loading later?

    requests_timeout = 6000  # seconds

    cache_dir = os.path.join(os.path.expanduser("~"), ".cache", "thesplicegirls-py")
    if not os.path.exists(cache_dir):
        os.makedirs(cache_dir)

    # Define file URLs
    host = "http://hgdownload.soe.ucsc.edu/"
    files = [
        "goldenPath/hg38/database/ucscGenePfam.txt.gz",
        "gbdb/hg38/uniprot/unipDomain.bb",
        "gbdb/hg38/uniprot/unipLocCytopl.bb",
        "gbdb/hg38/uniprot/unipLocExtra.bb",
        "gbdb/hg38/uniprot/unipLocSignal.bb",
        "gbdb/hg38/uniprot/unipLocTransMemb.bb"
    ]

    print("Downloading cache files:", end="")
    for f in files:
        print(f"\n\t{os.path.basename(f)}", end="")
        url = os.path.join(host, f)
        response = requests.get(url, timeout=requests_timeout)
        with open(os.path.join(cache_dir, os.path.basename(f)), 'wb') as file:
            file.write(response.content)

    # Download bigBedToBed binary based on OS
    print("\nConverting cache files:", end="")
    sys_id = "linux.x86_64" if os.name == "posix" else "macOSX.x86_64"
    bin_url = f"https://hgdownload.cse.ucsc.edu/admin/exe/{sys_id}/bigBedToBed"
    response = requests.get(bin_url)
    with open(os.path.join(cache_dir, os.path.basename(bin_url)), 'wb') as file:
        file.write(response.content)

    os.chmod(os.path.join(cache_dir, os.path.basename(bin_url)), 0o755)

    for f in files[1:]:
        f_base = os.path.basename(f).replace(".bb", "")
        print(f"\n\t{f_base}.bb to .bed", end="")
        os.system(f"{os.path.join(cache_dir, os.path.basename(bin_url))} {os.path.join(cache_dir, f_base + '.bb')} {os.path.join(cache_dir, f_base + '.bed')}")

    gencode_file = f"gencode.v{gencode}.annotation.gff3.gz"
    gencode_url = f"http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_{gencode}/{gencode_file}"
    data_file = f"gff_data.v{gencode}.pkl"
    if not os.path.exists(os.path.join(cache_dir, data_file)):
        print(f"\n Downloading Gencode GFF (v{gencode}).... (~3 mins)", end="")
        response = requests.get(gencode_url, timeout=requests_timeout)
        with open(os.path.join(cache_dir, gencode_file), 'wb') as file:
            file.write(response.content)
        print(f"\n\tConverting Gencode GFF (v{gencode}) to pickled pyrange.... (~2 mins)", end="")
        gff_data = pr.read_gff3(os.path.join(cache_dir, gencode_file))
        gff_data['Start'] = gff_data['Start']+1
        gff_data['len'] = gff_data['End']-gff_data['Start']
        gff_data = gff_data[gff_data['len']>0]
        valid_chr = ["chrX", "chrY"] + [f"chr{i}" for i in range(1, 23)]
        gff_data = gff_data[gff_data.Chromosome.isin(valid_chr)]
        with open(os.path.join(cache_dir, data_file), 'wb') as f: 
            pk.dump(gff_data,f)
    print("\nFinished downloading cache version: "+dt.datetime.fromtimestamp(os.path.getmtime(os.path.join(cache_dir, 'unipLocTransMemb.bed'))).strftime("%Y-%m-%d"))

    


def tsg_annotate(splices,build="hg38",gencode="46"):
    data_file = f"gff_data.v{gencode}.pkl"
    cache_dir = os.path.join(os.path.expanduser("~"), ".cache", "thesplicegirls-py")
    if not os.path.exists(os.path.join(cache_dir, data_file)):
        print("No cache found.. ",end="")
        tsg_update_cache(build=build,gencode=gencode)
        print()
    print(f"Starting annotation with cache version: { dt.datetime.fromtimestamp(os.path.getmtime(os.path.join(cache_dir, 'unipLocTransMemb.bed'))).strftime("%Y-%m-%d")}, Gencode v{gencode}", end="")

    print("\n\tLoading GFF data..", end="")
    with open(os.path.join(cache_dir, data_file), 'rb') as f: 
        gff_data = pk.load(f)
    
    valid_chr = ["chrX", "chrY"] + [f"chr{i}" for i in range(1, 23)]
    
    print("\n\tAnnotating splice site information..", end="")
    splices_df = pd.DataFrame(splices, columns=["id"])
    splices_split = splices_df['id'].str.split(":", expand=True)
    splices_df = pd.concat([splices_df, splices_split], axis=1)
    splices_df.columns = ["id", "chr", "start", "end", "clu"]
    splices_df['strand'] = splices_df['clu'].str.split("_").str[2]
    
    if not splices_df['chr'].str.contains("chr").all():
        splices_df['chr'] = "chr" + splices_df['chr']
    
    stranded = splices_df['strand'].isin(["+", "-"]).all()
    if not stranded:
        splices_df.loc[~splices_df['strand'].isin(["+", "-"]), 'strand'] = "*"
    
    # Set stranded coordinates
    splices_df['start_stranded'] = splices_df['start']
    splices_df['end_stranded'] = splices_df['end']
    if stranded:
        splices_df.loc[splices_df['strand'] == "-", 'start_stranded'] = splices_df['end']
        splices_df.loc[splices_df['strand'] == "-", 'end_stranded'] = splices_df['start']

    print("\n\tAnnotating gene body information..", end="")
    
    gff_genes = gff_data[gff_data.Feature == "gene"]
    gff_exons = gff_data[gff_data.Feature == "exon"]
    gff_cds = gff_data[gff_data.Feature == "CDS"]
    gff_3p_utr = gff_data[gff_data.Feature == "three_prime_UTR"]
    gff_5p_utr = gff_data[gff_data.Feature == "five_prime_UTR"]

    splices_ranges = pd.DataFrame(
    {'Chromosome':splices_df["chr"],
     'Start':[int(x)+1 for x in splices_df["start"]],
     'End':[int(x)-1 for x in splices_df["end"]],
     'Strand':splices_df["strand"],
     'title':splices_df["id"]})
    splices_ranges = pr.PyRanges(splices_ranges)
    gff_genes.loc[:,'gene_strand'] = gff_genes['Strand'].values
    
    ov = splices_ranges.join_ranges(gff_genes,report_overlap=True,strand_behavior="ignore")
    ov['splice_length'] = ov.End-ov.Start
    ov['same_strand'] = ov.Strand==ov.gene_strand
    ov['gene_full_match'] = (ov.Overlap==ov.splice_length) & ov.same_strand

    splices_df['gene_name'] = ""
    splices_df['gene_id'] = ""
    splices_df['gene_type'] = ""
    splices_df['gene_opposite_strand'] = False

    splices_df.index = splices_df['id']
    ov['is_protein_coding'] = ov['gene_type']=="protein_coding"
    ov = ov.sort_values('is_protein_coding',ascending=False)

    for sid in splices_df.index:
        splices_df.loc[sid,'gene_name'] = ov[ov['gene_full_match'] & (ov['title']==sid)]['gene_name'].str.cat(sep=",")
        splices_df.loc[sid,'gene_id'] = ov[ov['gene_full_match'] & (ov['title']==sid)]['gene_id'].str.cat(sep=",")
        splices_df.loc[sid,'gene_type'] = ov[ov['gene_full_match'] & (ov['title']==sid)]['gene_type'].str.cat(sep=",")


    if any(splices_df['gene_name']==""):
        for sid in splices_df.index:
            if splices_df.loc[sid,'gene_name']=="":
                splices_df.loc[sid,'gene_name'] = ov[ov['same_strand'] & (ov['title']==sid)]['gene_name'].str.cat(sep=",")
                splices_df.loc[sid,'gene_id'] = ov[ov['same_strand'] & (ov['title']==sid)]['gene_id'].str.cat(sep=",")
                splices_df.loc[sid,'gene_type'] = ov[ov['same_strand'] & (ov['title']==sid)]['gene_type'].str.cat(sep=",")


    if any(splices_df['gene_name']==""):
        for sid in splices_df.index:
            if splices_df.loc[sid,'gene_name']=="":
                splices_df.loc[sid,'gene_name'] = ov[ov['title']==sid]['gene_name'].str.cat(sep=",")
                splices_df.loc[sid,'gene_id'] = ov[ov['title']==sid]['gene_id'].str.cat(sep=",")
                splices_df.loc[sid,'gene_type'] = ov[ov['title']==sid]['gene_type'].str.cat(sep=",")
                splices_df.loc[sid,'gene_opposite_strand'] = True

    print("\n\tAnnotating exon matches..", end="")

    
    splices_df['overlaps_exon'] = ''
    splices_df['overlaps_cds'] = ''
    splices_df['overlaps_3p_utr'] = ''
    splices_df['overlaps_5p_utr'] = ''

    ov_exon = splices_ranges.join_ranges(gff_exons,report_overlap=True,strand_behavior="auto")
    ov_cds = splices_ranges.join_ranges(gff_cds,report_overlap=True,strand_behavior="auto")
    ov_3p = splices_ranges.join_ranges(gff_3p_utr,report_overlap=True,strand_behavior="auto")
    ov_5p = splices_ranges.join_ranges(gff_5p_utr,report_overlap=True,strand_behavior="auto")

    for sid in splices_df.index:
        splices_df.loc[sid,'overlaps_exon'] = sid in ov_exon['title'].values
        splices_df.loc[sid,'overlaps_cds'] = sid in ov_cds['title'].values
        splices_df.loc[sid,'overlaps_3p_utr'] = sid in ov_3p['title'].values
        splices_df.loc[sid,'overlaps_5p_utr'] = sid in ov_5p['title'].values


    splices_df['bp_excised'] = (splices_df['end'].astype('int')-splices_df['start'].astype('int'))-1

    gff_exons.loc[:,'start_stranded'] = gff_exons['Start']
    gff_exons.loc[gff_exons["Strand"]=="-",'start_stranded'] = gff_exons[gff_exons["Strand"]=="-"]['End']
    gff_exons.loc[:,'end_stranded'] = gff_exons['End']
    gff_exons.loc[gff_exons["Strand"]=="-",'end_stranded'] = gff_exons[gff_exons["Strand"]=="-"]['Start']
    gff_exons['chrstart'] = gff_exons['Chromosome'] + ":" + gff_exons['start_stranded'].astype('str')
    gff_exons['chrend'] = gff_exons['Chromosome'] + ":" + gff_exons['end_stranded'].astype('str')
    
    splices_df['p5_exon_overlap'] = [x in gff_exons['chrend'].values for x in (splices_df['chr'] + ":" + splices_df['start_stranded'])]
    splices_df['p3_exon_overlap'] = [x in gff_exons['chrstart'].values for x in (splices_df['chr'] + ":" + splices_df['end_stranded'])]
    splices_df['p5_and_p3_exon_overlap'] = splices_df['p5_exon_overlap'] & splices_df['p3_exon_overlap']
    splices_df['p5_or_p3_exon_overlap'] = splices_df['p5_exon_overlap'] | splices_df['p3_exon_overlap']

    def load_bed(bed_file,bump=0):
        bed_path = os.path.join(cache_dir, bed_file)
        bed_data = pd.read_table(bed_path,header=None,usecols=range(bump,6+bump))
        bed_data.columns = ["Chromosome","Start","End","Name","Score","Strand"]
        bed_data = pr.PyRanges(bed_data)
        bed_data['Start'] = bed_data['Start']+1
        bed_data['len'] = bed_data['End']-bed_data['Start']
        bed_data = bed_data[bed_data['len']>0]
        bed_data = bed_data[bed_data.Chromosome.isin(valid_chr)]
        return(bed_data)

    print("\n\tAnnotating protein localizations..", end="")

    bed_file = f"unipLocTransMemb.bed"
    bed_file = load_bed(bed_file)
    ov_tm = splices_ranges.join_ranges(bed_file,report_overlap=True,strand_behavior="auto")
    
    bed_file = f"unipLocCytopl.bed"
    bed_file = load_bed(bed_file)
    ov_cy = splices_ranges.join_ranges(bed_file,report_overlap=True,strand_behavior="auto")
    
    bed_file = f"unipLocExtra.bed"
    bed_file = load_bed(bed_file)
    ov_ex = splices_ranges.join_ranges(bed_file,report_overlap=True,strand_behavior="auto")
    
    bed_file = f"unipLocSignal.bed"
    bed_file = load_bed(bed_file)
    ov_si = splices_ranges.join_ranges(bed_file,report_overlap=True,strand_behavior="auto")
    
    splices_df['unipLocTransMemb'] = ""
    splices_df['unipLocCytopl'] = ""
    splices_df['unipLocExtra'] = ""
    splices_df['unipLocSignal'] = ""

    for sid in splices_df.index:
        splices_df.loc[sid,'unipLocTransMemb'] = (sid in ov_tm['title'].values) & splices_df.loc[sid,'overlaps_cds']
        splices_df.loc[sid,'unipLocCytopl'] = (sid in ov_cy['title'].values) & splices_df.loc[sid,'overlaps_cds']
        splices_df.loc[sid,'unipLocExtra'] = (sid in ov_ex['title'].values) & splices_df.loc[sid,'overlaps_cds']
        splices_df.loc[sid,'unipLocSignal'] = (sid in ov_si['title'].values) & splices_df.loc[sid,'overlaps_cds']


    print("\n\tAnnotating pfam and uniprot domains..", end="")

    splices_df['unipDomain'] = ""
    splices_df['pfamDomain'] = ""
    
    bed_file = f"unipDomain.bed"
    bed_file = load_bed(bed_file)
    ov_unip = splices_ranges.join_ranges(bed_file,report_overlap=True,strand_behavior="auto")
    
    bed_file = f"ucscGenePfam.txt.gz"
    bed_file = load_bed(bed_file,bump=1)
    ov_pfam = splices_ranges.join_ranges(bed_file,report_overlap=True,strand_behavior="auto")
    
    for sid in splices_df.index:
        splices_df.loc[sid,'unipDomain'] = ov_unip[(ov_unip['title']==sid)]['Name'].str.cat(sep=",")
        splices_df.loc[sid,'pfamDomain'] = ov_pfam[(ov_pfam['title']==sid)]['Name'].str.cat(sep=",")
    
    splices_df.loc[splices_df["overlaps_cds"]==False,"unipDomain"] = ""
    splices_df.loc[splices_df["overlaps_cds"]==False,"pfamDomain"] = ""

    print("\nDone!", end="")

    return(splices_df)




    