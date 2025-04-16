import os
import pandas as pd
import pyranges as pr
import requests
import pickle as pk
import time
import datetime as dt
import numpy as np

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

    # Should only add chr to those without
    # So first remove?
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
    
    ov = splices_ranges.join_ranges(gff_genes, report_overlap_column="Overlap", strand_behavior="ignore")
    ov['splice_length'] = ov.End - ov.Start
    ov['same_strand'] = ov.Strand == ov.gene_strand
    ov['gene_full_match'] = (ov.Overlap == ov.splice_length) & ov.same_strand
    ov['is_protein_coding'] = ov['gene_type'] == "protein_coding"
    # Sort by 'is_protein_coding' (protein coding genes take priority)
    ov = ov.sort_values('is_protein_coding', ascending=False)
    # Initialize gene-related columns in splices_df
    splices_df['gene_name'] = None
    splices_df['gene_id'] = None
    splices_df['gene_type'] = None
    splices_df['gene_full_match'] = None
    splices_df['gene_opposite_strand'] = None

    splices_df.index = splices_df['id']
    # Step 1: Apply gene_full_match to the entire DataFrame
    gene_full_match_df = ov[ov['gene_full_match']].groupby('title').agg({
        'gene_name': lambda x: ",".join(x),
        'gene_id': lambda x: ",".join(x),
        'gene_type': lambda x: ",".join(x),
    }).reset_index()
    gene_full_match_df['gene_full_match'] = True
    gene_full_match_df['gene_opposite_strand'] = False
    splices_df.update(gene_full_match_df.set_index('title'))
    
    # Step 2: For missing gene_full_match, fallback to same_strand
    same_strand_df = ov[ov['same_strand']].groupby('title').agg({
        'gene_name': lambda x: ",".join(x),
        'gene_id': lambda x: ",".join(x),
        'gene_type': lambda x: ",".join(x)
    }).reset_index()
    same_strand_df['gene_full_match'] = False
    same_strand_df['gene_opposite_strand'] = False
    # Fill missing gene info with same_strand matches
    splices_df.update(same_strand_df.set_index('title'),overwrite=False)
    
    # Step 3: For missing gene info, fallback to any overlap
    any_overlap_df = ov.groupby('title').agg({
        'gene_name': lambda x: ",".join(x),
        'gene_id': lambda x: ",".join(x),
        'gene_type': lambda x: ",".join(x)
    }).reset_index()
    any_overlap_df['gene_full_match'] = False
    any_overlap_df['gene_opposite_strand'] = True
    
    # Update remaining missing values
    splices_df.update(any_overlap_df.set_index('title'),overwrite=False)
    splices_df.fillna({'gene_name': "", 'gene_id': "", 'gene_type': "", 'gene_full_match': False, 'gene_opposite_strand': False}, inplace=True)

    print("\n\tAnnotating exon matches..", end="")

    # Collapse exon ranges before to speed up?
    # How much time overlap vs the loop layer?
    ov_exon = set(splices_ranges.join_ranges(gff_exons,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)

    # Shoulda also make another category for intra-cds, i.e. excises just a part of it
    ov_cds = set(splices_ranges.join_ranges(gff_cds,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)
    ov_cds_full = splices_ranges.join_ranges(gff_cds,report_overlap_column="Overlap",strand_behavior="auto")
    ov_cds_full['bp_excised'] = (ov_cds_full['End']-ov_cds_full['Start']) # Deliberately not +1 to match with pr "Overlap"
    ov_cds_intra = set(ov_cds_full[ov_cds_full['Overlap']==ov_cds_full['bp_excised']]['title'].values)

    ov_3p = set(splices_ranges.join_ranges(gff_3p_utr,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)
    ov_5p = set(splices_ranges.join_ranges(gff_5p_utr,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)

    splices_df['overlaps_exon'] = splices_df.index.isin(ov_exon)
    splices_df['overlaps_cds'] = splices_df.index.isin(ov_cds)
    splices_df['intra_cds_splice'] = splices_df.index.isin(ov_cds_intra)
    splices_df['overlaps_3p_utr'] = splices_df.index.isin(ov_3p)
    splices_df['overlaps_5p_utr'] = splices_df.index.isin(ov_5p)


    splices_df['bp_excised'] = (splices_df['end'].astype('int')-splices_df['start'].astype('int'))-1

    gff_exons['start_stranded'] = np.where(gff_exons["Strand"] == "-", gff_exons['End'], gff_exons['Start'])
    gff_exons['end_stranded'] = np.where(gff_exons["Strand"] == "-", gff_exons['Start'], gff_exons['End'])
    
    gff_exons['chrstart'] = gff_exons['Chromosome'] + ":" + gff_exons['start_stranded'].astype(str)
    gff_exons['chrend'] = gff_exons['Chromosome'] + ":" + gff_exons['end_stranded'].astype(str)

    chrend_set = set(gff_exons['chrend'])
    chrstart_set = set(gff_exons['chrstart'])

    chr_plus_start_stranded = splices_df['chr'] + ":" + splices_df['start_stranded'].astype(str)
    chr_plus_end_stranded = splices_df['chr'] + ":" + splices_df['end_stranded'].astype(str)
    
    splices_df['p5_exon_overlap'] = chr_plus_start_stranded.isin(chrend_set)
    splices_df['p3_exon_overlap'] = chr_plus_end_stranded.isin(chrstart_set)
    
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
    ov_tm = set(splices_ranges.join_ranges(bed_file,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)
    
    bed_file = f"unipLocCytopl.bed"
    bed_file = load_bed(bed_file)
    ov_cy = set(splices_ranges.join_ranges(bed_file,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)
    
    bed_file = f"unipLocExtra.bed"
    bed_file = load_bed(bed_file)
    ov_ex = set(splices_ranges.join_ranges(bed_file,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)
    
    bed_file = f"unipLocSignal.bed"
    bed_file = load_bed(bed_file)
    ov_si = set(splices_ranges.join_ranges(bed_file,report_overlap_column="Overlap",strand_behavior="auto")['title'].values)

    splices_df['unipLocTransMemb'] = splices_df.index.isin(ov_tm) & splices_df['overlaps_cds']
    splices_df['unipLocCytopl'] = splices_df.index.isin(ov_cy) & splices_df['overlaps_cds']
    splices_df['unipLocExtra'] = splices_df.index.isin(ov_ex) & splices_df['overlaps_cds']
    splices_df['unipLocSignal'] = splices_df.index.isin(ov_si) & splices_df['overlaps_cds']


    print("\n\tAnnotating pfam and uniprot domains..", end="")
    splices_df.index = splices_df['id']
    
    bed_file = f"unipDomain.bed"
    bed_file = load_bed(bed_file)
    ov_unip = splices_ranges.join_ranges(bed_file,report_overlap_column="Overlap",strand_behavior="auto")
    ov_unip = ov_unip.groupby('title')['Name'].apply(lambda x: ",".join(x)).reset_index()
    splices_df = splices_df.merge(ov_unip, left_index=True, right_on='title', how='left')
    splices_df['unipDomain'] = splices_df['Name'].fillna("")
    splices_df.drop(columns=['title', 'Name'], inplace=True)
    splices_df.index = splices_df['id']
    
    bed_file = f"ucscGenePfam.txt.gz"
    bed_file = load_bed(bed_file,bump=1)
    ov_pfam = splices_ranges.join_ranges(bed_file,report_overlap_column="Overlap",strand_behavior="auto")
    ov_pfam = ov_pfam.groupby('title')['Name'].apply(lambda x: ",".join(x)).reset_index()
    splices_df = splices_df.merge(ov_pfam, left_index=True, right_on='title', how='left')
    splices_df['pfamDomain'] = splices_df['Name'].fillna("")
    splices_df.drop(columns=['title', 'Name'], inplace=True)
    splices_df.index = splices_df['id']

    splices_df.loc[splices_df["overlaps_cds"]==False,"unipDomain"] = ""
    splices_df.loc[splices_df["overlaps_cds"]==False,"pfamDomain"] = ""

    print("\nDone!", end="")

    return(splices_df)

# A function to sum up to the cluster level? I.e. collect annots for any/all