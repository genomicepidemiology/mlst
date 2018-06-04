#!/usr/bin/env python3

import os, sys, re, time, pprint
import argparse, subprocess

from cgecore.alignment import extended_cigar
from CGE.cgefinder import CGEFinder
from CGE.blaster.blaster import Blaster
import json, gzip
from tabulate import tabulate

def get_read_filename(infiles):
    ''' Infiles must be a list with 1 or 2 input files.
        Removes path from given string and removes extensions:
        .fq .fastq .gz and .trim
        extract the common sample name i 2 files are given.
    '''
    # Remove common fastq extensions
    seq_path = infiles[0]
    seq_file = os.path.basename(seq_path)
    seq_file = seq_file.replace(".fq", "")
    seq_file = seq_file.replace(".fastq", "")
    seq_file = seq_file.replace(".gz", "")
    seq_file = seq_file.replace(".trim", "")
    if len(infiles) == 1:
        return seq_file.rstrip()

    # If two files are given get the common sample name
    sample_name = ""
    seq_file_2 = os.path.basename(infiles[1])
    for i in range(len(seq_file)):
        if seq_file_2[i] == seq_file[i]:
            sample_name += seq_file[i]
        else: 
            break
    if sample_name == "":
        sys.error("Input error: sample names of input files, {} and {}, \
                   does not share a common sample name. If these files \
                   are paired end reads from the same sample, please rename \
                   them with a common sample name (e.g. 's22_R1.fq', 's22_R2.fq') \
                   or input them seperately.".format(infiles[0], infiles[1]))

    return sample_name.rstrip("-").rstrip("_")

def get_file_format(input_files):
    """
    Takes all input files and checks their first character to assess
    the file format. Returns one of the following strings; fasta, fastq, 
    other or mixed. fasta and fastq indicates that all input files are 
    of the same format, either fasta or fastq. other indiates that all
    files are not fasta nor fastq files. mixed indicates that the inputfiles
    are a mix of different file formats.
    """

    # Open all input files and get the first character
    file_format = []
    for infile in input_files:
        if infile[-3:] == ".gz":
            f = gzip.open(infile, "rb")
            fst_char = f.read(1);
        else:
            f = open(infile, "rb")
            fst_char = f.read(1);
        f.close()
        # Assess the first character
        if fst_char == b"@":
            file_format.append("fastq")
        elif fst_char == b">":
            file_format.append("fasta")
        else:
            invalid_files.append("other")
    if len(set(file_format)) != 1:
        return "mixed"
    return ",".join(set(file_format))

def import_profile(database, species, loci_list):
    """Import all possible allele profiles with corresponding st's
    for the species into a dict. The profiles are stored in a dict 
    of dicts, to easily look up what st types are accosiated with 
    a specific allele number of each loci. 
    """
    # Open allele profile file from database
    profile_file = open("{0}/{1}/{1}.tsv".format(database, species), "r")
    profile_header = profile_file.readline().strip().split("\t")[1:len(loci_list)+1]
    
    # Create dict for looking up st-types with locus/allele combinations
    st_profiles = {}
    # For each locus initate make an inner dict to store allele and st's
    for locus in loci_list:
        st_profiles[locus] = {} 

    # Fill inner dict with allele no as key and st-types seen with the allele as value   
    for line in profile_file:
        profile = line.strip().split("\t")
        st_name = profile[0]
        allele_list = profile[1:len(loci_list)+1]

        # Go through all allele profiles. Save locus-allele combination with the st-type
        for i in range(len(allele_list)):
            allele = allele_list[i]
            locus = profile_header[i]
            if allele in st_profiles[locus]:
                st_profiles[locus][allele] += [st_name]
            else:
                st_profiles[locus][allele] = [st_name]
    profile_file.close()

    return st_profiles

def st_typing(st_profiles, allele_matches, loci_list):
    """
    Takes the path to a pickled dictionary, the inp list of the allele 
    number that each loci has been assigned, and an output file string
    where the found st type and similaity is written into it.  
    """

    # Find best ST type for all allele profiles
    st_output = ""
    note = ""

    # First line contains matrix column headers, which are the specific loci
    st_hits = []
    st_marks = []
    note = ""

    # Check the quality of the alle hits
    for locus in allele_matches:
        allele = allele_matches[locus]["allele"]
 
        # Check if allele is marked as a non-perfect match. Save mark and write note.
        if "?*" in allele:
            note += "* {}: Completely imperfect hit, ST can not be trusted!\n".format(locus)
            st_marks = ["?","*"]
        elif "?" in allele:
            note += "* {}: Uncertain hit, ST can not be trusted.\n".format(locus)
            st_marks.append("?")
        elif "*" in allele:
            note += "* {}: Novel allele, ST indicates nearest ST.\n".format(locus)
            st_marks.append("*")

        # Remove mark from allele so it can be used to look up nearest st types
        allele = allele.rstrip("*?")

        # Get all st's that have the alleles in it's allele profile
        st_hits += st_profiles[locus].get(allele, ["None"])

    # Save allele marks to be transfered to the ST
    st_mark = "".join(set(st_marks))
    # Add marks information to notes
    if len(st_mark) == 2:
        note = st_mark + " alleles with less than 100% identity and 100% coverages found\n" + note
    elif st_mark == "*":
        note = st_mark + " alleles with less than 100% identity found\n" + note
    elif st_mark == "?":
        note = st_mark + " alleles with less than 100% coverages found\n" + note

    # Find most frequent st in st_hits
    st_hits_counter = {}
    max_count = 0
    best_hit = ""
    for hit in st_hits:
        if hit is not "None":
            if hit in st_hits_counter:
                st_hits_counter[hit] += 1
            else:
               st_hits_counter[hit] = 1
            if max_count < st_hits_counter[hit]:
                max_count = st_hits_counter[hit]
                best_hit = hit

    # Check if allele profile match found st 100 %
    similarity = round(float(max_count)/len(loci_list)*100, 2)

    if similarity != 100:
        st = "Unknown"
        nearest_sts = []
        # If st is not perfect find nearest st's
        for st_hit, allele_score in sorted(st_hits_counter.items(), key=lambda x: x[1], reverse=True):
            if allele_score < max_count:
                break
            nearest_sts.append(st_hit)
        nearest_sts = ",".join(nearest_sts) #+ st_mark
    else:
        # allele profile has a perfect ST hit but the st marks given to the alleles might indicate imperfect hits 
        st = best_hit + st_mark
        nearest_sts = ""

    return st, note, nearest_sts

def make_aln(species, file_handle, allele_matches, query_aligns, homol_aligns, sbjct_aligns):
    for locus, locus_info in allele_matches.items():        
        allele_name = locus_info["allele_name"]
        if allele_name == "No hit found":
            continue
        hit_name = locus_info["hit_name"]

        query_seq = query_aligns[species][hit_name]
        homol_seq = homol_aligns[species][hit_name]
        sbjct_seq = sbjct_aligns[species][hit_name]

        file_handle.write("# {}".format(allele_name) + "\n")
        for i in range(0,len(query_seq),60):
            file_handle.write("%-10s\t%s\n"%("template:", sbjct_seq[i:i+60]))
            file_handle.write("%-10s\t%s\n"%("", homol_seq[i:i+60]))
            file_handle.write("%-10s\t%s\n\n"%("query:", query_seq[i:i+60]))

def text_table(headers, rows, empty_replace='-'):
   ''' Create text table
   
   USAGE:
      >>> from tabulate import tabulate
      >>> headers = ['A','B']
      >>> rows = [[1,2],[3,4]]
      >>> print(text_table(headers, rows))
      **********
        A    B
      **********
        1    2
        3    4
      ==========
   '''
   # Replace empty cells with placeholder
   rows = map(lambda row: map(lambda x: x if x else empty_replace, row), rows)
   # Create table
   table = tabulate(rows, headers, tablefmt='simple').split('\n')
   # Prepare title injection
   width = len(table[0])
   # Switch horisontal line
   table[1] = '*'*(width+2)
   # Update table with title
   table = ("%s\n"*3)%('*'*(width+2), '\n'.join(table), '='*(width+2))
   return table


#if __name__ == "__main__":
parser = argparse.ArgumentParser(description="")
# Posotional arguments
parser.add_argument("-i", "--infile",
                    help="FASTA or FASTQ files to do MLST on.",
                    nargs="+", required=True)
# Optional arguments
parser.add_argument("-o", "--outdir",
                    help="Output directory.",
                    default="MLST_out.txt")
parser.add_argument("-s", "--species",
                    help="species database used for MLST prediction", required=True)
parser.add_argument("-p", "--database",
                    help="Directory containing the databases.", required=True)#,
                    #default="/home/data1/services/MLST/database_repository/mlst-2.0_db/")
parser.add_argument("-t", "--tmp_dir",
                    help="Temporary directory for storage of the results\
                          from the external software.",
                    default="tmp_MLST")
parser.add_argument("-mp", "--method_path",
                    help="Path to the method to use (kma or blastn)\
                          if assembled contigs are inputted the path to executable blastn should be given,\
                          if fastq files are given path to executable kma should be given",
                    default="blastn")
parser.add_argument("-x", "--extented_output",
                    help="Give extented output with allignment files, template and query hits in fasta and\
                          a tab seperated file with allele profile results", action="store_true")
#parser.add_argument("-c", "--coverage",
#                    help="Minimum template coverage required", default = 0.6)
#parser.add_argument("-i", "--identity",
#                    help="Minimum template identity required", default = 0.9)
args = parser.parse_args()

#TODO what are the clonal complex data used for??

# TODO error handling
infile = args.infile
# Check that outdir is an existing dir...
outdir = os.path.abspath(args.outdir)
species = args.species
database = os.path.abspath(args.database)
tmp_dir = os.path.abspath(args.tmp_dir)
# Check if method path is executable
method_path = args.method_path
extented_output = args.extented_output

min_cov = 0.6   # args.coverage
threshold = 0.9 # args.identity

# Check file format (fasta, fastq or other format)
file_format = get_file_format(infile)

db_path = "{}/{}/".format(database, species)

config_file = open(database + "/config","r")

# Get loci list from config file
for line in config_file:
    line = line.split("\t")
    if line[0] == species:
        organism = line[1]
        loci_list = line[2].strip().split(",")

# Call appropriate method (kma or blastn) based on file format 
if file_format == "fastq":
    # Check the number of files
    if len(infile) == 1:
        infile_1 = infile[0]
        infile_2 = None
    elif len(infile) == 2:
        infile_1 = infile[0]
        infile_2 = infile[1]
    else:
        sys.exit("Only 2 input file accepted for raw read data,\
                  if data from more runs is avaliable for the same\
                  sample, please concatinate the reads into two files")
    
    sample_name = get_read_filename(infile)
    method = "kma"

    # Call KMA
    method_obj = CGEFinder.kma(infile_1, outdir, [species], db_path, min_cov=min_cov,
                                threshold=threshold, kma_path=method_path, sample_name=sample_name,
                                inputfile_2=infile_2, kma_mrs=0.75, kma_gapopen=-5,
                                kma_gapextend=-1, kma_penalty=-3, kma_reward=1)

elif file_format == "fasta":
    # Assert that only one fasta file is inputted
    assert len(infile) == 1, "Only one input file accepted for assembled data"
    infile = infile[0]
    method = "blast"

    # Call BLASTn
    method_obj = Blaster(infile, [species], db_path, tmp_dir, 
                         min_cov, threshold, method_path, cut_off=False)
else:
    sys.exit("Input file must be fastq or fasta format, not "+ file_format)

results      = method_obj.results
query_aligns = method_obj.gene_align_query
homol_aligns = method_obj.gene_align_homo
sbjct_aligns = method_obj.gene_align_sbjct

# Check that the results dict is not empty
if results[species] == "No hit found":
    sys.exit("Non of the MLST loci was found in the input data. \
              Make sure that the correct MLST scheme was chosen, \
              or check the quality of the sequencing/assembled data")

allele_matches = {}

# Get the found allele profile contained in the results dict
for hit, locus_hit in results[species].items():

    # Get allele number for locus
    allele_name = locus_hit["sbjct_header"]
    allele_obj  = re.search("(\w+)[_|-](\d+$)", allele_name)

    # Get variable to later storage in the results dict
    locus     = allele_obj.group(1)
    allele    = allele_obj.group(2)
    coverage  = float(locus_hit["perc_coverage"])
    identity  = float(locus_hit["perc_ident"])
    score     = float(locus_hit["cal_score"])
    gaps      = int(locus_hit["gaps"])    
    align_len   = locus_hit["HSP_length"]
    sbj_len   = int(locus_hit["sbjct_length"])
    sbjct_seq = locus_hit["sbjct_string"]
    query_seq = locus_hit["query_string"] 
    homol_seq = locus_hit["homo_string"]
    cigar     = extended_cigar(sbjct_aligns[species][hit], query_aligns[species][hit]) 

    # Check for perfect hits
    if coverage == 100 and identity == 100:
        # Overwrite alleles already saved, save the perfect match and break to go to next locus
        allele_matches[locus] = {"score":score, "allele":allele, "coverage":coverage,
                                 "identity":identity, "match_priority": 1, "align_len":align_len,
                                 "gaps":gaps, "sbj_len":sbj_len, "allele_name":allele_name,
                                 "sbjct_seq":sbjct_seq, "query_seq":query_seq, "homol_seq":homol_seq, 
                                 "hit_name":hit, "cigar":cigar} 
    else:
        # If no hit has yet been stored initialize dict variables that are looked up below
        if locus not in allele_matches:
            allele_matches[locus] = {"score":0, "match_priority": 4}

        # We weight full coverage higher than perfect identity match
        if coverage == 100 and identity != 100:
            # Check that better (higher prioritized) 100% coverage hit has not been stored yet
            if allele_matches[locus]["match_priority"] > 2 or (allele_matches[locus]["match_priority"] == 2 and score > allele_matches[locus]["score"]):
                allele_matches[locus] = {"score":score, "allele":allele+"*", "coverage":coverage,
                                         "identity":identity, "match_priority": 2, "align_len":align_len,
                                         "gaps":gaps, "sbj_len":sbj_len, "allele_name":allele_name,
                                         "sbjct_seq":sbjct_seq, "query_seq":query_seq, "homol_seq":homol_seq,
                                         "hit_name":hit, "cigar":cigar}
        elif coverage != 100 and identity == 100:
            # Check that higher prioritized hit was not already stored
            if allele_matches[locus]["match_priority"] > 3 or (allele_matches[locus]["match_priority"] == 3 and score > allele_matches[locus]["score"]):
                allele_matches[locus] = {"score":score, "allele":allele + "?", "coverage":coverage,
                                         "identity":identity, "match_priority": 3, "align_len":align_len,
                                         "gaps":gaps, "sbj_len":sbj_len, "allele_name":allele_name,
                                         "sbjct_seq":sbjct_seq, "query_seq":query_seq, "homol_seq":homol_seq,
                                         "hit_name":hit, "cigar":cigar}
        else: # coverage != 100 and identity != 100:
            if allele_matches[locus]["match_priority"] == 4 and score > allele_matches[locus]["score"]:
                allele_matches[locus] = {"score":score, "allele":allele + "?*", "coverage":coverage,
                                         "identity":identity, "match_priority": 4, "align_len":align_len,
                                         "gaps":gaps, "sbj_len":sbj_len, "allele_name":allele_name,
                                         "sbjct_seq":sbjct_seq, "query_seq":query_seq, "homol_seq":homol_seq,
                                         "hit_name":hit, "cigar":cigar}
for locus in loci_list:
    if locus not in allele_matches:
        allele_matches[locus] = {"identity":"", "coverage":"", "allele":"", "allele_name":"No hit found", "align_len":"", "gaps":"", "sbj_len":""}

# Import all possible st profiles into dict
st_profiles = import_profile(database, species,loci_list)

# Find st or neatest sts
st, note, nearest_sts = st_typing(st_profiles, allele_matches, loci_list)
    
# Get run info for JSON file
service = os.path.basename(__file__).replace(".py", "")
date = time.strftime("%d.%m.%Y")
time = time.strftime("%H:%M:%S")

# TODO find a system to show the database and service version using git

# Make JSON output file
data = {service:{}}
allele_results = {}
for locus, locus_info in allele_matches.items():
    allele_results[locus] = {}
    for (key, value) in locus_info.items():
        if key in ["identity", "coverage", "allele", "allele_name", "align_len", "gaps", "sbj_len"]:
            allele_results[locus][key] = value

userinput = {"filename":args.infile, "species":args.species, "organism":organism,"file_format":file_format}
run_info = {"date":date, "time":time}#, "database":{"remote_db":remote_db, "last_commit_hash":head_hash}}
server_results = {"sequence_type":st, "allele_profile": allele_results,
                             "nearest_sts":nearest_sts, "notes":note}

data[service]["user_input"] = userinput
data[service]["run_info"] = run_info
data[service]["results"] = server_results

pprint.pprint(data)

# Save json output
result_file = "{}/data.json".format(outdir) 
with open(result_file, "w") as outfile:  
    json.dump(data, outfile)

if extented_output:
    # Define extented output 
    table_filename  = "{}/results_tab.tsv".format(outdir)
    query_filename  = "{}/Hit_in_genome_seq.fsa".format(outdir)
    sbjct_filename  = "{}/MLST_allele_seq.fsa".format(outdir)
    result_filename = "{}/results.txt".format(outdir)
    table_file  = open(table_filename, "w")
    query_file  = open(query_filename, "w")
    sbjct_file  = open(sbjct_filename, "w")
    result_file = open(result_filename, "w")

    # Make results file
    result_file.write("{0} Results\n\n{0} Profile: {1}\n\n".format(service, species))
    result_file.write("Organism: {}\n\nSequence Type: {}\n".format(organism, st))
    # If ST is unknown report nearest ST
    if st == "Unknown" and nearest_sts != "":
        if len(nearest_sts.split(",")) == 1:
            result_file.write("Nearest ST: {}\n".format(nearest_sts))
        else:
            result_file.write("Nearest STs: {}\n".format(nearest_sts))

    # Write tsv table header
    table_header = ["Locus", "Identity", "Coverage", "Alignment Length", "Allele Length", "Gaps", "Allele"]
    table_file.write("\t".join(table_header) + "\n")
    rows = []
    for locus, allele_info in allele_matches.items():

        identity = str(allele_info["identity"])
        coverage = str(allele_info["coverage"])
        allele = allele_info["allele"]
        allele_name = allele_info["allele_name"]
        align_len = str(allele_info["align_len"])
        sbj_len = str(allele_info["sbj_len"])
        gaps = str(allele_info["gaps"])

        # Write alleles names with indications of imperfect hits
        if allele_name != "No hit found":
            allele_name_w_mark = locus+"_"+allele
        else:
            allele_name_w_mark = allele_name          
        
        # Write allele results to tsv table
        row = [locus, identity, coverage, align_len, sbj_len, gaps, allele_name_w_mark]
        rows.append(row)
        table_file.write("\t".join(row) + "\n")

        if allele_name == "No hit found":
            continue

        # Write query fasta output
        query_seq = allele_info["query_seq"]
        sbjct_seq = allele_info["sbjct_seq"] 
        homol_seq = allele_info["homol_seq"]
        
        if allele_info["match_priority"] == 1:
            match = "PERFECT MATCH"
        else:
            match = "WARNING"
        header = ">{}:{} ID:{}% COV:{}% Best_match:{}\n".format(locus, match, allele_info["identity"], 
                                                  allele_info["coverage"], allele_info["allele_name"])
        query_file.write(header)
        for i in range(0,len(query_seq),60):
            #print(query_seq[i:i+70])
            query_file.write(query_seq[i:i+60] + "\n")

        # Write template fasta output
        header = ">{}\n".format(allele_info["allele_name"])
        sbjct_file.write(header)
        for i in range(0,len(sbjct_seq),60):
            #print(sbjct_seq[i:i+70])
            sbjct_file.write(sbjct_seq[i:i+60] + "\n")

    # Write Allele profile results table 
    result_file.write(text_table(table_header, rows))

    # Write any notes
    if note != "":
       result_file.write("\nNotes: {}\n\n".format(note))

    # Write allignment output
    result_file.write("\n\nExtended Output:\n\n")
    make_aln(species, result_file, allele_matches, query_aligns, homol_aligns, sbjct_aligns)

    # Close all files
    query_file.close()
    sbjct_file.close()
    table_file.close()
    result_file.close()
