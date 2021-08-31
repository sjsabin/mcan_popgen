#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 26 2021

@author: susannasabin

This code calls a consensus sequence based on basic variant counts from custom SLiM 3
    fasta output. The input is specifically fasta files subsampled to 25 sequences per file.
    For this particular project, we are subsampling 25 sequences 17 times from the same 
    SLiM replicate.

Subsample directory refers to the set of subsamples taken from a single SLiM replicate. 
Do not mix subsamples from different replicates in the same directory.
    
Usage: python fasta_to_consensus.py SUBSAMPLE/DIR OUTPUT/DIR
"""

import pandas as pd
import sys, glob, collections, os, random

# Define the input and output directories based on arguments.
fadir = (sys.argv[1])
outdir = (sys.argv[2])


# Create a new dictionary to hold the consensus sequences for each subsample.
con = {}
# Make sure the dictionary can be appended to like a list.
con = collections.defaultdict(list)

# Get list of subsampled ms files in directory.
falist = glob.glob(fadir + "/*.sub*.fasta")

# STEP 1: Prepare parts of file name for later use.
filenm = os.path.basename(falist[0])
filenmsplit = filenm.split(".")


# STEP 2: Create a list of positions and save it for later!
#Extract the 'positions' line from the ms file
# Open the ms file
fa = open(falist[0])
# Read lines of the ms files to a variable
fa_lines = fa.readlines()
# Write header to a separate variable
head = fa_lines[0]
#STEP 3: Generate a consensus sequence!
# Loop through all the ms files in ms list.
for fafile in falist:    
    # Open the ms file
    fa = open(fafile)
    # Read lines of the ms files to a variable
    fa_lines = fa.readlines()
    # Remove header line from fa_lines variable
    fa_lines = fa_lines[1:]
    # Create a list to hold all sequence headers for later use.
    seqhead = []
    # Create a list to hold all sequence lines
    seqs = []
    # Create a dictionary to hold the unseparated sequences (i.e. one big string per sequence)
    seqdict = {}
    # Put all sequence lines into a single variable
    for line in fa_lines:
        if ">" not in line:
            seqs.append(line)
    print("Length of seqs list: ", len(seqs))
    # Transfer sequence headers into an indexed list for later use
    for line in fa_lines:
        if ">" in line:
            seqhead.append(line)
    # Transfer sequences from list to dictionary
    for i in range(len(seqs)):
        seqdict[i]=seqs[i]
    print("Length of seqdict: ", len(seqdict))
    # Create a new dictionary to hold the sequences separated allele by allele,
    #   and modify it so we can append items to elements in the dictionary as
    #   if they were lists.
    sepseq = {}
    sepseq = collections.defaultdict(list)
    # Separate the alleles in each sequence and place them in the new dictionary per sequence
    for i in range(len(seqdict)):
        tmp = seqdict[i]
        #print("Length of sequence ", i, ":", len(tmp))
        for c in range(len(tmp)):
            sepseq[i].append(tmp[c])
    # Turn the new dictionary into a data frame. Rows will be sites and columns
    #   will be sequences.
    dictdf = pd.DataFrame(sepseq)
    # Calculate the length of the data frame to establish the last row (new lines only).
    end = len(dictdf)
    # Calculate the range of lines that include useful information.
    cut = (end-1)
    # Creat a new dataframe with the last row cut off.
    dictdfcut = dictdf.iloc[0:cut,:]
    # Turn the data frame into a series of lists.
    dflist = dictdfcut.values.tolist()
    #dflist = dictdf.values.tolist()
    # Count the reference alleles. Based on there being 25 genomes in each subsample,
    #   we know that a majority allele will be represented by 13 or greater occurrences.
    for i in range(len(dflist)):
        countA = int(dflist[i].count('A'))
        countT = int(dflist[i].count('T'))
        countC = int(dflist[i].count('C'))
        countG = int(dflist[i].count('G'))
        counts = []
        counts = [countA, countT, countC, countG]
        nucs = []
        if max(counts) == countA:
            nucs.append('A')
        if max(counts) == countT:
            nucs.append('T')
        if max(counts) == countC:
            nucs.append('C')
        if max(counts) == countG:
            nucs.append('G')
        n = random.choice(nucs)
        con[i].append(n)
    # Report on progress...
    print(" ")
    print("One consensus sequence called for " + str(len(dictdf)) + " sites, based on " + str(len(seqs)) + " sequences...")
# Report on progress ...
print(str(len(falist)) + " consensus sequences called for model " + filenmsplit[0] + ", replicate " + filenmsplit[1] + ".")    

# Convert the consensus dictionary to a data frame.
condf = pd.DataFrame.from_dict(con)
# Create a list to house each printed line of sites.
conlist = []
# Create a separator string object to use for proper file formatting.
file_sep = str("\n")
# Loop through rows in the consensus sequence data frame to print sequences as strings.
for i in range(len(condf)):
    tmp = condf.iloc[i,:].to_string(header=False, index=False)
    tmp = tmp.replace('\n','')
    tmp = tmp.replace(' ','')
    conlist.append(tmp)


# Build the sites file header.
# Determine how many separate sequences are in the ms file by creating a new object excluding lines 0-2 (i.e. the "//" header, segsites, and positions lines) and calculating its length.
num_seqs = str(len(conlist))
#Create object with blank space separators as the interlocuter to create the sites header.
sites_seps = str(" ")
# Determine number of sites in each sequence
num_sites = str(len(conlist[0]))
#Join the strings representing the number of sequences, the number of sites, and the haplotype/genotype flag into a list.
sites_head_ls = [num_seqs, num_sites, "1"]
#Join the parts of the list by the blank space separator to create the header string for the sites file.
sites_head = sites_seps.join(sites_head_ls)

new_sites = str()

for num, line in enumerate(conlist, start=1):
    #Write line from seqs to a temporary object.
    tmp = line
    #Convert the index of the seqs list to a string instead of an integer. 
    num = str(num)
    #Write the sequence header.
    tmp_head = ">consensus_" + num
    #Combine the sequence header and sequence into a single object.
    tmp_dup = [tmp_head, tmp]
    #Write the header and sequence unit to a string object, joined by the space separator.
    add = file_sep.join(tmp_dup)
    #Add the combined string to the object that contains all other sequences and headers in the correct format.
    new_sites = new_sites + add + "\n"

original = sys.stdout
# Combine the sequence blocks with the sites file header.    
sites_ls = [sites_head, new_sites]
#Write all lines to a new text file.
with open(outdir + "/" + filenmsplit[0] + "." + filenmsplit[1] + ".consensus.fasta", "w") as cons_fasta:
    #Temporarily redirect the output stream to the new file.
    sys.stdout = cons_fasta
    for line in sites_ls:
        print(line)
    #Redirect the output stream back to the original.
    sys.stdout = original

sys.stderr.write("Done!")



    

    
