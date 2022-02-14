#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jun 30 10:24:11 2020

@author: susannasabin

Usage:
    python extract_pairwise_stats.py [PATH/TO/OUTDIR]
    
Summary:
    Extract information from the outfile.txt output from LDHat's 'pairwise' program and save to a tsv file with information on the model being tested.
    See: http://ldhat.sourceforge.net/manual.pdf for information on 'pairwise'
    
    Input:
        PATH/TO/OUTDIR is the only necessary input. This script assumes the following file structure: 
            - the specified output directory contains 
                a) subdirectories named for each simulation model
                    - each model directory contains the output of ldHat pairwise program, specifically prefixed outfiles
                b) a table containing model names, mu, and recombination rate (C) for each simulated model
        The output of this script will be deposited in the same ldhat output parent directory
                

"""
import glob, sys, re

# Import the table with model information, taking the model information table from the output directory, which is the provided input argument.
outdir = (sys.argv[1])
fullmodlist = (outdir + "/modlist.txt")
fullmodfile = open(fullmodlist)
allmodlist = fullmodfile.readlines()
modlist = (allmodlist[1:])

# Define default standard output to an object.
original = sys.stdout

# Create an empty list to hold all the information that you will ultimately put into a tsv file.
masterOb = []

# Loop through model types, parsing the mu and C values individually as you go.
for line in modlist:
    divline = (line.split())
    mod = (divline[0])
    mu = float(divline[1])
    oldC = float(divline[2])
    # Make a list of all output files from pairwise.
    outreps = glob.glob(outdir + "/" + mod + "/replicate_*_outfile.txt")
    # Loop through the ouput files.
    for outfile in outreps:
        out = open(outfile)
        outlines = out.readlines()
        #Extract the floating point theta from the LDHat output file.
        theta_ln = outlines[2]
        theta_out = re.findall(r"[-+]?\d*\.\d+|\d+", theta_ln)
        theta = float(theta_out[0])
        #Extract the line containing maximum rho (4Ner) over the region and its likelihood.
        rholk_ln = outlines[4]
        rholk = re.findall(r"[-+]?\d*\.\d+|\d+", rholk_ln)
        #Extract rho!
        rho = float(rholk[1])
        #Extract likelihood!
        lk = float(rholk[2])
        r = ( rho / theta )
        newC = ( r * mu )
        # Turn all the numeric values into strings for inclusion in the table.
        oldC = str(oldC)
        newC = str(newC)
        mu = str(mu)
        r = str(r)
        theta = str(theta)
        rho = str(rho)
        lk = str(lk)
        # Write all the important stuff to a vector
        vtr = [mod, oldC, newC, mu, r, theta, rho, lk]
        # Create an object containing tabs to properly format the vector you just created for ultimate inclusion in a tsv.
        col_sep = str("\t")
        tvtr = col_sep.join(vtr)
        # Append vector to master object.
        masterOb.append(tvtr)



#Write a header for the output file.
head = ["mod", "simRecomb", "ldHatRecomb", "simMu", "ldHat_r", "ldHat_theta", "ldHat_rho", "ldHat_lk"]
# Turn the header into a tab separated line
thead = col_sep.join(head)
# Combine header with information in masterOb
masterTab = [head, masterOb]
#Write all lines of masterTab to a new text file.
with open(outdir + "/pairwise_data_allMods.txt", "w") as file:
    # Temporarily redirect output stream to the new file.
    sys.stdout = file
    for line in masterTab:
        print(line)
    sys.stdout = original


    

        
        
        
        



