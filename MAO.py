#! /usr/bin/env python3

import sys
import os

if len(sys.argv) != 2:
    sys.exit('----- project_3.py v. 1.0, M. Buczek, 12th Feb 2021 -----\n'
             'This script produce four files, from COI zOTU table:\n' 
             '- barcode.txt that contains info about most abundand COI barcode, taxonomy and bacteria presence\n' 
             '- barcode.fasta, containing most abundand Eucaryotic COI barcode per each library/sample\n'
             '- euc5.txt, containing information about Eucaryotic COI sequences that represents at least 5% of total Eucaryotic reads per sample\n'
             '- bac.txt, containing information about bacterial COI sequences\n'
             '     e.g.,  project_3.py OTUs.txt all_zotu_table_expanded.txt\n')   
             
### Assigning script arguments to variables
Script, path_to_input_file = sys.argv

INPUT = open(path_to_input_file, "r")
OUTPUT_BAC = open("bac_" + os.path.splitext(path_to_input_file)[0] + ".txt", mode= "w") # bac.txt
OUTPUT_EUC = open("euc5_" + os.path.splitext(path_to_input_file)[0] + ".txt", mode= "w") # euc5.txt
OUTPUT = open("barcode_" + os.path.splitext(path_to_input_file)[0] + ".txt", mode= "w") # barcode.txt
OUTPUT2 = open("barcode_" + os.path.splitext(path_to_input_file)[0] + ".fasta", mode= "w") # barcode.fasta

# read input file as list of lists
TABLE = []
for line in INPUT:
    TABLE.append(line.strip().split())    
INPUT.close()

ROW_NO = len(TABLE) 
COL_NO = len(TABLE[0])

### Divide table to Bacteria (BAC) and Eucaryote (EUC)
HEAD = TABLE[0]

EUC = [HEAD]
BAC = [HEAD]
BAC_YES = []
NO_BAC_index = []
for row_no in range(1,ROW_NO):
    if TABLE[row_no][2].startswith("Bacteria"): # if taxonomy starts with "Bacteria" save whole row in output_bac
        BAC.append(TABLE[row_no])
        for col_no in range(5, COL_NO): # this applies only to "Bacteria" rows
            if int(TABLE[row_no][col_no]) > 0: # if number of reads is larger than 0, add sample name to the list BAC_YES
                BAC_YES.append(HEAD[col_no])

    else:
        EUC.append(TABLE[row_no]) # if taxonomy do not starts with "Bacteria" add the row to EUC table as a list

    BAC_YES = list(set(BAC_YES)) # list of samples that have at least one bacterial zotu, takes unique values from BAC_YES

for col_no in range(5, COL_NO):
    count = 0 # number of zOTUs with bacterial reads per sample
    for row_no in range(1, len(BAC)):
        if 0 < int(BAC[row_no][col_no]): # I do not know why it has to be int() here
            count += 1
    if count == 0:
        NO_BAC_index.append(col_no)

BAC_2 = [x[:] for x in BAC] # deep copy, believe me it was needed
for row_no in range(0,len(BAC)):
    for item in sorted(NO_BAC_index, reverse = True):  
        del BAC_2[row_no][item] # removes samples that do not contain bacteria
    for item in BAC_2[row_no][:-1]:
        print(item, file=OUTPUT_BAC, end="\t")
    print(BAC_2[row_no][-1], file=OUTPUT_BAC)

### Adding the TOTAL row to Eucaryotic table --- added-up counts of all reads for the same sample
TOTALS = TABLE[0][:4]
ROW_NO = len(EUC)
for col_no in range(4,COL_NO):
    Total = 0
    for row_no in range(1,ROW_NO):
        Total += int(EUC[row_no][col_no])
    TOTALS.append(Total)
EUC.append(TOTALS)
# print(EUC[-1])

### Translating counts into relative abundances, with keeping both tables
EUC_ABUND = [x[:] for x in EUC] # deep copy
for col_no in range(4, COL_NO):
    for row_no in range(1, ROW_NO):    # Note: does NOT include the TOTAL row!
        if EUC[-1][col_no] > 0:
            EUC_ABUND[row_no][col_no] = round(float(EUC[row_no][col_no])/EUC[-1][col_no], 4)
        else:
            EUC_ABUND[row_no][col_no] = 0

### Creating new table with most abundand barcode per library (sample) and fasta file
NEW_HEAD = ["Sample", "Sequence", "Abundance", "Abundance > 5%", "Bacteria", "Taxonomy"] # new col names
for item in NEW_HEAD[:-1]:
    print(item, file=OUTPUT, end="\t")
print(NEW_HEAD[-1], file=OUTPUT)

EUC5_YES = []
for col_no in range(5, COL_NO):
    first = 0 # item to store highest abundance value
    first_index = 0
    treshold = 0.05 # abundance treshold
    count = 0 # number of zOTUs with abundance above treshold value
    for row_no in range(1, ROW_NO):
        if first <= EUC_ABUND[row_no][col_no]:
            first = EUC_ABUND[row_no][col_no]
            first_index = row_no # stores row number of otu with highest abundance
        if treshold < EUC_ABUND[row_no][col_no]:
            count += 1
            EUC5_YES.append(row_no)

    if HEAD[col_no] in BAC_YES:
        bacteria = "yes"
    else:
        bacteria = "no"

    EUC5_YES = list(set(EUC5_YES)) # list of zotus that in at least one sample represents at least 5% of total Eucaryotic reads

    to_add = [HEAD[col_no],EUC[first_index][3], first, count, bacteria, EUC[first_index][2]]
    to_fasta = to_add[:2]
    print(">", to_fasta[0], "\n", to_fasta[1], file=OUTPUT2, sep="" )
    for item in to_add[:-1]:
        print(item, file=OUTPUT, end="\t")
    print(to_add[-1], file=OUTPUT)

for item in HEAD[:-1]:
    print(item, file=OUTPUT_EUC, end="\t")
print(HEAD[-1], file=OUTPUT_EUC)

for row_no in range(1,ROW_NO):
    if row_no in EUC5_YES:
        for item in EUC[row_no][:-1]:
            print(item, file=OUTPUT_EUC, end="\t")
        print(EUC[row_no][-1], file=OUTPUT_EUC)
