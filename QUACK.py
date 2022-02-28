#!/usr/bin/env python3

import sys

print("""				---------- WELCOME TO QUACK (Beta_version): ----\n
					Q -Quantification \n
					U - Utility \n
					A - And \n
					C - Contamination \n
					K - Killer \n
				-------------------------------------------------\n""")


if len(sys.argv) != 8:
    sys.exit("""ERROR! CHECK YOUR INPUT PARAMETERS!
This scripts takes a count_table generated with Symbiosis Evolution Group pipeline (zotu_table_expanded.txt) as well as a list of blank (negative control) library names, list of spikeins and conducts a quantification of obtained data.
It produces four files:
I) Table_with_classes.txt - where every zOTU is assigned to Symbiont, Other, PCR or Extraction Contaminant and PCR or Extraction Spikein class
II) Statistics_table.txt - with statistics about every library composition in terms of e.g contamination or spikein percentage 
III) Decontaminated_zOTU_table.txt - where all contaminants and spikeins are deleted, as well as libraries that sum of those were higher than ThresholdC
IV) Decontaminated_OTU_table.txt - table based on Decontaminated zOTU table and otus.tax file
Usage: QUACK.py <count_table> <list_of_blanks> <list_of_spikeins> <otus.tax> <ThresholdA; recommended value 10> <ThresholdB; recommended value 0.001> <ThresholdC; recommended value 30> 
Parameters:
<count_table>       Tab-delimited table produced by SEG pipeline. It must include experimental as well as blank libraries and taxonomy affiliation of each zOTU.
<list_of_blanks>    Text file with names of blank (negative control) libraries with description (PCR/Extraction_blank)
<list_of_spikeins> Text file with names of used spikeins with description (PCR/Extraction_spikein)
<otus_tax> One of the files produced by LSD script, used for taxonomy assignment in Decontaminated_OTU_table
<ThresholdA>        Value representing the thresholdA parameter; a unique genotype will be assigned as a contaminant UNLESS the maximum relative abundance it attains in at least one experimental library is more than ThresholdA * of the maximum relative abundance it attains in any blank library. Recommended value: 10
<ThresholdB>        Value representing the thresholdB parameter;  a unique genotype assigned previously as a symbiont will be assigned as "Other" UNLESS the maximum relative abundance it attains in at least one experimental library is more than ThresholdB. Recommended value: 0.001
<ThresholdC>        Value representing the thresholdC parameter; a library will be deleted UNLESS the % of contamination will be lower than ThresholdC. Recommended value: 30""")
   

Script, Input_count_table, List_of_blank_names, List_of_spikeins, otus_tax, ThresholdA, ThresholdB, ThresholdC  = sys.argv

##### Setting names of output files
Output_table = "Table_with_classes.txt"
Output_statistics = "Statistics_table.txt"
Decont_zotu_table = "Decontaminated_zOTU_table.txt"


print("Opening OTU table..................... ", end="")

COUNTS = open(Input_count_table, "r")

### Reading OTU_table and saving it as a list of lists
TABLE = []
for line in COUNTS:
    LINE = line.strip().split()
    TABLE.append(LINE)    
COUNTS.close()
print("OK!")   
### Reading Blank list and assigning it to one of the categories
print("Opening List of blanks..................... ", end="")
BLANKS = open(List_of_blank_names, "r")

Blank_PCR = []
Blank_Extr = []

print("OK!")  
for line in BLANKS:
    LINE_BLANK = line.strip().split()
    if len(LINE_BLANK) == 2 and LINE_BLANK[1] == "blank_PCR":
        Blank_PCR.append(LINE_BLANK[0])
    elif len(LINE_BLANK) == 2 and LINE_BLANK[1] == "blank_extr":
        Blank_Extr.append(LINE_BLANK[0])

        

if len(Blank_PCR) >= 1:
        print('Blank list proceeded succesfully,\nI am going to use %s of PCR blanks:' %len(Blank_PCR))
        for item in Blank_PCR:
            print(item, end=", ")
elif len(Blank_PCR) == 0:
    Blank_PCR.append("There is no PCR blanks")
    print("There is no PCR blanks")
print()


if len(Blank_Extr) >= 1:
        print('\nI am going to use %s of Extraction blanks:' %len(Blank_Extr))
        for item in Blank_Extr:
            print(item, end=", ")
elif len(Blank_Extr) == 0:
    Blank_PCR.append("There is no Extraction blanks")
    print("There is no Extraction blanks")
print()


###Reading SPIKEIN list and assigning it to one of the categories
print("Opening List of Spikeins..................... ", end="")
SPIKE = open(List_of_spikeins, 'r')

PCR_SPIKEIN = []
EXT_SPIKEIN = []


for line in SPIKE:
    LINE_SPIKE = line.strip().split()
    if len(LINE_SPIKE) == 2 and LINE_SPIKE[1] == "PCR_Spikein":
        PCR_SPIKEIN.append(LINE_SPIKE[0])

    elif len(LINE_SPIKE) == 2 and LINE_SPIKE[1] == "Extr_Spikein":
        EXT_SPIKEIN.append(LINE_SPIKE[0])
    else:
        print('This line could not be read: %s' % line)
        
print("OK!")      
print('\nSpikein list proceeded succesfully!')
### In case there is no spikeins:
if len(PCR_SPIKEIN) >= 1:
    for item in PCR_SPIKEIN:
        print("I am going to use following as PCR spikein: ", item, end = "")
elif len(PCR_SPIKEIN) == 0:
    PCR_SPIKEIN.append("There is no PCR spikein")
    print("There is no PCR spikein")
print()


if len(EXT_SPIKEIN) >= 1:
    for item in EXT_SPIKEIN:
        print("I am going to use following as Extraction spikein: ",item, end = "")
elif len(EXT_SPIKEIN) == 0:
    EXT_SPIKEIN.append("There is no Extraction spikein")
    print("There is no Extraction spikein")
print()


ROW_NO = len(TABLE) 
COL_NO = len(TABLE[0])

### Cleaning TABLE with deleting sequences belonging to Eukaryota, Chloroplast, Mitochondria or Archea #Add chimeras
print("Searching for Non-bacteria taxa..................... ", end="")
for row_no in range(ROW_NO - 1, 0, -1): # Iterates through the list backwards!
    if "Chimera" in TABLE[row_no][1]:
        TABLE[row_no][2] = "Non-Bacteria"
    if "Eukaryota" in TABLE[row_no][2]:
        TABLE[row_no][2] = "Non-Bacteria"
    if "Chloroplast" in TABLE[row_no][2]:
        TABLE[row_no][2] = "Non-Bacteria"
    if "Mitochondria" in TABLE[row_no][2]:
        TABLE[row_no][2] = "Non-Bacteria"
    if "Archaea" in TABLE[row_no][2]:
        TABLE[row_no][2] = "Non-Bacteria"
print("\nChimeras, Eukaryota, Chloroplast, Mitochondria and Archea reads has been recognized!")

ROW_NO = len(TABLE) ### re-creating variable ROW_NO again as number of rows has changed.

###Calculating the sum of reads without spikeins for decontamination
TOTALS = ["", "", "", "",""] ### filling first 5 elements in column
for col_no in range(5,COL_NO):
    Total = 0
    for row_no in range(1,ROW_NO):
        ###print(row_no, col_no, TABLE[row_no][col_no])### troubleshooting
        if not "Spikein" in TABLE[row_no][2]:
            if not "Non-Bacteria" in TABLE[row_no][2]:
                Total += int(TABLE[row_no][col_no])
    TOTALS.append(Total)
TABLE.append(TOTALS)

###### Finding the maximum relative abundance of each unique sequence (ASV, zOTU)
###### and assigning it to one of the classes.
STATS = ["Decision"]
for row_no in range(1, len(TABLE[:-1])):
    max_r_extr = 0
    max_r_PCR = 0
    max_r_real = 0
    Decision = 0
    for col_no in range(5, COL_NO):
        if TABLE[-1][col_no] == 0:
            relabund = 0
        else:
            relabund = int(TABLE[row_no][col_no])/TABLE[-1][col_no]
        ###print(row_no, col_no, TABLE[row_no][col_no])###
            if TABLE[0][col_no] in Blank_PCR:
                if max_r_PCR < relabund:
                    max_r_PCR = relabund
            elif TABLE[0][col_no] in Blank_Extr:
                if max_r_extr < relabund:
                    max_r_extr = relabund
            else:
                if max_r_real < relabund:
                    max_r_real = relabund

    if PCR_SPIKEIN[0] in TABLE[row_no][2]:
        Decision = "PCR_Spikein"
    elif EXT_SPIKEIN[0] in TABLE[row_no][2]:
        Decision = "Extraction_Spikein"
    elif "Brachybacterium" in TABLE[row_no][2]:
        Decision = "PCR_Contaminant"
    elif TABLE[row_no][2] == "Non-Bacteria":
        Decision = "Non-Bacteria"
    elif max_r_real > (max_r_extr + max_r_PCR) * float(ThresholdA):
        Decision = "Symbiont"
    elif max_r_extr > max_r_PCR * float(ThresholdA): 
        Decision = "Extraction_Contaminant"

    else:
        Decision = "PCR_Contaminant"
### Adding class "Other" for low abundant symbionts
    if Decision == "Symbiont":
        if max_r_real < float(ThresholdB):
            Decision = "Other"
        
    STATS.append(Decision)
TABLE[0].append("Class")
###Joining TABLE with STATS:
for row_no, item in zip(TABLE[1:], STATS[1:]): ### Iterates two lists and assignes item from one to another
    row_no.append(item)

###Deleting sum of reads without spikeins as we don't need this anymore
del(TABLE[-1])

print("Saving Table with classes....................... ", end = "")
###Printing the output_table
OUTPUT = open(Output_table, "w")
for LINE in TABLE:
    for item in LINE[:-1]:
        print(item, end="\t", file = OUTPUT)
    for item in LINE [-1:]:
        print(item, end='\n',file = OUTPUT)
print("OK!")  
   
ROW_NO = len(TABLE) 
COL_NO = len(TABLE[0])


### Computing statistics for our libraries:
Statistics = []
HEADINGS = TABLE[0]
HEADINGS.insert(5, "Statistics") ### Adding column "Statistics"
Statistics.append(HEADINGS[5:-1]) ###Only headings with names of libraries

###Calculating sum of all counts in a library 
TOTALS = ["Total_reads_per_library"]
for col_no in range(5,COL_NO-1):
    Total = 0
    for row_no in range(1,ROW_NO):
        if not "Non-Bacteria" in TABLE[row_no][-1]: ###Excluding non-Bacteria
            Total += int(TABLE[row_no][col_no])
    TOTALS.append(Total)
###Calculating a sum of each class reads in each library
OTHERS = ["Others_reads"]
SYM = ["Symbionts_reads"]
PCR_CON = ["PCR_Contaminants_reads"]
EXT_CON = ["Extraction_Contaminants_reads"]
PCR_SPIKE = ["PCR_Spikeins_reads"]
EXT_SPIKE = ["Extraction_Spikeins_reads"]
SUM_OT_SY = ["Symbionts and Other_reads"]
for col_no in range(5,COL_NO-1):
    others = 0
    symbiont = 0
    pcr_cont = 0
    extr_cont = 0
    pcr_spikein = 0
    ext_spikein = 0
    sum_other_sym = 0
    for row_no in range(1,len(TABLE)):
        if "Other" in TABLE[row_no][-1]:
            others += int(TABLE[row_no][col_no])
        if "Symbiont" in TABLE[row_no][-1]:
            symbiont += int(TABLE[row_no][col_no])
        if "PCR_Contaminant" in TABLE[row_no][-1]:
            pcr_cont += int(TABLE[row_no][col_no])
        if "Extraction_Contaminant" in TABLE[row_no][-1]:
            extr_cont += int(TABLE[row_no][col_no])
        if "PCR_Spikein" in TABLE[row_no][-1]:
            pcr_spikein += int(TABLE[row_no][col_no])
        if "Extraction_Spikein" in TABLE[row_no][-1]:
            ext_spikein += int(TABLE[row_no][col_no])
    OTHERS.append(others)
    SYM.append(symbiont)
    PCR_CON.append(pcr_cont)
    EXT_CON.append(extr_cont)
    PCR_SPIKE.append(pcr_spikein)
    EXT_SPIKE.append(ext_spikein)
    
Statistics.append(SYM)
Statistics.append(OTHERS)
Statistics.append(PCR_CON)
Statistics.append(EXT_CON)
Statistics.append(PCR_SPIKE)
Statistics.append(EXT_SPIKE)

### Summing all symbionts reads (others + symbionts)
for (item1, item2) in zip(OTHERS[1:], SYM[1:]):
    SUM_OT_SY.append(item1+item2)
Statistics.insert(3,SUM_OT_SY)
Statistics.append(TOTALS)


###Deleting unique genotypes (ASVs, zOTUS) belonging to "Other" class
###for row_no in range(ROW_NO - 1, 0, -1): # Iterates through the list backwards!
   ### if "Other" in TABLE[row_no][-1]:
      ###  del(TABLE[row_no])


### Our Statistics table looks like this:
### 1st row - Libraries headings
### 2nd row - sum of Symbionts reads
### 3rd row - sum of Others reads
### 4th row - sum of Symbionts and Others reads
### 5th row - sum of PCR contaminants reads
### 6th row - sum of Extraction contaminants reads
### 7th row - sum of PCR spikein reads
### 8ht row - sum of Extracion spikein reads
### 9th row - sum of reads in the library

        
SYMB = ["%_of_all_Symbionts_reads_in_lib"]
PCR_C = ["%_of_PCR_contamination_reads_in_lib"]
EXT_C = ["%_of_Extraction_contamination_reads_in_lib"]
CONT = ["%_of_all_contaminants_reads_in_lib"]
PCR_SPIKE = ["%_of_PCR_Spikein_reads_in_lib"]
EXT_SPIKE = ["%_of_Extraction_Spikein_reads_in_lib"]
SPIKE = ["%_of_all_Spikeins_in_lib"]
CONT_SPIKE = ["%_of_sum_of_Spikeins_and_Contaminants"]
SYM_EXT_SPIKE_RAT = ["Symbiont_to_Extraction_Spikein_Ratio"]
###In the future add Symbiont to Ext_Spike in ratio * copies of spikein added
for col_no in range(1, len(Statistics[0])):
    Symbionts = 0
    Con_PCR = 0
    Con_Ext = 0
    PCR_Spike = 0
    Extr_Spike = 0
    sym_ext_rat = 0

    for row_no in range(1, len(Statistics[:-1])):
        if Statistics[-1][col_no]!= 0:
            Symbionts = int(Statistics[3][col_no])/int(Statistics[-1][col_no]) * 100
            Con_PCR = int(Statistics[4][col_no])/int(Statistics[-1][col_no]) * 100
            Con_Ext = int(Statistics[5][col_no])/int(Statistics[-1][col_no]) * 100
            PCR_Spike = int(Statistics[6][col_no])/int(Statistics[-1][col_no]) * 100
            Extr_Spike = int(Statistics[7][col_no])/int(Statistics[-1][col_no]) * 100
        else:
            Symbionts = 0
            Con_PCR = 0
            Con_Ext = 0
            PCR_Spike = 0
            Extr_Spike = 0
            
        if Statistics[7][col_no] == 0:
            sym_ext_rat = 0
        else:
            sym_ext_rat = int(Statistics[3][col_no])/int(Statistics[7][col_no])
    SYMB.append(Symbionts)    
    PCR_C.append(Con_PCR)
    EXT_C.append(Con_Ext)
    PCR_SPIKE.append(PCR_Spike)
    EXT_SPIKE.append(Extr_Spike)
    SYM_EXT_SPIKE_RAT.append(sym_ext_rat)
for (item1, item2) in zip(PCR_C[1:], EXT_C[1:]): ###Summing % of contaminants
    CONT.append(round(item1+item2, 4))

for (item1, item2) in zip(PCR_SPIKE[1:], EXT_SPIKE[1:]): ###Summing % of Spikeins
    SPIKE.append(round(item1+item2, 4))
    
for (item1, item2) in zip(CONT[1:], SPIKE[1:]): ###Summing % of Spikeins and Contaminants
    CONT_SPIKE.append(round(item1+item2, 4))
    
EMPT = ["#########"]    
Statistics.append(EMPT)
Statistics.append(SYMB)
Statistics.append(PCR_C)
Statistics.append(EXT_C)
Statistics.append(CONT)
Statistics.append(PCR_SPIKE)
Statistics.append(EXT_SPIKE)
Statistics.append(SPIKE)
Statistics.append(CONT_SPIKE)
Statistics.append(SYM_EXT_SPIKE_RAT)

print("Saving Table with statistics....................... ", end = "")
STAT = open(Output_statistics, "w")
for LINE in Statistics:
    for item in LINE[:-1]:
        print(item, end="\t", file = STAT)
    for item in LINE [-1:]:
        print(item, end='\n',file = STAT) 
print("OK!")


   
###Creating the decontaminated Table
# Deleting everything that is not Symbiont:
TABLE[0].remove("Statistics") ###Deleting "Statistics"
for row_no in range(ROW_NO - 1, 0, -1):
    if not ("Other") in TABLE[row_no][-1]:
        if not ("Symbiont") in TABLE[row_no][-1]:
            del(TABLE[row_no])
            
DEC = open(Decont_zotu_table, "w")


###Adding information about % of total contamination in the library:
G = ["", "", "", "",] + Statistics[-6] + [""] ### Filling first 4 cells and last one in the last row
TABLE.append(G)



col_nos_to_keep = [0,1,2,3,4]
for col_no in range(5,(COL_NO-1)):
    if float(TABLE[-1][col_no]) < float(ThresholdC):
        col_nos_to_keep.append(col_no)
        
print("Saving Decontaminated zOTU Table....................... ", end = "") 
DEC = open(Decont_zotu_table, "w")
for row_no in range(len(TABLE)-1):
    for col_no in col_nos_to_keep:
        print(TABLE[row_no][col_no], end='\t', file = DEC)
    print(TABLE[row_no][-1], end='\n', file = DEC)
DEC.close()
print("OK!")
   

OTU = open("Decontaminated_zOTU_table.txt", "r")
OTU_TABLE = []
for line in OTU:
    LINE = line.strip().split()
    if line.startswith("OTU_ID"):
        COUNT_headings = line.strip().split()[4:-1]    ### Keeping the names of libraries
    else:
        OTU_TABLE.append(LINE)   
OTU.close()

otu_dict = {}
for row_no in range(0, len(OTU_TABLE)):
    otu_key = OTU_TABLE[row_no][1]
    if not otu_key in otu_dict.keys():
         otu_dict[otu_key] = OTU_TABLE[row_no][4:-1]
    else:
        otu_dict[otu_key] = [sum(map(int, i)) for i in list(zip(otu_dict[otu_key], OTU_TABLE[row_no][4:-1]))]

##### Adding taxonomy info to DICT
TAX = open(otus_tax, "r")
OTU_TAX = []
for line in TAX:
    LINE = line.strip().split()
    OTU_TAX.append(LINE)
    
for list in OTU_TAX:
    list[0] = list[0].lower()
        
        

for row_no in range(0, len(OTU_TAX)):   
    if OTU_TAX[row_no][0] in otu_dict.keys():
        if len(OTU_TAX[row_no]) > 1:
            otu_dict[OTU_TAX[row_no][0]].append(OTU_TAX[row_no][1])
        else:
            otu_dict[OTU_TAX[row_no][0]].append("unassigned")
TAX.close()


            
for row_no in range(0, len(OTU_TABLE)):
    if OTU_TABLE[row_no][1] in otu_dict.keys():
        if not "Symbiont" in otu_dict[OTU_TABLE[row_no][1]]:
            if not "Other" in otu_dict[OTU_TABLE[row_no][1]]:
                otu_dict[OTU_TABLE[row_no][1]].append(OTU_TABLE[row_no][-1])
                
print("Adding sequences...................... ", end="")
###We are adding 1 to the end of our dictionary to 
for row_no in range(0, len(OTU_TABLE)):
    if otu_dict[OTU_TABLE[row_no][1]][-1] != 1 and OTU_TABLE[row_no][1] in otu_dict.keys():
        otu_dict[OTU_TABLE[row_no][1]].append(OTU_TABLE[row_no][3])
        otu_dict[OTU_TABLE[row_no][1]].append(1)
print("OK!")        

                
COUNT_headings.insert(0,"#OTU")
COUNT_headings.insert(1,"Taxonomy")
COUNT_headings.insert(2,"Class")
COUNT_headings.insert(3,"Sequence")
data = []
data.append(COUNT_headings)      
for otu in otu_dict.keys():
    data.append([otu] + [otu_dict[otu][-4]] + [otu_dict[otu][-3]] + [otu_dict[otu][-2]] + otu_dict[otu][:-4])
    
print("Saving Decontaminated OTU Table....................... ", end = "")   


with open("Decontaminated_OTU_Table.txt", "w") as bigFile:
    for LINE in data:
        for item in LINE[:-1]:
            print(item, end="\t", file = bigFile)
        for item in LINE [-1:]:
            print(item, end='\n',file = bigFile)
bigFile.close()
print("OK!")


print("""
DONE! May QUACK bring you luck!
         ,-.
       ,--' ~.).
     ,'         `.
    ; (((__   __)))
    ;  ( (#) ( (#)
    |   \_/___\_/|
   ,"  ,-'    `__".
  (   ( ._   ____`.)--._        _
   `._ `-.`-' \(`-'  _  `-. _,-' `-/`.
    ,')   `.`._))  ,' `.   `.  ,','  ;
  .'  .     `--'  /     ).   `.      ;
 ;     `-        /     '  )         ;
 \                       ')       ,'
  \                     ,'       ;
   \               `~~~'       ,'
    `.                      _,'
      `.                ,--'
        `-._________,--') """)
