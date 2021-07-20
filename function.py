# import of package
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import rc
import pandas as pd
from datetime import datetime
from statistics import mean
from collections import Counter
import logomaker
from itertools import repeat, chain
import os  
import shutil
import re
import math 

from difflib import SequenceMatcher
import mappy as mp



# Classic functions 
#complement
# return a list containing the complement nucleotides of sequence s
#@param
#@s, a string corresponding to the sequence to reverse
def complement(s): 
    basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 
    letters = list(s) 
    letters = [basecomplement[base] for base in letters] 
    return ''.join(letters)

# return the reverse list of the list s
def revcom(s):
    return complement(s[::-1]) # s[::-1] means that the elements of s are sorted in  reverse order

# # Definition of variables

# Definition of the different treatments
concat = 0 # 0 is no concat or 1 if concatenation
raw = -1 # 1 : raw file (fastQ without cleaning); -1 : file after treatment (cleaning, concatenation, etc..)
treatment_type = 2 # 1 : removal of sequences that don't have GGG; 2 : removal of sequences that don't have [CG]{3}; -1 if any

if_controle = 0 # decide if we take into account a control batch 1 if yes o otherwise
#controle = C2 # controle proportions of nucléotides

UMI_pos = [0] # Position and size of the UMI:  UMI_pos=[x] if the UMI is the x first nucleotides 
# and UMI_pos=[x,y] if the UMI is the x first nucleotides and the y last
Nb_MU = 6 #number of MU plotted
MU_thresh = 26 # number of MU for stems only processed by the SequenceDiversity fonction
H_thresh = 16 # H_thresh : H_thresh is the number of combination allowed by the SequenceDiversity fonction
# ex : for Hx,y with x<=3 and y<=3, there are 4*4 combinations == 16 
MU_in_H_thresh = 6 # number of MU within hairpins between the two stems processed by the SequenceDiversity fonction

# definition of the sequences of interest
#substrate_name = 'A1' # file name, there has to be a folder named 'figures' and another named 'file' in 'figures' in the directory
#stem = "CGCGAATTAACGCGTTAACAT"
substrate_name = 'X4' 
stem = "CGCGAATTAACGCGACAACAT"
compl_stem = revcom(stem)
stem_RC = r"ATG%s" % (revcom(stem))
hairpin_stem1 = stem[-(len(stem)-2):]
hairpin_stem2 = stem[-(len(stem)-4):]
tag = "CAT"
MU = 2 # size of the standard mobile unit in this batch, tag excluded
compl_tag = revcom(tag)

# Other variables
polyA = (r"A{6,18}")
T_patch = (r"T{6,18}")

regex_NN = '[ATGC]*?CAT' # regex used to find all the mobile units
regex_NNCA = '[ATGC]*?CA' # regex used to find if there is a last mobile unit which tag has been cut 

regex_MU_cuted = 'AACA'
regex_MU = 'AACAT'
regex_repet_stem = r'^( + regex_MU + )*( + regex_MU_cuted + )?%s'%polyA
print(regex_repet_stem )



# position of the diverse files
path_to_file = r"C:\Users\Camille\Documents\GitHub\RNASeq_analysis\SequencingData\RNASeq_1/"
position_mothur = r"C:\Users\Camille\Documents\GitHub\RNASeq_analysis\Mothur\fastqtest" # position of the directory inside which mothur process the files (see the .batch file)
excel_pos = path_to_file + "%s.xlsx"%(substrate_name) # the name of the Excel file in which the datas are stored
path_to_figure = r"C:\Users\Camille\Documents\GitHub\RNASeq_analysis\Figures"


# # Processing Cell
# file_to_process = Preprocessing(path_to_file, substrate_name, concat, raw)
# print(file_to_process)

# seq, final_file_to_process = UnexpectedSeqRemover(path_to_file, file_to_process, treatment_type)
# S, NUS_Stem = SFishing(regex_repet_stem, seq, regex_AA, regex_AACA)
# H = HpFishing(seq, regex_AA, regex_AACA)
# X_S, Prop_S, X_Hxx, Prop_Hxx, X_H, Prop_H1, Prop_H2, X, Prop = Sequence_diversity(S, H, MU_thresh, H_thresh, MU_in_H_thresh)




# # FastQToFasta
# Function which convert a FastQ file into a Fasta file and a quality file.
    # @param
    # @path_to_fastq, path to the file to process, should finish by the .fastq extention
    # @read_identifier, first three letters of the read identifier from the raw fastq file.
def FastQToFasta(path_to_fastq, read_identifier = '@M0'):
     
    output_fasta = open("%s.fasta"%path_to_fastq.replace(".fastq",""), "w")
    output_qual = open("%s.qual"%path_to_fastq.replace(".fastq",""), "w")
    i = 0 
    for name, seq, qual in mp.fastx_read(path_to_fastq):

        i += 1 #count the number of read in the file
        newFastQID = name.replace(":","_").replace("@",">") 
        output_fasta.write(newFastQID + '\n' + seq + '\n')
        output_qual.write(newFastQID+ '\n' + qual + '\n')
    

    output_fasta.close()
    output_qual.close()


# # ReadCleaner
# Function which clean the fastq file, by cutting the nucleotides located after the polyA tail if any.
    # This function also count the number of PhiX within the file and the number of sequences trimmed during the treatment.
    # @param
    # @path_to_fastq, path to the file to process, should finish by the .fastq extention
    # @read, indicate if the file is a Read1 or Read2 file, can take 1 or 2 value
    # @read_identifier, first three letters of the read identifier from the raw fastq file.
def ReadCleaner(path_to_fastq, read_identifier, read):
    output = open(r"%s_cleaned.fastq" %(path_to_fastq.replace('.fastq',"")), "w")
    f = open(path_to_fastq, "r")

    i = 0 #count the number of read in the file
    PhiX = 0 #count the number of PhiX
    trimmed_read = 0 #count the number of trimmed reads
    
    for line in f:
        if line.startswith(read_identifier):
            i += 1 #count the number of read in the file
            ID = line.strip()
            seq = next(f).strip()
            plus = next(f).strip()
            quality = next(f).strip()
            
            if read == 1:
                if re.search(stem, seq) and re.search(polyA, seq):
                    end_poly = int(re.search(polyA, seq).end())
                    seq = seq[:end_poly]
                    quality = quality[:end_poly]
                    trimmed_read += 1 #count the number of trimmed reads
                elif seq == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                    PhiX += 1 #count the number of PhiX
                
                
            if read == 2: 
                if re.search(T_patch,seq)and re.search(stem_RC,seq):
                    end_stem = int(re.search(stem_RC,seq).end())
                    seq = seq[:end_stem]
                    quality = quality[:end_stem]
                    trimmed_read += 1 #count the number of trimmed reads
                elif seq == "NNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN":
                    PhiX += 1 #count the number of PhiX
                
        output.write(ID + '\n' + seq + '\n'+ plus + '\n'+ quality + '\n')
    print("There are %s reads in the raw file and %s : %s o/o PhiX" %(i, PhiX, str(PhiX/i*100)))
        
    f.close()
    output.close() 


# # Preprocessing
def Preprocessing(path_to_file, substrate_name, concat, raw):

    get_ipython().run_line_magic('cd', '"C:\\Users\\Camille\\Documents\\GitHub\\RNASeq_analysis\\SequencingData\\RNASeq_1"')

    #create the excel file which will compile the different treatments, with the first sheet regrouping general info
    df = pd.DataFrame([substrate_name, stem, MU, UMI_pos, Nb_MU, path_to_file],
                index=['fichier', 'stem', 'longueur_unite_mobile', 'position_UMI', 'Nb_Mu_plot', 'repertoire']) 


    if not os.path.isfile(excel_pos):
        with pd.ExcelWriter(excel_pos, engine = "openpyxl") as writer :  
            df.to_excel(writer, sheet_name = 'sequence_info')


    # Cleaning of the R1 and R2 files
    ReadCleaner(path_to_file + substrate_name + "_R1.fastq", '@M0', 1) # ouvre le fichier read1 et coupe après la queue polya
    ReadCleaner(path_to_file + substrate_name + "_R2.fastq", '@M0', 2) # idem pour le read2

    # In case of no concatenation, we don't need to use Mothur
    # We just need to convert the R1 FastQ file in Fasta file
    # And rename the fasta file to indicate that no concatenation was done


    if concat == 0 and raw == 1:
        FastQToFasta(path_to_file + substrate_name + '_R1.fastq', '@M0')
        if not os.path.isfile(r"%s\%s_R1_raw.fasta" % (path_to_file, substrate_name)):
            os.rename(r"%s\%s_R1.fasta" % (path_to_file, substrate_name) , r"%s\%s_R1_raw.fasta" % (path_to_file, substrate_name))
        else :
            os.remove(r"%s\%s_R1_raw.fasta" % (path_to_file, substrate_name))
            os.rename(r"%s\%s_R1.fasta" % (path_to_file, substrate_name) , r"%s\%s_R1_raw.fasta" % (path_to_file, substrate_name))
        file_to_process = substrate_name + '_R1_raw'
    
    elif concat == 0 and raw == -1:
        FastQToFasta(path_to_file + substrate_name + '_R1_cleaned.fastq', '@M0')
        if not os.path.isfile(r"%s\%s_without_concat.fasta" % (path_to_file, substrate_name)):
            os.rename(r"%s\%s_R1_cleaned.fasta" % (path_to_file, substrate_name) , r"%s\%s_without_concat.fasta" % (path_to_file, substrate_name))
        else :
            os.remove(r"%s\%s_without_concat.fasta" % (path_to_file, substrate_name))
            os.rename(r"%s\%s_R1_cleaned.fasta" % (path_to_file, substrate_name) , r"%s\%s_without_concat.fasta" % (path_to_file, substrate_name))   
        file_to_process = substrate_name + '_without_concat'
    
    # In case of concatenation, we need to move the cleaned R1 and R2 fastq file in the Mothur folder
    elif concat == 1 :  
        if raw == -1 :
            file_name1 = substrate_name + "_R1_cleaned"
            file_name2 = substrate_name + "_R2_cleaned"
            file_to_process = substrate_name + '_concat'
        if raw == 1 :
            file_name1 = substrate_name + "_R1"
            file_name2 = substrate_name + "_R2"
            file_to_process = substrate_name + '_raw_concat'
        shutil.move(r"%s\%s.fastq" % (path_to_file, file_name1),r"%s\%s.fastq" % (position_mothur, file_name1))
        shutil.move(r"%s\%s.fastq" % (path_to_file, file_name2),r"%s\%s.fastq" % (position_mothur, file_name2))

    
        #--------------------------
        # WARNING : path to change manually
    
        # Let's move in the Mothur folder
        get_ipython().run_line_magic('cd', '"C:\\Users\\Camille\\Documents\\GitHub\\RNASeq_analysis\\Mothur"')
        # Run Muther in batch mode
        get_ipython().system(' mothur RNA_seq.batch')
        # Let's move back to the working folder
        get_ipython().run_line_magic('cd', '"C:\\Users\\Camille\\Documents\\GitHub\\RNASeq_analysis\\SequencingData\\RNASeq_1"')
    
        #----------------------------

        # And move the Mothur output file called final.fasta into the working folder
        file_name3 = 'final'
        shutil.move(r"%s\%s.fasta" % (position_mothur, file_name3),r"%s\%s.fasta" % (path_to_file, file_name3))

    
        # We need to rename the Mothur output file
        #file_to_process = substrate_name + '_concat'
        if raw == -1 :
            if not os.path.isfile(r"%s\%s_concat.fasta" % (path_to_file, substrate_name)):
                os.rename(r"%s\%s.fasta" % (path_to_file, file_name3), r"%s\%s_concat.fasta" % (path_to_file, substrate_name))   
            else :
                os.remove(r"%s\%s_concat.fasta" % (path_to_file, substrate_name))
                os.rename(r"%s\%s.fasta" % (path_to_file, file_name3), r"%s\%s_concat.fasta" % (path_to_file, substrate_name))   
        elif raw == 1 :
            if not os.path.isfile(r"%s\%s_raw_concat.fasta" % (path_to_file, substrate_name)):
                os.rename(r"%s\%s.fasta" % (path_to_file, file_name3), r"%s\%s_raw_concat.fasta" % (path_to_file, substrate_name))   
            else :
                os.remove(r"%s\%s_raw_concat.fasta" % (path_to_file, substrate_name))
                os.rename(r"%s\%s.fasta" % (path_to_file, file_name3), r"%s\%s_raw_concat.fasta" % (path_to_file, substrate_name))   
    
        # And then empty the folder of Mothur 
    
        if raw == -1 :
            file_to_remove1 = substrate_name + "_R1_cleaned"
            file_to_remove2 = substrate_name + "_R2_cleaned"
        if raw == 1 :
            file_to_remove1 = substrate_name + "_R1"
            file_to_remove2 = substrate_name + "_R2"        

        shutil.move(r"%s\%s.fastq" % (position_mothur, file_to_remove1), r"%s\%s.fastq" % (path_to_file, file_to_remove1))
        shutil.move(r"%s\%s.fastq" % (position_mothur, file_to_remove2), r"%s\%s.fastq" % (path_to_file, file_to_remove2))
    
    return file_to_process


# # UnexpectedSeqRemover
def UnexpectedSeqRemover(path_to_file, file_to_process, treatment_type):
    # open the final_file_to_process
    f = open(r"%s\%s.fasta" % (path_to_file, file_to_process),"r")

    # create a file to store the sequences without stem
    # WARNING : if you run this program more than once, you will store the sequences twice
    sequence_removed_file = open('unexpected_no_stem_' + 'type_%s_' %(treatment_type) + file_to_process + '.fasta', "a")
    sequence_passing_test_file = open('sequence_with_stem_' + 'type_%s_' %(treatment_type) + file_to_process + '.fasta', "a")


    all_seq_list = [] # list of the odd line of f, which are the sequences (the even lines are the fasta identifier)
    seq = [] # list of the sequences that pass the test


    # Save in the all_seq_list all the sequences of the file f
    i = 0
    for line in f :
        if i %2 != 0 : # check if the line are even or odd. Select the odd lines here.
            all_seq_list.append(line.rstrip()) # Python count the /n as a character, it is then removed
        i += 1

    # For every element of all_seq_list, check if the element has the stem
    # If the element has the stem, it is added the the list seq
    # If not, the element is directly written in the sequence_removed_file

    unpass = 0
    pas = 0
    i = 0
    
    for x in all_seq_list :
        i += 1
        if treatment_type == 1:
            if x.find(stem) != -1 and x.find(tag) != -1 and x.find('GGG') != -1 :
                pas += 1
                sequence_passing_test_file.write(x + '\n')
                seq.append(x.rstrip()) 
            else:
                unpass += 1
                sequence_removed_file.write(x + '\n')

        elif treatment_type == 2:
            if re.search('[GC]{1,3}%s' %stem, x) or re.search('^%s' %stem, x) :
                pas += 1
                sequence_passing_test_file.write(x + '\n')
                seq.append(x.rstrip()) 
            else:
                unpass += 1
                sequence_removed_file.write(x + '\n')

        elif treatment_type == 3:  
            if x.find(stem) != -1 and x.find(tag) != -1 :
                pas += 1
                sequence_passing_test_file.write(x + '\n')
                seq.append(x.rstrip()) 
            else:
                unpass += 1
                sequence_removed_file.write(x + '\n')
            
    

    final_file_to_process = file_to_process + '_type_%s_' %(treatment_type)    
    print("The number of sequences passing the test %s is : %s : %s o/o in %s"%(treatment_type, str(pas/i*100), len(seq), file_to_process))
    print("There are %s seq : %s o/o not passing the test" %(unpass, str(unpass/i*100)))

    sequence_removed_file.close()
    sequence_passing_test_file.close()
    f.close()
    return seq, final_file_to_process


# # New Fishing functions
# Function which select sequences having only one stem
    # @param
    # seq, list of the sequences after removing the unexpected ones
    # regex, regex used to find all the mobile units NNCAT
    # regexca, regex used to find if there is a last mobile unit which tag has been cut NNCA
    # 
    # Returns
    # S : list of the following elements 
    # [list of all mobile units of the sequence  (ex : 'AA', 'TG'), number of mobile unit of the sequence, 
    # position of the sequence among all sequences substracted by the sequences without stem or without MU]
    
    # NUS_Stem : list of mobile unit with different length

def SFishing(select_seq, seq, regex, regexca):
    S = []
    NUS_Stem = []
    
    idx = 0
    nb_NNCA_MU = 0
    
    S_file = open('S_sequence_' + file_to_process + '_type_%s_' %(treatment_type) + '.fasta', "w")
    S_insert_file = open('S_sequence_with_insert_' + file_to_process + '_type_%s_' %(treatment_type) + '.fasta', "w")
    
    strict_cond = 0
    loose_cond = 0
    GAA_MU = 0
    broken_hp = 0
    
    for x in seq:
            
        # FISHING OF STEM SEQUENCES
        # -----------------------
        # Selection of sequences after the stem, if their is not a second stem afterwards (that would indicate an hairpin)  
        #if x.split(stem)[1].find(hairpin_stem2) == -1 or not re.search('%s[ATCG]*%s'%(stem, stem), x):

        
        if re.search(select_seq, x.split(stem)[1]):
            strict_cond += 1
            S_file.write(x + '\n')
            # Case where the UMI are only at the beginning of the read
            if len(UMI_pos) == 1:
                # Selection of the MU by looking for the regex, and excluding the UMI position and the last position
                S1 = re.findall(regex, x[UMI_pos[0]:-1])
            
            # Case where the UMI are located at the beginning and the end of the read
            # We need to exlude now the last position corresponding to the UMI when looking for mobile units
            elif len(UMI_pos) == 2: 
                S1 = re.findall(regex, x[UMI_pos[0]:len(x) - UMI_pos[1]])
                

            # If there is some MU :
            if len(S1) != 0: # on compte  les 'NNCA' si il y en a a la fin 
                # Selection of the last mobile unit in the sequence x
                last_MU = x.split(S1[-1])[-1]
                # Look if the last MU finishes by a NNCA instead of NNCAT by matching the regex regexca
                find_CA = re.findall(regexca, last_MU)
                # '#renvoi que NN' ?
                find_CA = [y for y in find_CA if len(y) == MU + 2] 

            # Storage of Mobile Unit (only the NN, not NNCAT) with correct length in S1
            S1 = [x[0:MU] for x in S1 if len(x) == MU + 3] 
            # Check if there are MU in the list storing the 'NNCA' mobile units 
            if len(find_CA) != 0:
                # Count the NNCA mobile units
                if find_CA[0] != '':
                    nb_NNCA_MU += 1
                    # Add the NNCA mobile unit to the S1 list
                    S1.append(find_CA[0])
                    
            # Store in the final S list [the list S1 of every standard mobile unit by sequence, 
            # the nb of mobile unit by sequence, position of the sequence among all sequences substracted by the sequences without stem or without MU]
            S.append([S1, len(S1), idx])
        
        elif re.search(stem, x) and not re.search('%s[ATCG]*%s'%(stem, hairpin_stem2), x):
            loose_cond += 1
            S_insert_file.write(x + '\n')
            if re.search(stem + '[ATCG]*(GAACA)', x):
                GAA_MU += 1 
            elif re.search('%s[CG]{3}'%(stem), x):
                broken_hp += 1
        idx += 1
        
    S_file.close()
    S_insert_file.close()
    
    print('The % of sequences finishing by AACAX for S is :' + str(nb_NNCA_MU/len(seq)))
    print('Among the %s sequences, \n %s : %s o/o of them have a strict stem composition and \n %s : %s o/o of them have a stem, MUs and insertions'%(idx, strict_cond, str(strict_cond/idx*100), (loose_cond), str((loose_cond)/idx*100)))
    print('GAACAT : %s and broken hp : %s'%(GAA_MU, broken_hp))
    return S, NUS_Stem 

#HpFishing
# Function which select the hairpins that are sequences having two stems
    # @param
    # seq, list of the sequences after removing the unexpected ones
    # regex, regex used to find all the mobile units
    # regexca, regex used to find if there is a last mobile unit which tag has been cut
    
    # Returns
    # H : list of the following elements
    # [site of attack, list of all mobile units located before the second stem, number of mobile unit before the second stem, 
    # list of all mobile units located after the second stem, number of mobile unit after the second stem,
    # position of the sequence among all sequences substracted by the sequences without stem or without MU]

def HpFishing(seq, regex, regexca, treatment_type):
    
    H = []
    idx = 0
    nb_NNCA_MU = 0
    
    H_file = open('H_sequence_' + file_to_process + '_type_%s_' %(treatment_type) + '.fasta', "w")
    H_insert_file = open('H_sequence_insert_' + file_to_process + '_type_%s_' %(treatment_type) + '.fasta', "w")
    no_H_file = open('no_H_sequence_' + file_to_process + '_type_%s_' %(treatment_type) + '.fasta', "w")
    multiple_H_file =open('multiple_H_sequence_' + file_to_process + '_type_%s_' %(treatment_type) + '.fasta', "w")
    strict_cond = 0
    loose_cond = 0
    Hp = 0
    Hp_like = 0
    multiple_hp = 0
    no_Hp = 0
    complete_stem_hp = 0
    
    hp_stem = r'(CG)?'+ hairpin_stem2
    regex_repet_hp1 = r'^(AACAT)*%s'%hp_stem
    regex_repet_hp2 = r'^(AACAT)*(AACA[ATCG])?%s'%polyA


    
    for x in seq:
        # Hairpin-like sequences having experienced TSA
        if len(re.findall(hairpin_stem2, x)) > 1 :
            Hp_like += 1
            # Multiple hairpins havong more than 2 stems
            if len(re.findall(hairpin_stem2, x)) > 2 :
                multiple_hp += 1
                multiple_H_file.write(x + '\n')
            # Hairpin with only 2 stem
            if len(re.findall(hairpin_stem2, x)) == 2 :
                Hp += 1
                if len(re.findall(stem, x)) == 2 :
                    complete_stem_hp += 1

                if re.search(regex_repet_hp1, x.split(stem)[1]) and re.search(regex_repet_hp2, x.split(hairpin_stem2)[2]):
                    strict_cond += 1
                    H_file.write(x + '\n')
                        
                    # Find the mobile units inside the hairpin, located between the two stems
                    H1 = re.findall(regex, x.split(stem)[1].split(hairpin_stem2)[0]) 
                    # Case where the UMI are only at the beginning of the read
                    if len(UMI_pos) == 1: 
                    # Find the mobile units after the second stem
                        H2 = re.findall(regex, x.split(hairpin_stem2)[2])
                        # Case where the UMI are located at the beginning and the end of the read
                    elif len(UMI_pos) == 2:
                        H2 = re.findall(regex, x.split(hairpin_stem2)[2][0:len(x.split(hairpin_stem2)[2]) - UMI_pos[1]])

                    if len(H2) != 0 :
                        last_MU = x.split(H2[-1])[-1]
                        find_CA = re.findall(regexca, last_MU)
                        find_CA = [y for y in find_CA if len(y) == MU + 2]
                    else:
                    # Case where the only MU at the end of hp is NNCA ??
                        last_MU_bis = x.split(hairpin_stem2[-1])[-1]
                        find_CA = re.findall(regexca, last_MU_bis) 
                        # or find_CA = re.findall(regexca, (x.split(second_stem)[-1])[-1])
                        find_CA = [y for y in find_CA if len(y) == MU + 2]
                        # We are added only the mobile units that have correct length to H1 and H2 respectively
                    #H1 = [x[0:MU] for z in H1 if len(z) == MU + 3] 
                    #H2 = [x[0:MU] for z in H2 if len(z) == MU + 3] 
                    # Now we are counting in the NNCA mobile units, and adding them in H2 list before adding them to the final H list
                    if len(find_CA) != 0:
                        if find_CA[0] != '' :
                            nb_NNCA_MU += 1
                            H2.append(find_CA[0][0:MU])
                    H.append([H1, len(H1), H2, len(H2), idx])
                else :
                    H_insert_file.write(x + '\n')
        elif len(re.findall(hairpin_stem2, x)) == 1 and re.search('[CG]{3}', x.split(stem)[1]):
            no_Hp += 1
            no_H_file.write(x + '\n')
        idx += 1
        
    H_file.close()
    H_insert_file.close()
    no_H_file.close()
    
    print('The % of sequences finishing by AACAX for Hp is :' + str(nb_NNCA_MU/len(seq)))
    print('There are %s : %s hairpain_like, %s : %s multiple-hp and %s : %s hairpins ' %(Hp_like, str(Hp_like/idx*100), multiple_hp, str(multiple_hp/idx*100), Hp, str(Hp/idx*100)))
    print('Among the %s hps, there are %s : %s complete stem hps, and %s : %s strict hps'%(Hp, complete_stem_hp, str(complete_stem_hp/idx*100), strict_cond, str(strict_cond/idx*100)))
    print('There are also %s : %s not hp'%(no_Hp, str(no_Hp/idx*100)))
    
    return H


# # Sequence_diversity 
def Sequence_diversity(S, H, MU_thresh, H_thresh, MU_in_H_thresh): 
    """"
    Function which count the number of sequences of each type (hairpin, stem with mobile unit, etc...)
    and store the results into an excel file.
    
    
    Parameters:
    ----------
    S : list of the following elements 
    [list of all mobile units of the sequence  (ex : 'AA', 'TG'), number of mobile unit of the sequence, 
    position of the sequence among all sequences substracted by the sequences without stem or without MU]
    
    H : list of the following elements
    [site of attack, list of all mobile units located before the second stem, number of mobile unit before the second stem, 
    list of all mobile units located after the second stem, number of mobile unit after the second stem,
    position of the sequence among all sequences substracted by the sequences without stem or without MU]

    MU_thresh : number of MU for stems only processed 
    H_thresh : H_thresh is the number of combination allowed
    ex : for Hx,y with x<=3 and y<=3, there are 4*4 combinations == 16 
    MU_in_H_thresh : number of MU within hairpins between the two stems processed

    
    Output :
    ---------
    substrate_name.xlsx : excel file 
    """
    
    # STEM S
    #-----------------------------
    Nb_MU_S = [] # list that store the number of mobile unit for each sequence in S
    # X is a list of digit from 0 to MU_thresh-1, will be used as axis when plotting the data
    X_S = [x for x in range(MU_thresh)] 
     
    
    i = 0 
    for element in S:
        # The number of mobile units on the substrate of type Stem is stored in the sub-list in second position ==> S[i][1] 
        # Here, we store in a new list called Nb_MU_S the number of mobile units for every sequences of S
        Nb_MU_S.append(S[i][1]) 
        i += 1
    
    nb = 0
    Prop_S = [] # list that store the number of the sequences having from 0 to MU_threshold mobile units
    for nb in range(MU_thresh):
        Prop_S.append(Nb_MU_S.count(nb)) #count the number of sequence in S with i mobile units
        #print(Nb_MU_S,Prop_S)

    # HAIRPINS all confounded Hxx
    # -------------------------------
    Prop_Hxx = [0 for x in range(H_thresh)]
    # list of label for plotting the data     
    X_Hxx = ['H' + str(x) + ',' + str(y) 
           for x in range(int(math.sqrt(H_thresh))) 
           for y in range(int(math.sqrt(H_thresh)))] 
    
    # for each sequence in H we decide its attribution (Hx,y)
    for element in H:
        i = 0 
        for x in range(int(math.sqrt(H_thresh))): # we check if this sequence is an Hx,y
            for y in range(int(math.sqrt(H_thresh))):
                if element[2] == x and element[4] == y: 
                    # H[i][2] = number of mobile unit before the second stem
                    # H[i][4] = number of mobile unit after the second stem
                    Prop_Hxx[i] += 1
                i += 1   
   
    #Prop_Hxx=[x/len(real_seq) for x in Prop_Hxx]
    
    
    # HAIRPINS H1 and H2
    #--------------------------
    # MU_in_H_thresh = 6
    X_H = [x for x in range(MU_in_H_thresh)]
    
    # 1. AFTER THE FIRST STEM = H1,y
    #---------------------------
    Nb_MU_H1 = []
    i = 0
    for element in H:
        Nb_MU_H1.append(H[i][2])
        i += 1
    
    Prop_H1 = []
    for i in range(MU_in_H_thresh):
        Prop_H1.append(Nb_MU_H1.count(i))#/len(real_seq))
    
    
         
    # 2. MU AFTER THE SECOND STEM Hx,2
    #--------------------------- 
    Nb_MU_H2 = []
    i = 0
    for element in H:
        Nb_MU_H2.append(H[i][4])
        i += 1
        
    Prop_H2 = []
    for i in range(MU_in_H_thresh):
        Prop_H2.append(Nb_MU_H2.count(i))#/len(real_seq))
        
    
    
    # STEM + HAIRPINS
    # ---------------------------
     
    Prop = Prop_S + Prop_Hxx
    X = ['S'+str(x) for x in range(MU_thresh)]
    X.extend(X_Hxx) # add X_H to X
    
    
    
    # SAVING DATA INTO EXCEL FILE
    # ------------------------
        
    # store the results in a panda dataframe to save it in an excel file
    df = pd.DataFrame(Prop, index = X, columns = ['Number of sequences'])
    
    with pd.ExcelWriter(excel_pos, engine = "openpyxl", mode ='a') as writer:  
        df.to_excel(writer, sheet_name = final_file_to_process)

    return X_S, Prop_S, X_Hxx, Prop_Hxx, X_H, Prop_H1, Prop_H2, X, Prop


# # PlotDiversity

# In[16]:


def PlotDiversity(x, y, title, xlabel, ylabel, yscale, file_format, path_to_figure):
    plt.bar(x,y) 
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.yscale(yscale)
    plt.xticks(fontsize=8)
    if file_format == 'png':
        save = path_to_figure + substrate_name + title + '.png'
        plt.savefig(save)
    if file_format =='svg':
        save = path_to_figure + substrate_name + title + '.svg'
        plt.savefig(save)
    plt.close()
    



