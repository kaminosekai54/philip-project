# Imports
from function import *
import os  
import pandas as pd


def main():
    print("hello and welcome")
    print("This program will guide you throw its execution")
    settings = loadSettings("./settings.json")
    if settings["workingDirectory"]:
        changeDir = str.upper(input("Your working directory is set as : " + settings["workingDirectory"] + "\n Do you want to change it ? \n Type Y / N"))
        while not re.match("^[YN]*$", changeDir) or not changeDir:
            print("Please enter Y to change your current working directory and N to keep it")
            changeDir = str.upper(input("Your working directory is set as : " + settings["workingDirectory"] + "\n Do you want to change it ? \n Type Y / N"))
        if changeDir == "Y":
            settings["workingDirectory"]= input("Please enter your working directory (path to the directory where your sequence file are stored) ")
            saveSettings(settings, "./settings.json")
    else:
        settings["workingDirectory"] = input("Please enter your working directory (path to the directory where your sequence file are stored) ")
        saveSettings(settings, "./settings.json")

    path_to_file = settings["workingDirectory"]
    # get the folder where to store the results , if not, creat it
    if not os.path.isdir(settings["workingDirectory"] + "Results"):
        os.makedirs(settings["workingDirectory"] + "Results")
    excel_pos = settings["workingDirectory"] + "Results/"
# get the folder where the figures are stored, if not, creat it
    if not os.path.isdir(settings["workingDirectory"] + "Figures"):
        os.makedirs(settings["workingDirectory"] + "Figures")
    path_to_figure = settings["workingDirectory"] + "Figures/"
# get the file where the stem/substrat info are stored, if not, creat it
    if not os.path.isfile("./Settings/substrates_list.csv"):
        tmp = pd.DataFrame(columns=["Substrat_name", "Substrat", "Stem"])
        tmp.to_csv("./Settings/substrates_list.csv", sep= ",", index= False)
    substratList = pd.read_csv("./Settings/substrates_list.csv", sep= ",")

    # get the position of mothur 
    # position_mothur 


    print("List of registered substrat : ")
    print("Substrat | Stem")
    # print(substratList)
    for index, row in substratList.iterrows():
        print(row["Substrat_name"] + " | " + row["Stem"])

    substrate_name = str.upper(input("Please enter the name of the substrat in the list to choose it, or type a new one to creat it \n Type 'enter' to continu")).replace(",", ";").replace(" ", "-")
    while not substrate_name:
        substrate_name = str.upper(input("Please enter a valide substrat :")).replace(",", ";").replace(" ", "-")

    exist = substratList.query('Substrat_name == "' + substrate_name + '"')
    if not exist.empty:
        print("existe")
        stem = str.replace(str.upper(exist["Stem"].iloc[0]), "U", "T")

    else:
        new_substrat= str.upper(input("Please enter the Full substratyou want to alocate to the new substrat name  (Stem will be automatically extract from it): \n For exemple CGCGAATTAACGCGACAACATAACAT \n and type 'Enter to continue \n Stem :")).replace("U", "T")
        while not re.match("^[ATCG]*$", stem) or not stem :
            print("error : Only nucleotyde (A,T,C,G) are autorised \n Please not that 'U' is automaticly replaced by 'T' \n Please try again")
            new_substrat= str.upper(input("Stem : ")).replace("U", "T")

        stem = new_substrat.split("CAT")[0] + "CAT"
        substratList = substratList.append(pd.DataFrame([(substrate_name, new_substrat, stem)], columns=["Substrat_name", "Substrat", "Stem"]), ignore_index= True)
        substratList.to_csv("./Settings/substrates_list.csv", sep= ",", index= False)

    compl_stem = revcom(stem)
    stem_RC = r"ATG%s" % (revcom(stem))
    hairpin_stem1 = stem[-(len(stem)-2):]
    hairpin_stem2 = stem[-(len(stem)-4):]

    tag = "CAT"
    # uncommant this if you want to choose a different tag
    # tag = str.upper(input("Please enter the tag separating your mobil unit \n For exemple 'CAT' \n Tag :")).replace("U", "T")
    # while not re.match("^[ATCG]*$", tag) or not tag:
    #         print("error : Only nucleotyde (A,T,C,G) are autorised \n Please not that 'U' is automaticly replaced by 'T' \n Please try again")
            # tag = str.upper(input("Tag : ")).replace("U", "T")
    
    compl_tag = revcom(tag)
    polyA = (r"A{6,18}")
    T_patch = (r"T{6,18}")

    MU = input("Please enter the length of the mobil unit (without the tag) you want to search : \n (Type 0 if you want to search for the mobil unit of any length) \n MU length :")
    while not re.match("^[0-9]*$", MU) or not MU:
        print("Error : Only number are allowed. \n Please try again")
        MU = input("Length of mobile unit : ")
        
    if int(MU) > 0:
        regex_MU = str.upper(input("Please enter the value of the mobil unit you want to search for \n The tag " + tag + " will be added at the end. \n For exemple type AA if your full mobile unit is AACAT")).replace("U","T")
        while not re.match("^[ATCG]*$", regex_MU ) or not regex_MU :
            print("error : Only nucleotyde (A,T,C,G) are autorised \n Please not that 'U' is automaticly replaced by 'T' \n Please try again")
            regex_MU = str.upper(input("Please enter the value of the mobil unit you want to search for \n The tag " + tag + " will be added at the end. \n For exemple type AA if your full mobile unit is AACAT")).replace("U","T")

        regex_MU = regex_MU + tag
    else:
        regex_MU= '[ATGC]*?%s'%tag # regex used to find all the mobile units 

    regex_MU_cuted = regex_MU [:-1]
    regex_repet_stem = r'^(' + re.escape(regex_MU) + ')*(' + re.escape(regex_MU_cuted) + ')?%s'%polyA

    # option for the treatement
    concat = 0 # 0 is no concat or 1 if concatenation
    raw = -1 # 1 : raw file (fastQ without cleaning); -1 : file after treatment (cleaning, concatenation, etc..)
    treatment_type = input("Please enter the type of treatement you want to process the sequences : n 1 : removal of sequences that don't have GGG; \n 2 : removal of sequences that don't have [CG] \n 3 : if any") 
    while not re.match("^[1-3]*$", treatment_type) or not treatment_type:
        print("Error : Only the number mentionned above are allowed. \n Please try again")
        treatment_type = input("Please enter the type of treatement you want to process the sequences : \n 1 : removal of sequences that don't have GGG; \n 2 : removal of sequences that don't have [CG] \n 3 : if any ")

    UMI_pos = [] # Position and size of the UMI:  UMI_pos=[x] if the UMI is the x first nucleotides 
    # and UMI_pos=[x,y] if the UMI is the x first nucleotides and the y last
    umi_start = input("Please enter the last position of the first UMI")
    while not re.match("^[0-9]*$", umi_start) or not umi_start:
        print("Error : Only number are allowed. \n Please try again")
        umi_start = input("Please enter the last position of the first UMI")

    UMI_pos.append(int(umi_start))
    umi_end= input("Please enter the first position of the second part of the UMI, type 0 if your UMI is in one part only (must be > than the start and can not be equal")
    while not re.match("^[0-9]*$", umi_end) or not umi_end :
        print("Error : Only number are allowed or the value you entered don't respect the above constrinct. \n Please try again")
        umi_end= input("Please enter the first position of the second part of the UMI, type 0 if your UMI is in one part only (must be > than the start and can not be equal")

    if int(umi_end) > 0:
        if int(umi_end) <= int(umi_start):
            print("Error : the first position of the second UMI must be > than the last position of the first UMI |n Must be > than : " + umi_start )
            while not re.match("^[0-9]*$", umi_end) or not umi_end or int(umi_end) <= int(umi_start):
                print("Error : the first position of the second UMI must be > than the last position of the first UMI |n Must be > than : " + umi_start )
                umi_end = input("Please enter the first position of the second UMI ( Must be > than the last position of the first UMI")

        UMI_pos.append(int(umi_end))

# config for the plot
    Nb_MU = 6 #number of MU plotted
    MU_thresh = 26 # number of MU for stems only processed by the SequenceDiversity fonction
    H_thresh = 16 # H_thresh : H_thresh is the number of combination allowed by the SequenceDiversity fonction
    # ex : for Hx,y with x<=3 and y<=3, there are 4*4 combinations == 16 
    MU_in_H_thresh = 6

    #Calling functions
    output_path, file_to_process = Preprocessing([path_to_file, excel_pos, substrate_name, stem, stem_RC, MU, UMI_pos, Nb_MU, polyA, T_patch, concat, raw])
    seq, final_file_to_process = UnexpectedSeqRemover(output_path, file_to_process, treatment_type, stem, tag)
    S, NUS_Stem = SFishing(output_path, final_file_to_process , regex_repet_stem, seq, regex_MU, regex_MU_cuted, treatment_type, stem, UMI_pos, int(MU), hairpin_stem1, hairpin_stem2)
    H = HpFishing(output_path, final_file_to_process , regex_repet_stem, seq, regex_MU, regex_MU_cuted, treatment_type, stem, UMI_pos, int(MU), hairpin_stem1, hairpin_stem2, polyA)
    X_S, Prop_S, X_Hxx, Prop_Hxx, X_H, Prop_H1, Prop_H2, X, Prop = Sequence_diversity(output_path, substrate_name, final_file_to_process , S, H, MU_thresh, H_thresh, MU_in_H_thresh)



if __name__ == '__main__':
    main()