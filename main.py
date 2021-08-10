# Imports
from re import sub
from function import *
import os  
import pandas as pd


def transforme():
    
def main():
    print("hello and welcome")
    print("This program will guide you throw its execution")
    if not os.path.isdir("./SequencingData"):
        os.makedirs("./SequencingData")

    path_to_file = "./SequencingData"

    if not os.path.isdir("./Figures"):
        os.makedirs("./Figures")

    path_to_figure = "./Figures"

    if not os.path.isfile("./SequencingData/substrates_list.csv"):
        tmp = pd.DataFrame(columns=["Substrat_name", "Substrat", "Stem"])
        tmp.to_csv("./SequencingData/substrates_list.csv", sep= ",", index= False)

    substratList = pd.read_csv("./SequencingData/substrates_list.csv", sep= ",")
    print("List of registered substrat : ")
    print("Substrat | Stem")
    # print(substratList)
    for index, row in substratList.iterrows():
        print(row["Substrat"] + " | " + row["Stem"])

    substrate_name = str.upper(input("Please enter the name of the substrat in the list to choose it, or type a new one to creat it \n Type 'enter' to continu")).replace(",", ";").replace(" ", "-")
    while not substrate_name:
        substrate_name = str.upper(input("Please enter a valide substrat :")).replace(",", ";").replace(" ", "-")

    exist = substratList.query('Substrat == "' + substrate_name + '"')
    if not exist.empty:
        print("existe")
        stem = exist["Stem"].iloc[0]

    else:
        stem = str.upper(input("Please enter the stem you want to alocate to the new substrat : \n For exemple CGCGAATTAACGCGACAACAT \n and type 'Enter to continue \n Stem :")).replace("U", "T")
        while not re.match("^[ATCG]*$", stem) or not stem :
            print("error : Only nucleotyde (A,T,C,G) are autorised \n Please not that 'U' is automaticly replaced by 'T' \n Please try again")
            stem = str.upper(input("Stem : ")).replace("U", "T")

        substratList = substratList.append(pd.DataFrame([(substrate_name, stem)], columns=["Substrat", "Stem"]), ignore_index= True)
        substratList.to_csv("./SequencingData/substrates_list.csv", sep= ",", index= False)

    hairpin_stem1 = stem[-(len(stem)-2):]
    hairpin_stem2 = stem[-(len(stem)-4):]

    # tag = str.upper(input("Please enter the tag separating your mobil unit \n For exemple 'CAT' \n Tag :")).replace("U", "T")
    # while not re.match("^[ATCG]*$", stem) or not tag:
            # print("error : Only nucleotyde (A,T,C,G) are autorised \n Please not that 'U' is automaticly replaced by 'T' \n Please try again")
            # tag = str.upper(input("Tag : ")).replace("U", "T")

    tag = "CAT"
    
    compl_tag = revcom(tag)
    polyA = (r"A{6,18}")
    T_patch = (r"T{6,18}")
    regex_NN = '[ATGC]*?CAT' # regex used to find all the mobile units 
    regex_NNCA = '[ATGC]*?CA' # regex used to find if there is a last mobile unit which tag has been cut

    MU = input("Please enter the length of the mobil unit (without the tag) you want to search : \n (Type 0 if you want to search for the mobil unit of any length) \n MU length :")
    while not re.match("^[0-9]*$", MU) or not MU:
        print("Error : Only number are allowed. \n Please try again")
        MU = input("Length of mobile unit : ")
        
    if int(MU) > 0:
        regex_MU = input("Please enter the value of the full mobil unit you want to search for")
        regex_MU_cuted = input("Please enter the same mobile unit as before, but without the last T(")
        regex_repet_stem = r'^(' + re.escape(regex_MU) + ')*(' + re.escape(regex_MU_cuted) + ')?%s'%polyA
        print(regex_repet_stem)








if __name__ == '__main__':
    main()