#!/usr/bin/python2.7
# Import modules for CGI handling


import pandas as pd
import re
import os
import argparse
#import glob
#import string
#import urllib
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",help = "Input file")
    parser.add_argument("--property", help = """Enter 'AAP' for Amino acid properties""")
#    parser.add_argument("--SUB", help = "Substitution matrices")
#    parser.add_argument("--CON", help = "Pairwise properties and contact potential")
#    parser.add_argument("--ALL", help = "ALL properties")
    args = parser.parse_args()
newDF_49 = pd.DataFrame()
newDF_491 = pd.DataFrame()
f1 = open(args.input_file,'r').readlines()

#print(f1)
################ Left right pref#################
def user_prop(mutation,AA_dist):
    mut = mutation[-1]
    wild = mutation[0]
    return(int(AA_dist[wild])-int(AA_dist[mut]))


def left_pref(seq , mutation):
    pos =int(''.join(re.findall('\d+',mutation)))
    if pos <=1 or pos>= len(seq):
        left_pref = ''
        right_pref = ''
        return(left_pref,right_pref)
        out_df = pd.DataFrame({'Mutation':[mutation],'N_Terminal':[left_pref],'C_Terminal':[right_pref]},columns= ['Mutation','N_Terminal','C_Terminal'])
    else:
        left_pref = seq[pos-2]
        right_pref = seq[pos]
        out_df = pd.DataFrame({'Mutation':[mutation],'N_Terminal':[left_pref],'C_Terminal':[right_pref]},columns= ['Mutation','N_Terminal','C_Terminal'])
        return(out_df)

def prop_49(seq ,mutation):
    df1 = pd.read_csv('./data/49_properties_numerical_Values.csv', sep =',' )
#    file1 = pd.read_csv('./data/prop_49_list.csv')
    mut = mutation[-1]
    wild = mutation[0]
    wrt = (df1[wild])
#    out = pd.concat([file1,wrt],axis =1)
    return(wrt.transpose())
def prop_491(seq ,mutation):
    df_ = pd.read_csv('./data/49_properties_normalizedValues.csv', sep =',' )
#    file1 = pd.read_csv('./data/prop_49_list.csv')
    mut = mutation[-1]
    wild = mutation[0]
    wrt1 = (df_[wild])
#    out = pd.concat([file1,wrt],axis =1)
    return(wrt1.transpose())
path = "./out_protein"
directory = os.path.dirname(path)
if not os.path.exists(directory):
    os.makedirs(directory)
newDF = pd.DataFrame()
newDF_ = pd.DataFrame()
i = 1
for li in f1:
    #print("calculating for seq no" +str(i))
    i = i+1
    # Get data from fields
    seq = li.strip().split('\t')[0].replace('X','').replace('Z','')
    mutation_list = []
    AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    #mutation_list = form.getvalue('mutation').split(",")

    for number in range(len(seq)):
        mutation_list.append(seq[number])
#        for aa in AA:
#            if seq[number]+str(number)+aa not in mutation_list:
#                mutation_list.append(seq[number]+str(number)+aa)
    se = pd.Series(mutation_list)
    df_mut_list= pd.DataFrame()
    for mut in mutation_list:
        if(len(mut)<3):
            pass
#                continue

        mutation = mut.upper()

        if(str(args.property) == "AAP"):
#                    print("calculating... amino acid properties")
            pf = prop_49(str(seq) , str(mutation))
            pf1 = prop_491(str(seq) , str(mutation))
            newDF_49 = newDF_49.append(pf, ignore_index=True)
            newDF_491 = newDF_491.append(pf1, ignore_index=True)
        else:
            print("!!!Wrong choice!!!")
            print("Please select at least one property ...!!!")
        #print"%s"%(newDF)

    if(str(args.property) == "AAP"):
#        cols = ['Properties']
#        for ij in mutation_list:
#            cols.append(ij)
#    #
        file1 = pd.read_csv('./data/prop_49_list.csv').transpose()
        newDF1 = pd.concat([file1,newDF_49],axis =0)
        newDF1_1 = pd.concat([file1,newDF_491],axis =0)
        new = newDF1.transpose()
        new1 = newDF1_1.transpose()
        new = pd.concat([file1,newDF_49.mean(axis = 0).to_frame(name=None).transpose()]).transpose()
        new1 = pd.concat([file1,newDF_491.mean(axis = 0).to_frame(name=None).transpose()]).transpose()
        newDF = newDF.append(new.T,ignore_index = False)
        newDF_ = newDF_.append(new1.T,ignore_index = False)
if(str(args.property) == "AAP"):
        newDF.drop_duplicates(keep= 'first').to_csv("./out_protein/avg_prop_change_AAP.csv",index = False,header = False)
        newDF_.drop_duplicates(keep= 'first').to_csv("./out_protein/normalizedValues_avg_prop_change_AAP.csv",index = False,header = False)
