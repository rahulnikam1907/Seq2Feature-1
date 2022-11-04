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
    parser.add_argument("--property", help = """Enter 'AAP' for Amino acid properties,'SUB' for Substitution matrices, 'CON' for Pairwise properties and contact potential or 'ALL' for all properties""")
#    parser.add_argument("--SUB", help = "Substitution matrices")
#    parser.add_argument("--CON", help = "Pairwise properties and contact potential")
#    parser.add_argument("--ALL", help = "ALL properties")
    args = parser.parse_args()
newDF_49 = pd.DataFrame()
newDF_AA = pd.DataFrame()
newDF_con = pd.DataFrame()
newDF_neg = pd.DataFrame()
newDF_user_prop = pd.DataFrame()
new_user = []
f1 = open(args.input_file,'r').readlines()[1:]
################ Left right pref#################
def user_prop(mutation,AA_dist):
    mut = mutation[-1]
    wild = mutation[0]
    return(int(AA_dist[mut])-int(AA_dist[wild]))
    
    
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
    wrt = (df1[mut]-df1[wild])
#    out = pd.concat([file1,wrt],axis =1)
    return(wrt.transpose())
   
#####contact potential 47*2###############

def contact_potential(seq, mutation):
    pos =int(''.join( re.findall('\d+',mutation)))
    left = seq[pos-2]
    right = seq[pos]
    dN_wild = mutation[0] + left
    dN_mutant = mutation[-1]+ left
    dC_wild = mutation[0] + right
    dC_mutant = mutation[-1]+ right
    df1 = pd.read_csv('./data/data_contact_potential_diagonal.csv')
    if dN_wild not in df1:
        dN_wild = left + mutation[0]
    if dN_mutant not in df1:
        dN_mutant = left + mutation[-1]
    if dC_wild not in df1:
        dC_wild = right + mutation[0]
    if dC_mutant not in df1:
        dC_mutant = right + mutation[-1]
    potential_dN1 = df1[dN_mutant] -df1[dN_wild]
    potential_dC1 = df1[dC_mutant]- df1[dC_wild]
    file2 = pd.read_csv('./data/data_contact_potential_square.csv')
    potential_dN =  file2[dN_mutant] -file2[dN_wild]
    potential_dC =  file2[dC_mutant] -file2[dC_wild]
#    file3 = pd.read_csv('./data/contact_potential_prop.csv',header = None)
#    index = ["ZHAC000102","ZHAC000103","ZHAC000105","ZHAC000102","ZHAC000103","ZHAC000105","BASU010101","BETM990101","BONM030101","BONM030102","BONM030103","BONM030104","BONM030105","BONM030106","BRYS930101","GODA950101","KESO980101","KESO980102","KOLA930101","LIWA970101","MICC010101","MIRL960101","MIYS850102","MIYS850103","MIYS960101","MIYS960102","MIYS960103","MIYS990106","MIYS990107","MOOG990101","PARB960101","PARB960102","ROBB790102","SIMK990101","SIMK990102","SIMK990103","SIMK990104","SIMK990105","SKOJ000101","SKOJ000102","SKOJ970101","TANS760101","TANS760102","THOP960101","TOBD000101","TOBD000102","VENM980101","ZHAC000101","ZHAC000104","ZHAC000106","BASU010101","BETM990101","BONM030101","BONM030102","BONM030103","BONM030104","BONM030105","BONM030106","BRYS930101","GODA950101","KESO980101","KESO980102","KOLA930101","LIWA970101","MICC010101","MIRL960101","MIYS850102","MIYS850103","MIYS960101","MIYS960102","MIYS960103","MIYS990106","MIYS990107","MOOG990101","PARB960101","PARB960102","ROBB790102","SIMK990101","SIMK990102","SIMK990103","SIMK990104","SIMK990105","SKOJ000101","SKOJ000102","SKOJ970101","TANS760101","TANS760102","THOP960101","TOBD000101","TOBD000102","VENM980101","ZHAC000101","ZHAC000104","ZHAC000106"]
    out =  pd.concat([potential_dN, potential_dC,potential_dN1, potential_dC1],ignore_index = True)
    df_out = out.to_frame()
#    con = pd.concat([file3,df_out],axis = 1)
    return df_out.transpose()
    ####################### AA 94 properties######################
def AAIndex_94(mutation):
    file1 = pd.read_csv('./data/94_AA_Index.csv')
    mut1 = mutation[0]+mutation[-1]
    out = file1[file1['mutation'].str.contains(mut1, na = False)]
#    print(out.loc[:,'ALTS910101':])    
#    print(type(out))    
    return out
path = "./out_protein"
directory = os.path.dirname(path)
if not os.path.exists(directory):
    os.makedirs(directory)
newDF = pd.DataFrame()
clas = []
for li in f1:
    # Get data from fields
    seq = li.strip().split('\t')[0]
    mutation_list = []
    temp = li.strip().split('\t')[1]
    temp = ''.join(temp.split())
    clas.append(li.strip().split('\t')[2])
    #seq1=  seq.translate(string.maketrans("\t\r", "  "))
    #var =  seq1.splitlines()
    #flag = False
    #seq = ''
    #for lin in var:
    #    if(flag):
    #        seq += lin.strip()
    #    if lin.startswith(">"):
    #        flag = True
    #flag = False
    ##print seq
    AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']
    #mutation_list = form.getvalue('mutation').split(",")
    
#    number = random.randint(0,9)
    #prop_1  = 'properties'
    #num_value1 = form.getvalue('num_value1')
        
    
    
    if(True):
    
    #    seq = ''
    #    fileitem = ''
    #    #input file
    #    fileitem = 'filename'
    #    if(len(seq1)>1 or fileitem.filename):
    #        if(len(seq1)>1):       
    #            var = seq1.splitlines()
    #            # print("%s")%var
    #            flag = False
    #            for lin in var:
    #                if(flag):
    #                    seq += lin.strip()
    #                if lin.startswith(">"):
    #                    flag = True
    #    #        print("%s")%lin
    #            flag = False
    #        elif(fileitem.filename):
    #            
    #            fn = os.path.basename(fileitem.filename.replace("\\", "/" ))
    #            var = fileitem.file.readlines()
    #            flag = False
    #            for line in var:
    #                if(flag):
    #                    seq += line.strip()
    #                if ">" in line:
    #                    flag = True
    #            flag = False
    #        else:
    #            print("check seq")
        for number in range(len(seq)):
            for aa in AA:
                if seq[number]+str(number)+aa not in mutation_list:
                    mutation_list.append(seq[number]+str(number)+aa)
        se = pd.Series(mutation_list)
        df_mut_list= pd.DataFrame()
        df_mut_list["mutation"] = se.values
        #print("%s")%(df_mut_list.to_html())
        
        
    #    print"%s",num_value1
        #mutation input
        
        mutation_list = temp.split(",")
        #mutation file input
    #    fileitem1 = 'filename1'
    #    var1 = fileitem1.file.readlines()
        
    #    if(len(mutation_list)):
    #        
    #        
    #        if(fileitem1.filename):
    #            for i in var1:
    #                if i.strip() not in mutation_list:
    #                    mutation_list.append(i.strip())
    #        else:
    #            print("Enter mutation with position")
          
        for mut in mutation_list:
            if(len(mut)<3):
                continue
            
            mutation = mut.upper()
        #    print(mutation)
            
            
            #prop1 = "49_amino_acid_Physicochemical_properties"
        #    p = type(prop_1)
            # Get filename here.
        #    fileitem = form['filename']
            
            
                
                
                
    #            print("""<!-- content -->
    #             <div class="wrapper row2">
    #             <div id="container" class="clear">
    #             <div id="intro">
    #        	    </br></br>
    #        	   <section class="clear">""")
            
        #    print("%s")%seq
            #print"""<td bgcolor = "#eee" width = "100" height = "200">"""
            if(len(seq) <= int(mutation[1:-1])-1):
                print("!!! Short sequence length !!!")
                print("!!! Kindly check the sequence and mutation!!!%s"%(seq))
                exit()
            elif(seq[int(mutation[1:-1])-1] != mutation[0]):
                 print("At given position mutant not found...!!!")
                 print("Kindly check the respective position of mutation...!!!%s"%(mutation))
                 exit()
            else:
                
                if(str(args.property) == "AAP"):
#                    print("calculating... amino acid properties")
                    pf = prop_49(str(seq) , str(mutation))
                    newDF_49 = newDF_49.append(pf, ignore_index=True)
        #                prop1 = prop.transpose()
        #                print"%s"% (prop.to_html(header = False,index = False))
        #                prop.to_csv("../../html/SBFE/out_protein/out_prop_49.csv",header = False ,index = False,index_label = False)
        #                print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/out_prop_49.csv';">""")
                elif(str(args.property) == "SUB"):
#                    print("calculating...Sunstitution matrix properties")
                    AA_Index = AAIndex_94(str(mutation))
                    newDF_AA = newDF_AA.append(AA_Index, ignore_index=True)
        #                print("<strong>amino acid index properties:</strong>")
        #                print("%s")%(newDF.to_html(header = True,index = False))
        #                AA_Index.to_csv("../../html/SBFE/out_protein/out_AA_Index.csv",index = False)
        #                print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/out_AA_Index.csv';">""")
                elif(str(args.property) =="CON"):
#                    print("calculating...contact potential")
                    contact_p = contact_potential(str(seq), str(mutation))
                    newDF_con = newDF_con.append(contact_p, ignore_index=True)
        #                print"%s"%(newDF.to_html(header = False, index = False))
        #                contact_p.to_csv("../../html/SBFE/out_protein/out_contact_potential.csv",header = False ,index = False,index_label = False)
        #                print ("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/out_contact_potential.csv';">""")
                elif(str(args.property) == "neighboring residues of mutation position"):
        #                print("finding...neighboring residues")
                    left_pre = left_pref(str(seq) , str(mutation))
                    newDF_neg = newDF_neg.append(left_pre, ignore_index=True)
        #                print("%s"%(newDF.to_html(header = True, index = False)))
        #            print( """<button type="submit" onclick="window.open('file.doc')">Download!</button>""")
        #                left_pre.to_csv("../../html/SBFE/out_protein/neighbour.csv",index_label = False,index = False)
        #                print( """<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/neighbour.csv';">""")
                elif(str(args.property) =="ALL"):
#                    print("calculating...all properties")
                    contact_p = contact_potential(str(seq), str(mutation))
                    newDF_con = newDF_con.append(contact_p, ignore_index=True)
                    AA_Index = AAIndex_94(str(mutation))
                    newDF_AA = newDF_AA.append(AA_Index.loc[:,:], ignore_index=True)
                    prop = prop_49(str(seq) , str(mutation))
                    newDF_49 = newDF_49.append(prop, ignore_index=True)
                    left_pre = left_pref(str(seq) , str(mutation))
                    newDF_neg = newDF_neg.append(left_pre, ignore_index=True)
        #                print("<strong>49_properties:</strong>")        
        #                print("%s")% (prop.to_html(header = False, index = False))
        #                print("<strong>contact potential:</strong>")
        #                print("%s")% (contact_p.to_html(header = False, index = False))
        #                print("<strong>AA_Index_94 properties:</strong>")       
        #                print("%s")% (AA_Index.to_html(header = True, index = False))
        #                print("<strong>neighboring residues of mutation position:</strong>")        
        #                print("%s")% (left_pre.to_html(header = True , index = False))
        #                contact_p.to_csv("../../html/SBFE/out_protein/out_contact_potential.csv",header = False ,index = False,index_label = False)
        #                prop.to_csv("../../html/SBFE/out_protein/out_prop_49.csv",header = False ,index = False,index_label = False)
        #                AA_Index.to_csv("../../html/SBFE/out_protein/out_AA_Index.csv",index = False)
        #                left_pre.to_csv("../../html/SBFE/out_protein/neighbour.csv",index_label = False,index = False)
        #                cmd = "paste -d, "+"../../html/SBFE/out_protein/*.csv"+ " > ../../html/out.csv"
        #                os.system(cmd)
        #                print ("""<input type="button" value="Download Now!" onclick="window.location = '../../out.csv';">""")
            
                    
                else:
                    print("!!!Wrong choice!!!")
                    print("Please select at least one property ...!!!")
        #print"%s"%(newDF)
class1 = pd.Series(clas)
if(str(args.property) == "AAP"):
    cols = ['Properties']
    for ij in mutation_list:
        cols.append(ij)
#        print('%s')%type(cols)
#        print("""<p style = "color:#660214;"><strong>Your input seq is:</strong>""")
#        print("""<div style = "word-wrap: break-word;" > %s</div>""")%(seq)
#        print("""<p style = "color:#660214;"><strong>Amino acid properties:</strong>""")
#        print("""<p align = "justify"><br>We have considered physical, chemical, energetic and conformational <a href = "../../SBFE/49_prop.html">properties</a> and <a href = "../../SBFE/81_AA_property.html">list of properties</a> from AAindex database  <br> The change in property value upon mutation is calculated using following formula:<br><div class="boxed"><p style="text-align:center"><b> &Delta; P<sub>(mutation)</sub>   =   P<sub>(Snp)</sub> - P<sub>(wild-type)</sub></b><br>Where  P<sub>(wild-type)</sub> and P<sub>(Snp)</sub> are the property values of wild-type and mutant residues, respectively and &Delta;P<sub>(mutation)</sub> is the change in property due to mutation</div>
#""")
    file1 = pd.read_csv('./data/prop_49_list.csv').transpose()
    col = open('./data/prop_49_list.csv','r').readlines()
    col1 = []
    for ii in col:
        col1.append(ii.strip())
    newDF1 = pd.concat([file1,newDF_49],axis =0)
    new = newDF1.transpose()
    
#    print(newDF1)
#    col1 = file1.values
#    col1.append('class')
#    print(col1)
    newDF1.columns = col1
    newDF1['Class'] = class1
#            new.columns = cols
#        print"%s"% (new.to_html(columns= cols,index = False))
    newDF1[1:].round(3).to_csv("./out_protein/AAP_out.csv",header = True , index = False,index_label = False)
#        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/"""+ str(number)+"""out_prop_49.csv';">""")
elif(str(args.property) == "SUB"):
#        print("""<p style = "color:#660214;"><strong>Your input seq is:</strong>""")
#        print("""<div style = "word-wrap: break-word;" > %s</div>""")%(seq)
#        print("""<p style = "color:#660214;"><strong>Substitution matrix properties:</strong>""")
#        print("""<p align = "justify"><br>Amino acid mutation <a href = "../../SBFE/substitution.html">matrices</a> are collected from AAIndex2 database and the mutation value is directly obtained from matrices.
#""")
    cols = []
    for ij in mutation_list:
        cols.append(ij)
    new = newDF_AA.transpose()[1:]
    newDF_AA['Class'] = class1
#            new.columns = cols
#        print("%s")%(new.to_html())
    newDF_AA.round(3).to_csv("./out_protein/SUB_out.csv",index = False)
#        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/"""+ str(number)+"""out_AA_Index.csv';">""")
elif(str(args.property) =="CON"):
#        print("""<p style = "color:#660214;"><strong>Your input seq is:</strong>""")
#        print("""<div style = "word-wrap: break-word;" > %s</div>""")%(seq)
    file3 = pd.read_csv('./data/contact_potential_prop.csv',header = None).transpose()
#    col = open('./data/contact_potential_prop.csv','r').reaadlines()
#    newDF_con['Class'] = class1
    con = pd.concat([file3,newDF_con],axis = 0)
    
#        print("""<p style = "color:#660214;"><strong>Contact potential properties:</strong>""")
#        print("""<p  align = "justify">Pair-wise contact potential <a href = "../../SBFE/contact_p.html">matrices</a> are collected from AAIndex3 database and difference of amino acid contact potential for a mutation is obtained by subtracting contact potential value of N-/C-neighbour of mutation position to wild type residue from N-/C-neighbour to mutant residue. 
#""")
    cols = ['Properties']
    for ij in mutation_list:
        cols.append(ij)
    new = con.transpose()
    con['Class'] = class1
#            new.columns = cols
#        print"%s"%(new.to_html( index = False))
#        print("Here '_N' and '_C' represents N-terminal and C-terminal Amino acid respectively")
    con.round(3).to_csv("./out_protein/out_contact_potential.csv",header = False ,index = False,index_label = False)
#        print ("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/"""+ str(number)+"""out_contact_potential.csv';">""")
elif(str(args.property) == "neighboring residues of mutation position"):
#    print("finding...neighboring residues")
#        print("%s"%(newDF_neg.transpose().to_html(header = False , index = True)))
    newDF_neg.to_csv("./out_protein/neighbour.csv",index_label = False,index = False)
#        print( """<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_protein/"""+ str(number)+"""neighbour.csv';">""")
elif(str(args.property) =="ALL"):
#    print("calculating...all properties")
#    print("49_properties:")
    file1 = pd.read_csv('./data/prop_49_list.csv').transpose()
    newDF1 = pd.concat([file1,newDF_49],axis =0)
    col = open('./data/prop_49_list.csv','r').readlines()
#        print"%s"% (newDF1.to_html(header = False,index = False))
#    newDF1.round(3).to_csv("./out_protein/AAP_out.csv",header = False ,index = False,index_label = False)
#    print("amino acid index properties:")
#        print("%s")%(newDF_AA.to_html(header = True,index = False))
    col1 = []
    for ii in col:
        col1.append(ii.strip())
    newDF1 = pd.concat([file1,newDF_49],axis =0)
    new = newDF1.transpose()
    
#    print(newDF1)
#    col1 = file1.values
#    col1.append('class')
#    print(col1)
    newDF1.columns = col1
    newDF1['Class'] = class1
#            new.columns = cols
#        print"%s"% (new.to_html(columns= cols,index = False))
    newDF1[1:].round(3).to_csv("./out_protein/AAP_out.csv",header = True , index = False,index_label = False)
    newDF_AA.drop(columns = 'mutation',inplace = True)
    newDF_AA.round(3).to_csv("./out_protein/SUB_out.csv",index = False)
#    print("Contact potential:")
    file3 = pd.read_csv('./data/contact_potential_prop.csv',header = None).transpose()
    con = pd.concat([file3,newDF_con],axis = 0)
#        print"%s"%(con.to_html(header = False, index = False))
#        print("Here '_N' and '_C' represents N-terminal and C-terminal Amino acid respectively")
    
    con.round(3).to_csv("./out_protein/out_contact_potential.csv",header = False ,index = False,index_label = False)
    
#    df111 = pd.DataFrame()
#    df111 = pd.concat([newDF1[1:].round(3),newDF_AA.round(3),con.round(3)],axis = 1,ignore_index = True)
#    df = df111.fillna(0)
#    df.to_csv("./output.csv",index = False)
#        print("<strong>neighboring residues of mutation position:</strong>")
#        print("%s"%(newDF_neg.to_html(header = True, index = False)))
#    con['Class'] = class1
#    newDF_neg['Class'] = class1
#    newDF_neg.to_csv("./out_protein/neighbour.csv",index_label = False,index = False)
    cmd = "paste -d, "+"./out_protein/*.csv"+ " > ./output.csv"
    os.system(cmd)


    
   
