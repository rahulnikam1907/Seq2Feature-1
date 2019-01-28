#!/usr/bin/python
import pandas as pd
import re
import sys,os,errno
import glob
def prog(inFile,mutation):
#    inFile = sys.argv[1]
#    mutation = sys.argv[2]
    #outFile = sys.argv[3]
    seq = ''
    if(inFile.strip().split(".")[1]=="csv"):
        with open(inFile,'r') as j:
            d = j.readlines()
            for line1 in d:
                l_list = line1.strip("\n").split(",")
                mut_with_pos = l_list[0].strip("\n")
                seq = l_list[1].strip("\n")
    elif(inFile.strip().split(".")[1]=="fasta"):
        with open(inFile,'r') as i:
            r = i.readlines()
            flag = False
            for line in r:
                if(flag):
                    seq += line.strip()        
                if ">" in line:
                    flag = True
    else:
        with open(inFile,'r') as i:
            r = i.readlines()
            for line in r:
                seq = seq +line.strip()
    # 
    #with open(inFile,'r') as i:
    #    r = i.readlines()
    #    for line in r:
    #        seq = seq +line.strip()
    #print seq
    ##print mutation
    ##
    #seq = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD'
    #mutation = 'G361R'
    
    if(len(seq) <= int(mutation[1:-1])-1):
        print("!!! Short sequence length !!!")
        print("!!! Kindly check the sequence and mutation!!!")
    elif(seq[int(mutation[1:-1])-1] != mutation[0]):
         print("At given position mutant not found...!!!")
         print("Kindly check the respective position of mutation...!!!")
    
    else:
        ################ Left right pref#################
        def left_pref(seq , pos):
            if pos <=1 or pos>= len(seq):
                left_pref = ''
                right_pref = ''
                return(left_pref,right_pref)
                out_df = pd.DataFrame({'left_neighbour':[left_pref],'right_neighbour':[right_pref]})
            else:
                left_pref = seq[pos-2]
                right_pref = seq[pos]
                out_df = pd.DataFrame({'left_neighbour':[left_pref],'right_neighbour':[right_pref]})
                return(out_df)
        #        print out_df
            
        def prop_49(seq ,mutation):
            df1 = pd.read_csv('./data/49_properties_normalizedValues.csv', sep ='\t' )
            file1 = pd.read_csv('./data/prop_49_list.csv')
            mut = mutation[-1]
            wild = mutation[0]
            wrt = (df1[mut]-df1[wild])
            out = pd.concat([file1,wrt],axis =1)
        #    out_df =  out.to_frame()
            return(out.transpose())
        #    print type(out)
        
        
        #####contact potential 47*2###############
        
        def contact_potential(seq, mut_with_pos):
            mutation = mut_with_pos[0] + mut_with_pos[-1]
            pos =int(''.join( re.findall('\d+',mut_with_pos)))
            left = seq[pos-2]
            right = seq[pos]
            dN_wild = mut_with_pos[0] + left
            dN_mutant = mut_with_pos[-1]+ left
            dC_wild = mut_with_pos[0] + right
            dC_mutant = mut_with_pos[-1]+ right
            file1 = pd.read_csv('./data/data_contact_potential_diagonal.csv')
            if dN_wild not in file1:
                dN_wild = left + mut_with_pos[0]
            if dN_mutant not in file1:
                dN_mutant = left + mut_with_pos[-1]
            if dC_wild not in file1:
                dC_wild = right + mut_with_pos[0]
            if dC_mutant not in file1:
                dC_mutant = right + mut_with_pos[-1]
            potential_dN1 = file1[dN_wild] - file1[dN_mutant]
            potential_dC1 = file1[dC_wild] - file1[dC_mutant] 
            file2 = pd.read_csv('./data/data_contact_potential_square.csv')
            potential_dN = file2[dN_wild] - file2[dN_mutant]
            potential_dC = file2[dC_wild] - file2[dC_mutant] 
        #    return(potential_dN, potential_dC,potential_dN1, potential_dC1)
        #    print potential_dN.append(potential_dC,ignore_index=True)
            file2 = pd.read_csv('./data/contact_potential_prop.csv',header = None)
        #    index = ["ZHAC000102","ZHAC000103","ZHAC000105","ZHAC000102","ZHAC000103","ZHAC000105","BASU010101","BETM990101","BONM030101","BONM030102","BONM030103","BONM030104","BONM030105","BONM030106","BRYS930101","GODA950101","KESO980101","KESO980102","KOLA930101","LIWA970101","MICC010101","MIRL960101","MIYS850102","MIYS850103","MIYS960101","MIYS960102","MIYS960103","MIYS990106","MIYS990107","MOOG990101","PARB960101","PARB960102","ROBB790102","SIMK990101","SIMK990102","SIMK990103","SIMK990104","SIMK990105","SKOJ000101","SKOJ000102","SKOJ970101","TANS760101","TANS760102","THOP960101","TOBD000101","TOBD000102","VENM980101","ZHAC000101","ZHAC000104","ZHAC000106","BASU010101","BETM990101","BONM030101","BONM030102","BONM030103","BONM030104","BONM030105","BONM030106","BRYS930101","GODA950101","KESO980101","KESO980102","KOLA930101","LIWA970101","MICC010101","MIRL960101","MIYS850102","MIYS850103","MIYS960101","MIYS960102","MIYS960103","MIYS990106","MIYS990107","MOOG990101","PARB960101","PARB960102","ROBB790102","SIMK990101","SIMK990102","SIMK990103","SIMK990104","SIMK990105","SKOJ000101","SKOJ000102","SKOJ970101","TANS760101","TANS760102","THOP960101","TOBD000101","TOBD000102","VENM980101","ZHAC000101","ZHAC000104","ZHAC000106"]
            out =  pd.concat([potential_dN, potential_dC,potential_dN1, potential_dC1],ignore_index = True)
            df_out = out.to_frame()
            con = pd.concat([file2,df_out],axis = 1)
            return con.transpose()
            ####################### AA 94 properties######################
        def AAIndex_94(mutation):
            file1 = pd.read_csv('./data/94_AA_Index.csv')
            mut1 = mutation[0]+mutation[-1]
            out = file1[file1['mutation'].str.contains(mut1, na = False)]
            return out
        ############## print     
        print("####SEQUENCE BASED FEATURES####")
        print("Enter your choice.....")
        print("1. 49 amino acid Physicochemical properties")
        print("2. amino acid index properties")
        print("3. contact potential")
        print("4. neighboring residues of mutation position")
        print("5. for all")
    ######make directory out_put
        path = "./out_put"
        directory = os.path.dirname(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
    #### Input from user 
        input_user = input()
    #    print("you enterd..",input_user)
        if(input_user ==1):
            print("calculating...49 amino acid Physicochemical properties")
            prop = prop_49(seq , mutation)
            prop.to_csv("./out_put/out_prop_49.csv",header = False ,index = False,index_label = False)
        elif(input_user == 2):
            print("calculating...amino acid index properties")
            AA_Index = AAIndex_94(mutation)
            AA_Index.to_csv("./out_put/out_AA_Index.csv")
        elif(input_user ==3):
            print("calculating...contact potential")
            contact_p = contact_potential(seq, mutation)
            contact_p.to_csv("./out_put/out_contact_potential.csv",header = False ,index = False,index_label = False)
            
        elif(input_user == 4):
            print("finding...neighboring residues")
            left_pre = left_pref(seq , int(mutation[1:-1]))
            left_pre.to_csv("./out_put/neighbour.csv",index_label = True)
        elif(input_user ==5):
            print("calculating...all properties")
            contact_p = contact_potential(seq, mutation)
            AA_Index = AAIndex_94(mutation)
            prop = prop_49(seq , mutation)
            left_pre = left_pref(seq , int(mutation[1:-1]))
            contact_p.to_csv("./out_put/out_contact_potential.csv",header = False ,index = False,index_label = False)
            prop.to_csv("./out_put/out_prop_49.csv",header = False ,index = False,index_label = False)
            AA_Index.to_csv("./out_put/out_AA_Index.csv")
            left_pre.to_csv("./out_put/neighbour.csv",index_label = True)
            cmd = "paste -d, "+"./out_put/*.csv"+ ">out.csv"
            os.system(cmd)
        else:
            print("!!!Wrong choice!!!")
    
        
        
        #    out.to_csv(outFile+'_AA_Index.csv') 
        #        
        #contact_p = contact_potential(seq, mutation)
        #AA_Index = AAIndex_94(mutation)
        #prop = prop_49(seq , mutation)
        #left_pre = left_pref(seq , int(mutation[1:-1]))
def prog1(inFile,mutation,outFile):
#    inFile = sys.argv[1]
#    mutation = sys.argv[2]
    #outFile = sys.argv[3]
    
    
    
    seq = ''
    if(inFile.strip().split(".")[1]=="csv"):
        with open(inFile,'r') as j:
            d = j.readlines()
            for line1 in d:
                l_list = line1.strip("\n").split(",")
                mut_with_pos = l_list[0].strip()
                seq = l_list[1]
    elif(inFile.strip().split(".")[1]=="fasta"):
        with open(inFile,'r') as i:
            r = i.readlines()
            flag = False
            for line in r:
                if(flag):
                    seq += line.strip()        
                if ">" in line:
                    flag = True
    else:
        with open(inFile,'r') as i:
            r = i.readlines()
            for line in r:
                seq = seq +line.strip()
    # 
    #with open(inFile,'r') as i:
    #    r = i.readlines()
    #    for line in r:
    #        seq = seq +line.strip()
    #print seq
    ##print mutation
    ##
    #seq = 'MEEPQSDPSVEPPLSQETFSDLWKLLPENNVLSPLPSQAMDDLMLSPDDIEQWFTEDPGPDEAPRMPEAAPPVAPAPAAPTPAAPAPAPSWPLSSSVPSQKTYQGSYGFRLGFLHSGTAKSVTCTYSPALNKMFCQLAKTCPVQLWVDSTPPPGTRVRAMAIYKQSQHMTEVVRRCPHHERCSDSDGLAPPQHLIRVEGNLRVEYLDDRNTFRHSVVVPYEPPEVGSDCTTIHYNYMCNSSCMGGMNRRPILTIITLEDSSGNLLGRNSFEVRVCACPGRDRRTEEENLRKKGEPHHELPPGSTKRALPNNTSSSPQPKKKPLDGEYFTLQIRGRERFEMFRELNEALELKDAQAGKEPGGSRAHSSHLKSKKGQSTSRHKKLMFKTEGPDSD'
    #mutation = 'G361R'
    
    if(len(seq) <= int(mutation[1:-1])-1):
        print("!!! Short sequence length !!!")
        print("!!! Kindly check the sequence and mutation!!!")
    elif(seq[int(mutation[1:-1])-1] != mutation[0]):
         print("At given position mutant not found...!!!")
         print("Kindly check the respective position of mutation...!!!")
    
    else:
        ################ Left right pref#################
        def left_pref(seq , pos):
            if pos <=1 or pos>= len(seq):
                left_pref = ''
                right_pref = ''
                return(left_pref,right_pref)
                out_df = pd.DataFrame({'left_neighbour':[left_pref],'right_neighbour':[right_pref]})
            else:
                left_pref = seq[pos-2]
                right_pref = seq[pos]
                out_df = pd.DataFrame({'left_neighbour':[left_pref],'right_neighbour':[right_pref]})
                return(out_df)
        #        print out_df
            
        def prop_49(seq ,mutation):
            df1 = pd.read_csv('./data/49_properties_normalizedValues.csv', sep ='\t' )
            file1 = pd.read_csv('./data/prop_49_list.csv')
            mut = mutation[-1]
            wild = mutation[0]
            wrt = (df1[mut]-df1[wild])
            out = pd.concat([file1,wrt],axis =1)
        #    out_df =  out.to_frame()
            return(out.transpose())
        #    print type(out)
        
        
        #####contact potential 47*2###############
        
        def contact_potential(seq, mut_with_pos):
            mutation = mut_with_pos[0] + mut_with_pos[-1]
            pos =int(''.join( re.findall('\d+',mut_with_pos)))
            left = seq[pos-2]
            right = seq[pos]
            dN_wild = mut_with_pos[0] + left
            dN_mutant = mut_with_pos[-1]+ left
            dC_wild = mut_with_pos[0] + right
            dC_mutant = mut_with_pos[-1]+ right
            file1 = pd.read_csv('./data/data_contact_potential_diagonal.csv')
            if dN_wild not in file1:
                dN_wild = left + mut_with_pos[0]
            if dN_mutant not in file1:
                dN_mutant = left + mut_with_pos[-1]
            if dC_wild not in file1:
                dC_wild = right + mut_with_pos[0]
            if dC_mutant not in file1:
                dC_mutant = right + mut_with_pos[-1]
            potential_dN1 = file1[dN_wild] - file1[dN_mutant]
            potential_dC1 = file1[dC_wild] - file1[dC_mutant] 
            file2 = pd.read_csv('./data/data_contact_potential_square.csv')
            potential_dN = file2[dN_wild] - file2[dN_mutant]
            potential_dC = file2[dC_wild] - file2[dC_mutant] 
        #    return(potential_dN, potential_dC,potential_dN1, potential_dC1)
        #    print potential_dN.append(potential_dC,ignore_index=True)
            file2 = pd.read_csv('./data/contact_potential_prop.csv',header = None)
        #    index = ["ZHAC000102","ZHAC000103","ZHAC000105","ZHAC000102","ZHAC000103","ZHAC000105","BASU010101","BETM990101","BONM030101","BONM030102","BONM030103","BONM030104","BONM030105","BONM030106","BRYS930101","GODA950101","KESO980101","KESO980102","KOLA930101","LIWA970101","MICC010101","MIRL960101","MIYS850102","MIYS850103","MIYS960101","MIYS960102","MIYS960103","MIYS990106","MIYS990107","MOOG990101","PARB960101","PARB960102","ROBB790102","SIMK990101","SIMK990102","SIMK990103","SIMK990104","SIMK990105","SKOJ000101","SKOJ000102","SKOJ970101","TANS760101","TANS760102","THOP960101","TOBD000101","TOBD000102","VENM980101","ZHAC000101","ZHAC000104","ZHAC000106","BASU010101","BETM990101","BONM030101","BONM030102","BONM030103","BONM030104","BONM030105","BONM030106","BRYS930101","GODA950101","KESO980101","KESO980102","KOLA930101","LIWA970101","MICC010101","MIRL960101","MIYS850102","MIYS850103","MIYS960101","MIYS960102","MIYS960103","MIYS990106","MIYS990107","MOOG990101","PARB960101","PARB960102","ROBB790102","SIMK990101","SIMK990102","SIMK990103","SIMK990104","SIMK990105","SKOJ000101","SKOJ000102","SKOJ970101","TANS760101","TANS760102","THOP960101","TOBD000101","TOBD000102","VENM980101","ZHAC000101","ZHAC000104","ZHAC000106"]
            out =  pd.concat([potential_dN, potential_dC,potential_dN1, potential_dC1],ignore_index = True)
            df_out = out.to_frame()
            con = pd.concat([file2,df_out],axis = 1)
            return con.transpose()
            ####################### AA 94 properties######################
        def AAIndex_94(mutation):
            file1 = pd.read_csv('./data/94_AA_Index.csv')
            mut1 = mutation[0]+mutation[-1]
            out = file1[file1['mutation'].str.contains(mut1, na = False)]
            return out
        ############## print     
        print("####SEQUENCE BASED FEATURES####")
        print("Enter your choice.....")
        print("1. 49 amino acid Physicochemical properties")
        print("2. amino acid index properties")
        print("3. contact potential")
        print("4. neighboring residues of mutation position")
        print("5. for all")
    ######make directory out_put
        path = "./out_put"
        directory = os.path.dirname(path)
        if not os.path.exists(directory):
            os.makedirs(directory)
    #### Input from user 
        input_user = input()
    #    print("you enterd..",input_user)
        if(input_user ==1):
            print("calculating...49 amino acid Physicochemical properties")
            prop = prop_49(seq , mutation)
            prop.to_csv("./out_put/out_prop_49.csv",header = False ,index = False,index_label = False)
        elif(input_user == 2):
            print("calculating...amino acid index properties")
            AA_Index = AAIndex_94(mutation)
            AA_Index.to_csv("./out_put/out_AA_Index.csv")
        elif(input_user ==3):
            print("calculating...contact potential")
            contact_p = contact_potential(seq, mutation)
            contact_p.to_csv("./out_put/out_contact_potential.csv",header = False ,index = False,index_label = False)
            
        elif(input_user == 4):
            print("finding...neighboring residues")
            left_pre = left_pref(seq , int(mutation[1:-1]))
            left_pre.to_csv("./out_put/neighbour.csv",index_label = True)
        elif(input_user ==5):
            print("calculating...all properties")
            contact_p = contact_potential(seq, mutation)
            AA_Index = AAIndex_94(mutation)
            prop = prop_49(seq , mutation)
            left_pre = left_pref(seq , int(mutation[1:-1]))
            contact_p.to_csv("./out_put/out_contact_potential.csv",header = False ,index = False,index_label = False)
            prop.to_csv("./out_put/out_prop_49.csv",header = False ,index = False,index_label = False)
            AA_Index.to_csv("./out_put/out_AA_Index.csv")
            left_pre.to_csv("./out_put/neighbour.csv",index_label = True)
            cmd = "paste -d, "+"./out_put/*.csv"+ ">out.csv"
            os.system(cmd)
        else:
            print("!!!Wrong choice!!!")
    
        
        
        #    out.to_csv(outFile+'_AA_Index.csv') 
        #        
        #contact_p = contact_potential(seq, mutation)
        #AA_Index = AAIndex_94(mutation)
        #prop = prop_49(seq , mutation)
        #left_pre = left_pref(seq , int(mutation[1:-1]))
if(len(sys.argv)<2):
    print("Please check input aurguments...!!!")
    print("It should be look like <python sequence_based_features.py> <input_file_name> <mutation_with_position>")
elif(len(sys.argv)==3):
        prog(sys.argv[1],sys.argv[2])
elif(len(sys.argv)==4):
    prog1(sys.argv[1],sys.argv[2],sys.argv[3])
else:
    print"!!!Something is wrong with input!!!"

    
    
    
    
    
   