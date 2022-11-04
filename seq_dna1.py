#!/usr/bin/python2.7
# Import modules for CGI handling 
"""
Created on Mon Apr 16 09:38:25 2018

@author: rahul
"""
import re
import pandas as pd
#import sys
import argparse
#import glob 
#import string
#import urllib
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--input_file",help = "Input file")
    parser.add_argument("--property", help = """Enter 'PHY' for Physicochemical properties,'CON' for Conformational properties, 'NA' for Nucleotide content and 'ALL' to calculate all properties""")
#    parser.add_argument("--SUB", help = "Substitution matrices")
#    parser.add_argument("--CON", help = "Pairwise properties and contact potential")
#    parser.add_argument("--ALL", help = "ALL properties")
    args = parser.parse_args()
# Create instance of FieldStorage 

#input from form
f1 = open(args.input_file,'r').readlines()[1:]

#seq1 = form.getfirst("input_seq1","0")
#mut1= ''
#mut2= ''
#mut3= ''
#mut4 = ''
#prop1 = form.getElementById('cb1')
df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()
df_all = pd.DataFrame()
df11 = pd.DataFrame()
df22 = pd.DataFrame()
df33 = pd.DataFrame()
df_all_1 = pd.DataFrame()
n1 = pd.DataFrame()
n2 = pd.DataFrame()
n3 = pd.DataFrame()
prop = []
if args.property == 'PHY':
    prop = ['P']
elif args.property == 'CON':
    prop = ['C']
elif args.property == 'NA':
    prop = ['L']
elif args.property == 'ALL':
    prop = ['P','C','L']
#print(prop)
class letter_based:
    def __inti__(self,seq):
        self.seq = seq
    def gc_content(self):
        GC_content = ((float(seq.count('G'))+ float(seq.count('C')))*100/len(seq))
        print(GC_content)
    def purine_AG(self):
        Purine_AG_content = ((float(seq.count('G'))+ float(seq.count('A')))*100/len(seq))
        print(Purine_AG_content)
    def Keto_GT(self):
        Keto_GT_content = ((float(seq.count('G'))+ float(seq.count('T')))*100/len(seq))
        print(Keto_GT_content)
    def A_content(self):
        Adenine_content = float(seq.count('A'))*100/len(seq)
        print(Adenine_content)
    def Guanine_content(self):
        Guanine_content = float(seq.count('G'))*100/len(seq)
        print(Guanine_content)
    def Cytosine_content(self):
        Cytosine_content = float(seq.count('C'))*100/len(seq)
        print(Cytosine_content)
    def Thymine_content(self):
        Thymine_content = float(seq.count('T'))*100/len(seq)
        print(Thymine_content)
    def all_content(self):
#        GC_content = ((float(seq.count('G'))+ float(seq.count('C')))*100/len(seq))
#        Purine_AG_content = ((float(seq.count('G'))+ float(seq.count('A')))*100/len(seq))
#        Keto_GT_content = ((float(seq.count('G'))+ float(seq.count('T')))*100/len(seq))
#        Adenine_content = float(seq.count('A'))*100/len(seq)
#        Guanine_content = float(seq.count('G'))*100/len(seq)
#        Cytosine_content = float(seq.count('C'))*100/len(seq)
#        Thymine_content = float(seq.count('T'))*100/len(seq)
        d = {'GC_content':[((float(seq.count('G'))+ float(seq.count('C')))*100/len(seq))],'Purine_AG_content' : [((float(seq.count('G'))+ float(seq.count('A')))*100/len(seq))],'Keto_GT_content' : [((float(seq.count('G'))+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : [float(seq.count('A'))*100/len(seq)],'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
        df1 = pd.DataFrame(data = d)

        return df1        
    
class seq_windowbased:
    # init method or constructor 
    def __init__(self, seq, position,ref,alt,window):
        self.seq = seq
        self.position = int(position)
        self.ref =ref
        self.alt = alt
        self.window = int(window)
    def win3(self):
        var1 = seq[self.position-2:self.position-1] +ref+ seq[self.position:self.position+1]
        var2 = seq[self.position-2:self.position-1] +alt+ seq[self.position:self.position+1]
#                print var1,var2 
#                mut1 = var1[0]+var1[1]
#                mut2 = var1[1]+var1[2]
#                mut3 = var2[0]+var2[1]
#                mut4 = var2[1]+var2[2]
#                return mut1,mut2,mut3,mut4,var1,var2
        return var1,var2
#                
    def win5(self):
        var1 = seq[self.position-3:self.position-1] +ref+ seq[self.position:self.position+2]
        var2 = seq[self.position-3:self.position-1] +alt+ seq[self.position:self.position+2]
        return var1,var2
    def win7(self):
        var1 = seq[self.position-4:self.position-1] +ref+ seq[self.position:self.position+3]
        var2 = seq[self.position-4:self.position-1] +alt+ seq[self.position:self.position+3]
        return var1,var2
    def win9(self):
        var1 = seq[self.position-5:self.position-1] +ref+ seq[self.position:self.position+4]
        var2 = seq[self.position-5:self.position-1] +alt+ seq[self.position:self.position+4]
        return var1,var2
    def win11(self):
        var1 = seq[self.position-6:self.position-1] +ref+ seq[self.position:self.position+5]
        var2 = seq[self.position-6:self.position-1] +alt+ seq[self.position:self.position+5]
        return var1,var2

class feature_cal:
    def __init__(self,mut):
        self.mut = mut
#                self.mut1 = mut1
#                self.mut2 = mut2
#                self.mut3 = mut3
#                self.mut4 = mut4
#            
    #conformational properties calculation
    def conformational(self):
        f1 = pd.read_csv("./data/dna_conformational.csv")
        cols= ['Properties','Scaleunit',self.mut]
        out = f1[cols]
        
#        outdf = out.to_frame()
        return out.transpose()
        
    def physico(self):
        f1 = pd.read_csv("./data/dna_physicochemical.csv")
        cols= ['Properties','Scaleunit',self.mut]
        out = f1[cols]
#                df_temp = pd.DataFrame()
#                df_temp['value'] = f1[self.mut1]+f1[self.mut2]-f1[self.mut3]-f1[self.mut4]
#                print(df_temp.to_html())
#                out = pd.concat([f1[cols],df_temp],axis = 1)
#                print out
        return out.transpose()
#                return out
            
for li in f1:
    seq = li.strip().split('\t')[0]
#    print(seq)
    mutaa1 = ''.join(li.strip().split('\t')[1])
    mutaa = mutaa1.split(",")
    #    ref = muta[0]
    #    alt = muta[-1]
    #    position = int(''.join(re.findall('\d+',muta)))
    #    prop = form.getlist('properties')
    #    window = form.getvalue('window')
#    prop = form.getlist('properties')
    window = 3
    for muta in mutaa:
        
    #    muta = form.getvalue('position')
        ref = muta[0].upper()
        alt = muta[-1].upper()
        position = int(''.join(re.findall('\d+',muta)))
#        print(muta)
        
        
        #parse fasta file 
        #var =  seq1.splitlines()
        #flag = False
        
        #for lin in var:
        #    if(flag):
        #        se += lin.strip(' \n\t\r')
        #    if lin.startswith(">"):
        #        name = lin
        #        flag = True
        #flag = False
        ######
        
#        if(len(prop)<1):
#            print("Kindly select at least one property...!!!</p>")
##        elif(len(seq)<1):
#            print("<p>Kindly enter or upload valid input sequence...!!!</p>")
        if(True):
        
            
            p = seq_windowbased(seq,position,ref,alt,window)
            if(int(window) ==3):
                var1,var2 = p.win3()
            elif(int(window) ==5):
                var1,var2 = p.win5()
            elif(int(window)== 7):
                var1,var2 = p.win7()
            elif(int(window) ==9):
                var1,var2 = p.win9()
            elif(int(window) ==11):
                var1,var2 = p.win11()
            else:
                print("<p>WRONG WINDOW SIZE</p>")
            #print("""<p>var1:%s,var2:%s,window:%s</p>""")%(var1,var2,window) 
            lis_pair_ref = []
            lis_pair_alt = []
            for i in range(0,len(var1)-1):
                lis_pair_ref.append(var1[i]+var1[i+1])
            for i1 in range(0,len(var2)-1):
                lis_pair_alt.append(var2[i1]+var2[i1+1])
            
        fA = pd.DataFrame()
        fT = pd.DataFrame()
        fG = pd.DataFrame()
        fC = pd.DataFrame()
        fA1 = pd.DataFrame()
        fT1 = pd.DataFrame()
        fG1 = pd.DataFrame()
        fC1 = pd.DataFrame()
#        
        for item in lis_pair_ref:
            for i in prop:
                if i=='P':
                    p1 = feature_cal(item)
                    d1 = p1.physico()
#                    df1= df1.append(d1) 
        #    out1 = pd.concat([df1,d1],axis = 1,ignore_index = True)
                    df1 = df1.append(d1)
                if i =='C':
                    c1 = feature_cal(item)
                    e1 = c1.conformational()
                    df2 = df2.append(e1)
                if i == 'L':
#                    lb = letter_based()
                    if alt =="A" and ref =="T":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT = pd.DataFrame(data = d)
                        df3 = fT-fA
                    if alt =="A" and ref =="G":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG = pd.DataFrame(data = d)
                        df3 = fG-fA
                    if alt =="A" and ref =="C":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC = pd.DataFrame(data = d)
                        df3 = fC-fA
                    if alt =="T"and ref =="A":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT = pd.DataFrame(data = d)
                        df3 = fA-fT
                    if alt =="T"and ref =="G":
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT = pd.DataFrame(data = d)
                        df3 = fG-fT
                    if alt =="T"and ref =="C":
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT = pd.DataFrame(data = d)
                        df3 = fC-fT
                    if alt =="C" and ref =="A":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC = pd.DataFrame(data = d)
                        df3 = fA-fC
                    if alt =="C" and ref =="G":
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC = pd.DataFrame(data = d)
                        df3 = fG-fC
                    if alt =="C" and ref =="T":
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC = pd.DataFrame(data = d)
                        df3 = fT-fC
                    if alt =="G" and ref =="A":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG = pd.DataFrame(data = d)
                        df3 = fA-fG
                    if alt =="G" and ref =="T":
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG = pd.DataFrame(data = d)
                        df3 = fT-fG
                    if alt =="G" and ref =="C":
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG = pd.DataFrame(data = d)
                        df3 = fT-fG
                    
#                    if ref =="A":
#                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
#                        fA = pd.DataFrame(data = d)
##                        df3 = df3.append(f1)
#                    if ref =="T":
#                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
#                        fT = pd.DataFrame(data = d)
#                    
#                    if ref =="C":
#                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
#                        fC = pd.DataFrame(data = d)
#                    if alt =="C":
#                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
#                        fC = pd.DataFrame(data = d)
#                    if alt =="G":
#                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
#                        fG = pd.DataFrame(data = d)
#                    if ref =="G":
#                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
#                        fG = pd.DataFrame(data = d)
                    
                if i == 'A':
                    p1 = feature_cal(item)
                    d1 = p1.physico()
                    df1 = df1.append(d1)
                    c1 = feature_cal(item)
                    e1 = c1.conformational()
                    df2 = df2.append(e1)
                    lb = letter_based()
                    f1 = lb.all_content()
                    df3 = df3.append(f1)
#        print"<p>%s</p>"%(df1.to_html())
#        
                    
#        df1['value']= df1.transpose()[lis_pair_ref].sum()
#        print"<p>%s</p>"%(df1.transpose().to_html())
                    
        for itm in lis_pair_alt:
            for i in prop:
                if i=='P':
                    p1 = feature_cal(itm)
                    d1 = p1.physico()
        #    out1 = pd.concat([df1,d1],axis = 1,ignore_index = True)
                    df11 = df11.append(d1)
                    
                if i =='C':
                    c1 = feature_cal(itm)
                    e1 = c1.conformational()
                    df22 = df22.append(e1)
                if i == 'L':
                    lb = letter_based()
                    if alt =="A" and ref =="T":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA1 = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT1 = pd.DataFrame(data = d)
                        df33 = fT1-fA1
                    if alt =="A" and ref =="G":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG1 = pd.DataFrame(data = d)
                        df33 = fG1-fA1
                    if alt =="A" and ref =="C":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC1 = pd.DataFrame(data = d)
                        df33 = fC1-fA1
                    if alt =="T"and ref =="A":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA1 = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT1 = pd.DataFrame(data = d)
                        df33 = fA1-fT1
                    if alt =="T"and ref =="G":
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG1 = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT1 = pd.DataFrame(data = d)
                        df33 = fG1-fT1
                    if alt =="T"and ref =="C":
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC1 = pd.DataFrame(data = d)
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT1 = pd.DataFrame(data = d)
                        df33 = fC1-fT1
                    if alt =="C" and ref =="A":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC1 = pd.DataFrame(data = d)
                        df33 = fA1-fC1
                    if alt =="C" and ref =="G":
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC1 = pd.DataFrame(data = d)
                        df33 = fG1-fC1
                    if alt =="C" and ref =="T":
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC1 = pd.DataFrame(data = d)
                        df33 = fT1-fC1
                    if alt =="G" and ref =="A":
                        d = {'GC_content':0,'Pyrimidine_TC_content':0,'Purine_AG_content':[float(seq.count('A'))*100/len(seq)],'Keto_GT_content':0,'Adenine_content':[float(seq.count('A'))*100/len(seq)],'Guanine_content':0,'Cytosine_content':0,'Thymine_content' : 0}
                        fA1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG1 = pd.DataFrame(data = d)
                        df33 = fA1-fG1
                    if alt =="G" and ref =="T":
                        d = {'GC_content':0,'Pyrimidine_TC_content':((float(0)+ float(seq.count('T')))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : [((float(0)+ float(seq.count('T')))*100/len(seq))],'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : 0,'Thymine_content' : [float(seq.count('T'))*100/len(seq)]}
                        fT1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG1 = pd.DataFrame(data = d)
                        df33 = fT1-fG1
                    if alt =="G" and ref =="C":
                        d = {'GC_content':[((float(0)+ float(seq.count('C')))*100/len(seq))],'Pyrimidine_TC_content':((float(seq.count('C'))+ float(0))*100/len(seq)),'Purine_AG_content' : 0,'Keto_GT_content' : 0,'Adenine_content' : 0,'Guanine_content' : 0,'Cytosine_content' : [float(seq.count('C'))*100/len(seq)],'Thymine_content' : 0}
                        fC1 = pd.DataFrame(data = d)
                        d = {'GC_content':[((float(seq.count('G'))+ 0)*100/len(seq))],'Pyrimidine_TC_content':0,'Purine_AG_content' : [((float(seq.count('G'))))*100/len(seq)],'Keto_GT_content' : [((float(seq.count('G'))+ 0))*100/len(seq)],'Adenine_content' : 0,'Guanine_content' : [float(seq.count('G'))*100/len(seq)],'Cytosine_content' : 0,'Thymine_content' : 0}
                        fG1 = pd.DataFrame(data = d)
                        df33 = fT1-fG1
#                    f1 = lb.all_content()
#                    df33 = df33.append(f1)
                if i == 'A':
                    p1 = feature_cal(itm)
                    d1 = p1.physico()
                    df11 = df11.append(d1)
                    c1 = feature_cal(itm)
                    e1 = c1.conformational()
                    df22 = df22.append(e1)
                    lb = letter_based()
                    f1 = lb.all_content()
                    df33 = df33.append(f1)
        
    for j in prop:
        if j == 'P':
            print("Physicochemical Properties:")
            
            df_111 = df1.transpose()
    #        print(df1)
            df_111['sum'] = df_111[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df = df_111.transpose().drop_duplicates().transpose()[cols]
            
            df_222 = df11.transpose()
            df_222['sum'] = df_222[lis_pair_alt].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df_9 = df_222.transpose().drop_duplicates().transpose()[cols]        
    #        print"<p>%s</p>"% (df11.drop_duplicates().to_html(index = False))
            
            
            df['Value'] = df['sum']-df_9['sum']
            cols = ['Properties','Scaleunit','Value']
    #        print("<p>%s</p>")%(df_111.transpose().drop_duplicates().transpose().to_html())          
#            df1[cols].T.to_csv("./out_dna/Physicochemical.csv", index = True,header = False)
            n1 = n1.append(df[cols].T)
    #        print"<p>%s</p>"% (df1[cols].drop_duplicates().to_html(index = False))
    #        print"<p>%s</p>"% ((df1.drop_duplicates(keep = 'first')).transpose().to_html(index = False))
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+str(number)+"""temp1.csv';">""")
        if j == 'C':
            print("Conformational Properties:")
            df_111 = df2.transpose()
            df_111['sum'] = df_111[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
    #        print"<p>%s</p>"% (df2.to_html())
            df = df_111.transpose().drop_duplicates().transpose()[cols]        
    #        print"<p>%s</p>"% (df2.to_html())
            df_222 = df22.transpose()
            df_222['sum'] = df_222[lis_pair_alt].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df9 = df_222.transpose().drop_duplicates().transpose()[cols]         
            df['Value'] = df['sum']-df9['sum']
            cols = ['Properties','Scaleunit','Value']
    #        
#            df2[cols].T.to_csv("./out_dna/Conformational.csv", index = True,header = False)
            n2 = n2.append(df[cols].T)
    #        print"<p>%s</p>"% (df2.drop_duplicates(keep = 'first').to_html(index = False))
    #        print"<p>%s</p>"% (df2[cols].drop_duplicates().to_html(index = False))
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp2.csv';">""")
        if j =='L':
            print("Nucleotide Content:")
            df3.drop_duplicates(keep = 'first').to_csv("./out_dna/Nucleotide.csv",index = False)
            n3 = n3.append(df3.drop_duplicates(keep = 'first'))
    #        new = df3.drop_duplicates(keep = 'first').transpose()[1:]
    #        new.columns = ['Value']
    #        print"<p>%s</p>"% (new.to_html())
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp3.csv';">""")
        if j == 'A':
            print("Physico-chemical and Conformational properties:")
            df_1112 = df2.transpose()
            df_1112['sum'] = df_1112[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df2 = df_1112.transpose().drop_duplicates().transpose()[cols]           
            df_111 = df1.transpose()
            df_111['sum'] = df_111[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df1 = df_111.transpose().drop_duplicates().transpose()[cols]        
            
            df_all = df1.append(df2)
            df = df_all.loc[:,~df_all.columns.duplicated()]
            df.to_csv("./out_dna/all.csv",index = False)
        
#        print("<p>%s</p>"%(df.to_html(index = False)))
#        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp_all.csv';">""")
#        
#                print("<p>%s</p>"%(df_all.to_html(index = False)))
#        print("<h5>Letter Based properties</h5>")
#        print(df3.drop_duplicates(keep = 'first').transpose())
#        print"<p>%s</p>"% (df3.drop_duplicates(keep = 'first').transpose().to_html(columns = cols))
#        df3.drop_duplicates(keep = 'first').transpose().to_csv("../../html/SBFE/out_dna/"+ str(number)+"temp3.csv", header = False)
#        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp3.csv';">""")
for j in prop:
        if j == 'P':
            print("Physicochemical Properties:")
            
            n1.drop_duplicates().to_csv("./out_dna/Physicochemical.csv", index = True,header = False)
#            n1 = n1.append(df1[cols].T)
    #        print"<p>%s</p>"% (df1[cols].drop_duplicates().to_html(index = False))
    #        print"<p>%s</p>"% ((df1.drop_duplicates(keep = 'first')).transpose().to_html(index = False))
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+str(number)+"""temp1.csv';">""")
        if j == 'C':
            print("Conformational Properties:")
            
    #        
            n2.drop_duplicates().to_csv("./out_dna/Conformational.csv", index = True,header = False)
    #        print"<p>%s</p>"% (df2.drop_duplicates(keep = 'first').to_html(index = False))
    #        print"<p>%s</p>"% (df2[cols].drop_duplicates().to_html(index = False))
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp2.csv';">""")
        if j =='L':
            print("Nucleotide Content:")
            n3.to_csv("./out_dna/Nucleotide.csv",index = False)
    #        new = df3.drop_duplicates(keep = 'first').transpose()[1:]
    #        new.columns = ['Value']
    #        print"<p>%s</p>"% (new.to_html())
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp3.csv';">""")
        if j == 'A':
            print("Physico-chemical and Conformational properties:")
            df_1112 = df2.transpose()
            df_1112['sum'] = df_1112[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df2 = df_1112.transpose().drop_duplicates().transpose()[cols]           
            df_111 = df1.transpose()
            df_111['sum'] = df_111[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df1 = df_111.transpose().drop_duplicates().transpose()[cols]        
            
            df_all = df1.append(df2)
            df = df_all.loc[:,~df_all.columns.duplicated()]
            df.to_csv("./out_dna/all.csv",index = False)
