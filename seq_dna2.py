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
f1 = open(args.input_file,'r').readlines()[1:3]
df1 = pd.DataFrame()
df2 = pd.DataFrame()
df3 = pd.DataFrame()
df_all = pd.DataFrame()
df11 = pd.DataFrame()
df22 = pd.DataFrame()
df33 = pd.DataFrame()
df_all_1 = pd.DataFrame()


    
prop = []
if args.property == 'PHY':
    prop = ['P']
elif args.property == 'CON':
    prop = ['C']
elif args.property == 'NA':
    prop = ['L']
elif args.property == 'ALL':
    prop = ['P','C','L']
#    prop = ['A']
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
        GC_content = ((float(seq.count('G'))+ float(seq.count('C')))*100/len(seq))
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
        return var1,var2
#        print var1,var2
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
    
    #conformational properties calculation
    def conformational(self):
        f2 = pd.read_csv("./data/dna_conformational.csv")
        cols= ['Properties','Scaleunit',self.mut]
        out = f2[cols]
        return out.transpose()
        
    def physico(self):
        f3 = pd.read_csv("./data/dna_physicochemical.csv")
        cols= ['Properties','Scaleunit',self.mut]
        out = f3[cols]
        return out.transpose()
df_n1 = pd.DataFrame()
df_n2 = pd.DataFrame()
df_n3 = pd.DataFrame()




for li in f1:
    seq = li.strip().split('\t')[0]
    seq = seq.upper()
#    print(seq)
    position = ''.join(li.strip().split('\t')[1])
#    mutaa = mutaa1.split(",")



#    position = form.getvalue('position')
    alt1 = ['A','T','G','C']
#    prop = form.getlist('properties')
    window = 3
#    print position, window
#    print  seq
    for position in range(1,len(seq)):
        
        for alt in alt1:
            ref = seq[position]
#            print alt+str(position)+ref
            position = position
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
#            for i1 in range(0,len(var2)-1):
#                lis_pair_alt.append(var2[i1]+var2[i1+1])
            
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
#                    if i == 'L':
#                        lb = letter_based()
#                        f1 = lb.all_content()
#                        df3 = df3.append(f1)
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
                        
                        
#            for itm in lis_pair_alt:
#                for i in prop:
#                    if i=='P':
#                        p1 = feature_cal(itm)
#                        d1 = p1.physico()
#    #    out1 = pd.concat([df1,d1],axis = 1,ignore_index = True)
#                        df11 = df11.append(d1)
#                    if i =='C':
#                        c1 = feature_cal(itm)
#                        e1 = c1.conformational()
#                        df22 = df22.append(e1)
#                    if i == 'A':
#                        p1 = feature_cal(itm)
#                        d1 = p1.physico()
#                        df1 = df1.append(d1)
#                        c1 = feature_cal(itm)
#                        e1 = c1.conformational()
#                        df2 = df2.append(e1)
#                        lb = letter_based()
#                        f1 = lb.all_content()
#                        df3 = df3.append(f1)
    for i in prop:
        if i == 'L':
            lb = letter_based()
            f1 = lb.all_content()
            df3 = df3.append(f1)
    for j in prop:
        if j == 'P':
    #            print("<h5>Physicochemical Properties:</h5>")             
    #            df1.drop_duplicates(keep = 'first').transpose().to_csv("../../html/SBFE/out_dna/" + str(number)+"temp1.csv", index = False)
    #            print"<p>%s</p>"% ((df1.drop_duplicates(keep = 'first')).transpose().to_html(index = False))
    #            print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+str(number)+"""temp1.csv';">""")
#            print("Physicochemical Properties:")
        
            df_111 = df1.transpose()
            df_111['sum'] = df_111[lis_pair_ref].astype(float).sum(axis = 1)
            cols = ['Properties','Scaleunit','sum']
            df = df_111.transpose().drop_duplicates().transpose()[cols]
            
            df_222 = df11.transpose()
#            df_222['sum'] = df_222[lis_pair_alt].astype(float).sum(axis = 1)
#            cols = ['Properties','Scaleunit','sum']
#            df_9 = df_222.transpose().drop_duplicates().transpose()[cols]        
#    #        print"<p>%s</p>"% (df11.drop_duplicates().to_html(index = False))
#            df['Value1'] = df['sum']-df_9['sum']
            df['Average value'] = df_111['sum']/len(seq)
            cols = ['Properties','Scaleunit','Average value']
    #        print("<p>%s</p>")%(df_111.transpose().drop_duplicates().transpose().to_html())          
#            df[cols].drop_duplicates(keep = 'first').T.to_csv("./out_dna/Avg_Physicochemical.csv", index = False,header = False)
            df_n1 = df_n1.append(df[cols].T.drop_duplicates(keep = 'first'))             
#            print("%s")% (df1[cols].drop_duplicates().to_html(index = False))
    #        print"<p>%s</p>"% ((df1.drop_duplicates(keep = 'first')).transpose().to_html(index = False))
    #            print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+str(number)+"""temp1.csv';">""")
     
        if j == 'C':
#                print("Conformational Properties:")
                df_111 = df2.transpose()
                df_111['sum'] = df_111[lis_pair_ref].astype(float).sum(axis = 1)
                cols = ['Properties','Scaleunit','sum']
        #        print"<p>%s</p>"% (df2.to_html())
                df = df_111.transpose().drop_duplicates().transpose()[cols]        
        #        print"<p>%s</p>"% (df2.to_html())
#                df_222 = df22.transpose()
#                df_222['sum'] = df_222[lis_pair_alt].astype(float).sum(axis = 1)
#                cols = ['Properties','Scaleunit','sum']
#                df_9 = df_222.transpose().drop_duplicates().transpose()[cols]         
#                df['Value1'] = df['sum']-df_9['sum']
                df['Average value']  = df_111['sum']/len(seq)
                cols = ['Properties','Scaleunit','Average value']
                df_n2 = df_n2.append(df[cols].drop_duplicates(keep = 'first').T)                
#                print"<p>%s</p>"% (dfzzzzzzzzzzzzzzzz2.drop_duplicates(keep = 'first').to_html(index = False))
    #                print"<p>%s</p>"% (df2[cols].drop_duplicates().to_html(index = False))
    #                print("""<input type="button" value="Download Now!" onclick="window.location = './out_dna/temp2.csv';">""")
        if j =='L':
#            print("Nucleotide Content:")
#            df3.drop_duplicates(keep = 'first').to_csv("./out_dna/avg_Nucleotide.csv",index = False)
            new = df3
            
#            new.columns = ['% value'] 
            
#            new = new.append(new) 
#            df_n3 = df_n3.append(new)
            df_n3 = df3
    #            print"<p>%s</p>"% (new.to_html())
    #            print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp3.csv';">""")
        if j == 'A':
            print("Physico-chemical and Conformational properties:")
            df_all = df1.transpose().append(df2.transpose())
            df = df_all.loc[:,~df_all.columns.duplicated()]
#            df.to_csv("./out_dna/temp_all.csv",index = False)
    #            print("<p>%s</p>"%(df.to_html(index = False)))
    #            print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp_all.csv';">""")
            
    #            print("<p>%s</p>"%(df_all.to_html(index = False)))
#            print("Letter Based properties")
            new = df3.drop_duplicates(keep = 'first').transpose()
#            new.columns = ['% value'] 
    #            print"<p>%s</p>"% (new.to_html(header = False))
    #            df3.drop_duplicates(keep = 'first').transpose().to_csv("../../html/SBFE/out_dna/"+ str(number)+"temp3.csv", header = False)
    #            print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp3.csv';">""")
                    
                
for j in prop:
        if j == 'P':
            print("Physicochemical Properties calculated check out_dna directory")
            
            df_n1.drop_duplicates().to_csv("./out_dna/avg_Physicochemical.csv", index = True,header = False)
#            n1 = n1.append(df1[cols].T)
    #        print"<p>%s</p>"% (df1[cols].drop_duplicates().to_html(index = False))
    #        print"<p>%s</p>"% ((df1.drop_duplicates(keep = 'first')).transpose().to_html(index = False))
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+str(number)+"""temp1.csv';">""")
        elif j == 'C':
            print("Conformational Properties calculated check out_dna directory")
            
    #        
            df_n2.drop_duplicates().to_csv("./out_dna/avg_Conformational.csv", index = True,header = False)
    #        print"<p>%s</p>"% (df2.drop_duplicates(keep = 'first').to_html(index = False))
    #        print"<p>%s</p>"% (df2[cols].drop_duplicates().to_html(index = False))
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp2.csv';">""")
        elif j =='L':
            print("Nucleotide Content calculated check out_dna directory")
            df_n3.to_csv("./out_dna/avg_Nucleotide.csv",index = False)
    #        new = df3.drop_duplicates(keep = 'first').transpose()[1:]
    #        new.columns = ['Value']
    #        print"<p>%s</p>"% (new.to_html())
    #        print("""<input type="button" value="Download Now!" onclick="window.location = '../../SBFE/out_dna/"""+ str(number)+"""temp3.csv';">""")
        elif j == 'A':
            print("Physico-chemical and Conformational properties calculated check out_dna directory")
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
            df.T.to_csv("./out_dna/all.csv",index = False)
        else:
            print("somrthing is wrong")