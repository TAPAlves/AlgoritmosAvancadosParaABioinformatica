# -*- coding: utf-8 -*-
"""
Created on Thu Jun 04 15:54:57 2015

@author: Tiago
"""

from Bio import Entrez
from Bio import SeqIO
Entrez.email='tiago_alves26@hotmail.com'

class CreateTree:
    
    def __init__(self):
        self.nodes={0:(-1,{})}#root node
        self.num=0 #numero de no

##imprimir a arvore
        
    def print_tree(self):
        for k in self.nodes.keys():
            if self.nodes[k][0]<0:
                print str(k)+"->"+str(self.nodes[k][1])
            else:
                print str(k)+":"+str(self.nodes[k][0])

###########################################################
###Construcao da arvore de compactacao###
###########################################################
    
    def addNode(self,origin,symbol,leafnum=-1):
        self.num+=1
        self.nodes[origin][1][symbol]=self.num
        self.nodes[self.num]=(leafnum,{})
        
    def addSufix(self,p,sufnum):
        pos=0
        node=0
        while pos<len(p):
            if not self.nodes[node][1].has_key(p[pos]):
                if pos==len(p)-1:
                    self.addNode(node,p[pos],sufnum)
                else:
                    self.addNode(node,p[pos])
            node=self.nodes[node][1][p[pos]]
            pos+=1

    def suffixTrieFromSeq(self,text):
        t=text+"$"
        for i in range (len(t)):
            self.addSufix(t[i:],i)
            
#para saber se é biforcado ou nao, dado o no respetivo
                
    def biforcado(self, node):
        try:
            if len(self.nodes[node][1].keys())>1:
                return True
            elif len(self.nodes[node][1].keys())==0:
                return -1
            else:
                return False
        except:
            print ("errado")
            return False
        

#funcao que retirao texto, dados os nos iniciais e finais do ramo em questao
    
    def text_with_cut (self,init,final):#distancia entre nos inicial e final
        text=""
        keys=self.nodes[init][1].keys()
        if len(keys)>1:
            for key in keys:
                next_node=self.nodes[init][1][key]
                final_node=self.getLastNodeBelow(next_node)#vai ver se o no final é o mesmo que é passado pela funcao anterior. Se o no seguinte for logo um biforcado não ocorre união de texto uma vez que so podiamos juntar uma letra apenas
                if final_node[0]==final:
                    text+=key
                    while next_node!=final:
                        letter2=self.nodes[next_node][1].keys()[0]
                        text+=letter2
                        destin=self.nodes[next_node][1].values()
                        destin_node=destin[0]
                        del self.nodes[next_node]            
                        next_node=destin_node
        elif len(keys)==1:
            text+=keys[0]
            next_node=self.nodes[init][1][keys[0]]
            while next_node!=final:
                letter2=self.nodes[next_node][1].keys()[0]
                text+=letter2
                destin=self.nodes[next_node][1].values()
                destin_node=destin[0]
                del self.nodes[next_node]            
                next_node=destin_node
        return text
            
#nodos sao eliminados na funcao anterior mas temos de ligar o no inicial ao no final            
    def encurtar(self,nod):
        dicio={}        
        try:
            nos=self.getLastNodeBelow_main(nod)
            for stop in range(len(nos)):
                stop1=int(stop)
                text=self.text_with_cut(nod,nos[stop1])
                dicio[text]=nos[stop1]
                #eliminar o dicionario existente    
            self.nodes[nod][1].clear()
            for key in dicio.keys():
                self.nodes[nod][1][key]=dicio[key]
        except:
            return False
        
            
    def getLastNodeBelow(self,node):#so aceita nos nao biforcados iniciais
        res=[]
        if self.nodes[node][0]>=0:
            res.append(node)#ultimo no (correspondente a uma folha)
            
        elif not self.biforcado(node):
            newnode=self.nodes[node][1].values()#lista
            stop=self.getLastNodeBelow(newnode[0])
            res.extend(stop)         
        else:#se for biforcado
            res.append(node)
        return res
        

#imprime ou um ou varios nos finais
    
    def getLastNodeBelow_main(self,no):#biforcado ou nao
        res=[]
        keys=self.nodes[no][1].keys()
        if len(keys)==1:#se nao for biforcado
            final_node=self.getLastNodeBelow(no)
            res.extend(final_node)
        elif len(keys)>1:#biforcado
            for key in keys:
                next_node=self.nodes[no][1][key]
                print next_node                
                final_node_bifor=self.getLastNodeBelow(next_node)
                res.extend(final_node_bifor)
        return res

    def s_to_c_tree(self):
        num=0
        while num is not self.num:
            bifor=self.biforcado(num)
            try:
                if bifor==-1:
                    pass
                else:
                    self.encurtar(num)
            except:
                pass
            num+=1

########################################################
###Algoritmo de pesquisa###
########################################################

    def findPattern(self,pattern):
        pos=0
        node=0
        while pos<len(pattern):
            for key in self.nodes[node][1].keys():
                lenght=len(key)
                if self.nodes[node][1].has_key(pattern[pos:pos+lenght]):
                    node=self.nodes[node][1][key]
                    pos+=(lenght-1)
                else:
                    pass
        if pos!=0:
            posic=self.getLeafesBelow(node)
            print "O padrao foi encontado na posição xxxx"#######################
            return posic

            if pos==len(pattern)-1:
                posic=self.getLeafesBelow(node)
                print 'o padrao foi encontrado na sua totalidade na sequencia'                
                return posic
    
    
    def getLeafesBelow(self,node):
        res=[]
        if self.nodes[node][0]>=0:
            res.append(self.nodes[node][0])
        else:
            for k in self.nodes[node][1].keys():
                newnode=self.nodes[node][1][k]
                leafes=self.getLeafesBelow(newnode)
                res.extend(leafes)         
        return res
        

###############################################################################
###Retirar as sequencias das bases de dados (padrão-gene) e sufixos (genoma)###
###############################################################################        
      
    def get_seq_from_genome(self,especie=0,idn=0):
        if idn==0:
            gf=self.get_genome_file(especie)
        elif especie==0:
            gf=self.get_genome_file(0,idn)
            print gf
        f=open(gf, 'r')
        sequence = SeqIO.read(f, 'fasta')
        print "A sequencia usada para construir a arvore sera:\n"+str(sequence.id)#para informar qual sera a sequencia usada
        return str(sequence.seq).strip(" ")
        
       
    def get_genome_file(self,especie=0,idn=0):#cria o ficheiro e da o nome a funcao anterior para retirar a sequencia do ficheiro
        return_filename=""        
        if idn==0:
            hand=Entrez.esearch(db='nucleotide',term=especie+"[ORGN]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            handle=Entrez.efetch(db='nucleotide',rettype="fasta",retmode="text",id=idnum)
            read=SeqIO.read(handle,"fasta")
            name="genome_"+str(especie).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
        elif especie==0:
            handle=Entrez.efetch(db='nucleotide',rettype="fasta",retmode="text",id=idn)
            read=SeqIO.read(handle,"fasta")
            name="genome_"+str(idn).strip(" ")+".fasta"
            SeqIO.write(read,name,"fasta")
            handle.close()
            return_filename+=name
        return return_filename

    def get_seq_from_gene(self,especie=0,idn=0):
        if idn==0:
            gf=self.get_genome_file(especie)
        elif especie==0:
            gf=self.get_genome_file(0,idn)
            print gf
        f=open(gf, 'r')
        sequence = SeqIO.read(f, 'fasta')
        print "A sequencia usada para construir a arvore sera:\n"+str(sequence.id)#para informar qual sera a sequencia usada
        return str(sequence.seq).strip(" ")
       
    def get_gene_file(self,gene=0,idn=0):#cria o ficheiro e da o nome a funcao anterior para retirar a sequencia do ficheiro
        return_filename=""        
        
        if idn==0:
            hand=Entrez.esearch(db='gene',term=gene+"[sym]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            handl=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=idnum)
            record=Entrez.read(handl)
            to=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para

            desde=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de

            identif=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                        
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(gene).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
        elif gene==0:
            handl=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=idn)
            record=Entrez.read(handl)
            to=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para

            desde=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de

            identif=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                        
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )
           
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(idn).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name

###############################################################################
################################  M  E  N  U  #################################         
###############################################################################        
       
def main():
    st=CreateTree()    
    while True:
        print("1 - Contrução de árvore de sufixos apartir de uma sequência já conhecida")
        print("2 - Construção da árvore de sufixos apartir de um ficheiro FASTA")
        print("3 - Extração de uma sequência da base de dados")
        print("4 - Pesquisa de uma sequência no genoma")
        print("5 - Imprimir a árvore criada") 
        #print("6 - " ) 
        print("9 - Sair")
        op=input("O que pretende fazer? ")

        if op==1:
            seq=input('\nInsira a sequência para a construção da árvore:\n')
            st.suffixTrieFromSeq(str(seq))
            st.s_to_c_tree()
            print('\nÁrvore criada com sucesso\n')
        elif op=="2":
            opt=input('\nPretende a pesquisa de um genoma na base de dados?(S ou N):\n').upper()
            if opt=='S':
                  opt1=input('\nPossui o GI?(S ou N):\n').upper()
                  if opt1=='S':
                      idn=input('\n Introduza o GI:\n').upper()
                      seq=st.get_seq_from_genome(0,str(idn))
                      st.suffixTrieFromSeq(seq)
                      st.s_to_c_tree()
                      print('\nÁrvore criada com sucesso\n')
                     
                  if opt1=='N':
                      especie=input('\n Introduza o nome da espécie:\n')
                      seq=st.get_seq_from_genome(especie)
                      st.suffixTrieFromSeq(seq)
                      st.s_to_c_tree()
                      print('\nÁrvore criada com sucesso\n')

            if opt=='N':
                file_fasta=input('\nIntroduza o nome do ficheiro FASTA (sem extensão):\n')
                seq=st.get_seq_from_genome(str(file_fasta)+'.fasta')
                st.suffixTrieFromSeq(seq)
                st.s_to_c_tree()
                print('\nÁrvore criada com sucesso\n')
            
        elif op=="3":
            word=input('\nPretende extrair um gene ou um genoma?').upper()
            if word=='GENE':
                opt1=input('\nPossui o GI?(S ou N):\n').upper()
                if opt1=='S':
                    idn=input('\n Introduza o GI:\n').upper()
                    st.get_seq_from_gene(0,str(idn))
                    print ('\nGene extraido com sucesso para a diretoria de trabalho\n')
                if opt1=='N':
                    gene=input('\n Introduza o nome do gene:\n')
                    st.get_seq_from_gene(gene)
                    print ('\nGene extraido com sucesso para a diretoria de trabalho\n')
            
            if word=='GENOMA':
                
                opt1=input('\nPossui o GI?(S ou N):\n').upper()
                if opt1=='S':
                    idn=input('\n Introduza o GI:\n').upper()
                    st.get_seq_from_genoma(0,str(idn))
                    print ('\nGenoma extraido com sucesso para a diretoria de trabalho\n')
                if opt1=='N':
                    especie=input('\n Introduza o nome da espécie:\n')
                    st.get_seq_from_genome(especie)
                    print ('\nGenoma extraido com sucesso para a diretoria de trabalho\n')
                    
        elif op=="4":
            opt=input('\nPretende a pesquisa de um gene na base de dados?(S ou N):\n').upper()
            if opt=='S':
                  opt1=input('\nPossui o GI?(S ou N):\n').upper()
                  if opt1=='S':
                      idn=input('\n Introduza o GI:\n').upper()
                      seq=st.get_seq_from_gene(0,str(idn))
                      pos=st.findPattern(seq)
                      print('\nPadrão encontado na posição:\n')+pos
                     
                  if opt1=='N':
                      gene=input('\n Introduza o nome do gene:\n')
                      seq=st.get_seq_from_gene(gene)
                      st.suffixTrieFromSeq(seq)
                      pos=st.findPattern(seq)
                      print('\nPadrão encontado na posição:\n')+pos
            if opt=='N':
                inp_seq=input('\nIntroduza a sequencia a pesquisar?\n')
                pos=st.findPattern(str(inp_seq))
                print('\nPadrão encontado na posição:\n')+pos               
        elif op=="5":
            st.print_tree()
        elif op=="9":
            print("\n Adeus!")
            break
        else:
            print("\n Opção não valida! Tente novamente.")
            

if __name__=='__main__':
    
    def test():
        seq="TACTA"
        st=CreateTree()
        st.suffixTrieFromSeq(seq)
        st.print_tree()
  #  test()
    
    def test2():
        seq="TACT"
        st=CreateTree()
        st.suffixTrieFromSeq(seq)
    #    print st.findPattern("TA")
      #  print st.InitialSearch()
        #print st.getPredecessors(9)#st.print_tree()       
        #print st.text_with_cut(6,9)
        #print st.getLastNodeBelow(2)
        #print st.encurtar(2)
        #print st.s_to_c_tree()
        #print st.text_with_cut(2,5)
       # print st.biforcado(2)
        #print st.findPattern("")
       # print st.pesquisa_genoma()
        #a=Entrez.einfo()
        #print (Entrez.read(a))
        #print st.get_genome_file("escherichia coli")
        #print st.get_gene_file("idh1")#id tem de ser da bd gene
        #print st.getLastNodeBelow_main(1)
    #test2()

    def test3():
        seq="TACT"
        st=CreateTree()
        st.suffixTrieFromSeq(seq)
        #print(st.biforcado(2))##
        #print st.text_with_cut(1,13)##
        #print st.print_tree()##corta bem os ficheiros 
        print st.encurtar(1)
        #print st.print_tree()
        #print type('drt')
        #print st.getLastNodeBelow(5)
        #print st.getLastNodeBelow_main(1)
        #st.s_to_c_tree()
        #print st.print_tree()
        #print st.findPattern('TA')
        #print st.getLeafesBelow(1)
        #print st.get_genome_file('escherichia coli')
        
        
        
    #main()
    test3()