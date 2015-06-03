# -*- coding: utf-8 -*-
"""
Created on Tue May 26 20:02:54 2015

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
            
######################################################################################################3
    def nobiforca(self,node=0):#procura por caminhos nao bifurcados dado um nodo.
        if len(self.nodes[node][1].keys())==0:#if len(self.nodes[node][1].keys())==0:#atingimos uma folha
            res=int(self.nodes[node][0])#numero da folha colocado na lista
            print res
            return res
            #cada no retornara apenas uma folha no caso de ser um caminho nao biforcado
        #for k in self.nodes[0][1].keys():
        elif len(self.nodes[node][1].keys())==1:#tem um no destino
            new_node=self.nodes[node][1].values()
            new_node2=new_node[0]#vai buscar o elemento a lista dos valores (neste caso é apenas um, por isso retira-se o primeiro)
            #print (new_node2)
            self.nobiforca(new_node2)
        else:
            return False#encontra uma bifurcação
        #return res
        
        
    def InitialSearch(self):#método para procurar em todos os caminhos provenientes da raiz
        res=[]
        for k in self.nodes[0][1].keys():#raiz possui muitos nodos ligados
            #print k
            #print self.nodes[0][1][k]
            search=self.nobiforca(self.nodes[0][1][k])
            #print search
            if search is not False:
                print search
                res.append(search)
        return res
                
#ok, até aqui temos as folhas onde existem caminhos diretos da raiz ate elas... e agora? 
                
            
    def encurta(self):
        root_to_leafs=self.InitialSearch()
        word=""
        for node in self.nodes.keys():
            for leaf in root_to_leafs:
                if self.nodes[node][0]>=0 and self.nodes[node][0]==leaf:#corresponde a uma folha, sendo aceites apenas as folhas sem biforcações
                    self.nodes             
            
            
            ###########################
    def getPredecessors(self,leafnode):
        word=""
        num_node=-1 

        while num_node!=0:
            for node in self.nodes.keys():
                for letter in self.nodes[node][1]:
                    destin_node=self.nodes[node][1][letter]
    
                    if destin_node==leafnode:
                        num_node=self.nodes[node]
                        word=str(self.nodes[node][0])
                        #insere a letra na ultima posicao da lista (o texto vai do fim para o inicio)
        print word
          ####################################################################   ############################
                
                
#para saber se é biforcado ou nao, dado o no respetivo
                
    def biforcado(self, node):
        #print len(self.nodes[node][1].keys())        
        try:
            if len(self.nodes[node][1].keys())>1:
                return True
            elif len(self.nodes[node][1].keys())==0:
                return -1
            else:
                return False
        except:
            print "errado"
            return False###########
        
#funcao que retirao texto, dados os nos iniciais e finais do ramo em questao
        
    def text_with_cut (self,init,final):#distancia entre nos inicial e final
        text=""
        #sdicio={}
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
 #acho que isto nao sera preciso uma vez que a seguinte apaga tudo mas confirmar  #del self.nodes[init][1][key]#apaga a chave correspondente ao caminho do no biforcado em estudo
                    #self.nodes[init][1][]=
                
        #if len(self.nodes[final][1].keys())==0:
        elif len(keys)==1:
            text+=keys[0]
            next_node=self.nodes[init][1][keys[0]]
     #########       del self.nodes[init]#apaga o no inicial
            while next_node!=final:
                letter2=self.nodes[next_node][1].keys()[0]
                text+=letter2
                destin=self.nodes[next_node][1].values()
                destin_node=destin[0]
                del self.nodes[next_node]            
                next_node=destin_node
         #if len(self.nodes[final][1].keys())==0:
         #   text+="$"
        #texto=text.strip(",")#em principio nao sera necessario uma vez que foi usada a concatenacao de strings
        #if len(texto)>=2:
        #if self.text_with_cut(init,final))=="s:
        return text
        #else:
            #return self.nodes[init][1].keys[0]
            
        
            
            
#nodos sao eliminados na funcao anterior mas temos de ligar o no inicial ao no final            
    def encurtar(self,nod):
        dicio={}        
        #if self.biforcado(node):
        try:
            nos=self.getLastNodeBelow_main(nod)
            #print nos
            for stop in range(len(nos)):
                
                print nos[stop]
                print nod
                text=self.text_with_cut(nod,nos[stop])
                #print nod                
                #print nos[stop]
                #print text
                #if text not False:                
                dicio[text]=nos[stop]
                #else:
                    #dicio[node][1][]
                #eliminar o dicionario existente    
            self.nodes[nod][1].clear()
            
            for key in dicio.keys():
                #print key
                self.nodes[nod][1][key]=dicio[key]
                #print self.nodes[node][1][key]
        #atribuicao do dicionario ao no inicial correspondente
        #else:
         #   text=self.text_with_cut(init,final)        
          #  dicio[text]=final
           # self.nodes[init][1]=dicio
        except:
            return False
        
            
    def getLastNodeBelow(self,node):#so aceita nos nao biforcados iniciais
        res=[]
        if self.nodes[node][0]>=0:
           # print node
            res.append(node)#ultimo no (correspondente a uma folha)
     #   elif len(self.nodes[node][1].keys())>1:
            
        elif not self.biforcado(node):
   #for k in self.nodes[node][1].keys():
            newnode=self.nodes[node][1].values()#lista
          #  print newnode
            stop=self.getLastNodeBelow(newnode[0])
            res.extend(stop)         
        else:#se for biforcado
            res.append(node)
        #print res
        return res
        

#imprime ou um ou varios nos finais
    
    def getLastNodeBelow_main(self,no):#biforcado ou nao
        res=[]
        keys=self.nodes[no][1].keys()
        print keys        
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
    
    
   # def direct_line_to_leaf(self,node):
        
   #     res=[]        
   #     if len(self.nodes[node][1].keys())==0:
    #        #atingiu-se uma folha
   #        res.append(node)
   #     elif len(self.nodes[node][1].keys())>1:
            
        
  #  def valida tamnho 2
        
        
        
    
    def s_to_c_tree(self):
        num=0
        while num is not self.num:
            print num
            bifor=self.biforcado(num)
            print bifor            
            try:
                if bifor==-1:
                    pass
                else:
                    self.encurtar(num)
#                if bifor:
 #                   for key in self.nodes[num][1].keys():
  #                      node=self.nodes[num][1][key]
   #                     print node
    #                    self.encurtar(node)
                #elif not bifor:
#                    self.encurtar(num)
 #               else:
  #                  pass
            except:
                pass
            num+=1
        ########
            
    #ate agora faz o pretendido mas falta compor funcoes para ficarem a juntar os caracteres das bifurcacoes. 
            #uma vez que so juntam o texto apartir do no seguinte de cada bifurcacao, 
            #deixando um no que podia ser eliminado para tras 
            #(fica o seguinte a bifurcacao e o final e podia ficar logo que se da a bifurcacao ligado ao final)
#para isso tem de se adaptar a text_with_cut ou criar outra que chame essa e junte o caracter inicial ao texto 
            #que sai dessa- ter atencao que se tem de eliminar o nodo)
        #a adaptar a funcao getlastnodebelow, uma vez que tem de permitir nos iniciais biforcados e dar um output de
            #quais os destinos da funcao cut
            
            #cut (initial, path,final)#ideia
########################################################
###Algoritmo de pesquisa###
########################################################





    def findPattern(self,pattern):
#        pos=0        
#        node=0
#        for k in self.nodes[0][1].keys():
#            if key[pos]==pattern[pos]:
#                
#            while not self.biforcado(node) or self.biforcado(node)==-1:
#                for key in self.nodes[node][1].keys():
#                    while pos<len(key):
#                        if key[pos]==pattern[pos]:
#                            pos+=1     
#                        else:
#                            return False
#                        
#            node=self.nodes[node][1][key]
#                       
#                
#                
#
#      # pos_suf=0
#        for nodes in self.nodes:
#            
#            for pos_pad in range(len(pattern)):
#                for key in self.nodes[node][1].keys():
#                    pos_suf=0
#                    if key[0]==pattern[0]:
#                        pos_suf+=1
#                        
#                        
#                    
#                    
#                    
#            if self.nodes[node][1].has_key(pattern[pos]):
#                node=self.nodes[node][1][pattern[pos]]
#                pos+=1
#            else:
#                return None
#        return self.getLeafesBelow(node)

        pos=0
        node=0
        while pos<len(pattern):
            for key in self.nodes[node][1].keys():
                lenght=len(key)
                    
                #for pos in range(len(pattern)):
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
      
#vamos supor que a arvore esta a funcionar direitinho... PAsso seguinte:
      
      
      
    def get_seq_from_genome(self,especie=0,idn=0):
        if idn==0:
            gf=self.get_genome_file(especie)
         #   f=open(gf+".fasta", 'r+')
        #    sequence = SeqIO.read(f, 'fasta')
       #     print str(sequence.seq).strip(" ")
        elif especie==0:
            gf=self.get_genome_file(0,idn)
            print gf
        f=open(gf, 'r')
        sequence = SeqIO.read(f, 'fasta')
        #print str(sequence.seq).strip(" ")
        print "A sequencia usada para construir a arvore sera:\n"+str(sequence.id)#para informar qual sera a sequencia usada
        return str(sequence.seq).strip(" ")
        
        

       
       
       
    def get_genome_file(self,especie=0,idn=0):#cria o ficheiro e da o nome a funcao anterior para retirar a sequencia do ficheiro
        return_filename=""        
        
        if idn==0:
            hand=Entrez.esearch(db='nucleotide',term=especie+"[ORGN]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            print idnum
            
            handle=Entrez.efetch(db='nucleotide',rettype="fasta",retmode="text",id=idnum)
     #       record=Entrez.read(handle)
            #print(record[0].keys())
            #print (record[0]["TSeq_defline"])
            #print handle.read()
            
            
            read=SeqIO.read(handle,"fasta")
            name="genome_"+str(especie).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
  ### #         record = SeqIO.read("genome_escherichia coli.fasta", "fasta")
  ###  #        print record             
           # filename = "genome_"+especie+".gb"
        
        elif especie==0:
            handle=Entrez.efetch(db='nucleotide',rettype="fasta",retmode="text",id=idn)
            read=SeqIO.read(handle,"fasta")
            name="genome_"+str(idn).strip(" ")+".fasta"
            SeqIO.write(read,name,"fasta")
            handle.close()
            return_filename+=name
        #print return_filename
        return return_filename
        

        #f=open(filename,"w+")
        #f.write(handle.read())    
        #print handle.read()
        #return record["IdList"]#especie
        #rint record
        
    
        
        
       
       
       
    def get_seq_from_gene(self,especie=0,idn=0):
        if idn==0:
            gf=self.get_genome_file(especie)
         #   f=open(gf+".fasta", 'r+')
        #    sequence = SeqIO.read(f, 'fasta')
       #     print str(sequence.seq).strip(" ")
        elif especie==0:
            gf=self.get_genome_file(0,idn)
            print gf
        f=open(gf, 'r')
        sequence = SeqIO.read(f, 'fasta')
        #print str(sequence.seq).strip(" ")
        print "A sequencia usada para construir a arvore sera:\n"+str(sequence.id)#para informar qual sera a sequencia usada
        return str(sequence.seq).strip(" ")
        
        

       
       
       
    def get_gene_file(self,gene=0,idn=0):#cria o ficheiro e da o nome a funcao anterior para retirar a sequencia do ficheiro
        return_filename=""        
        
        if idn==0:
            hand=Entrez.esearch(db='gene',term=gene+"[sym]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            print idnum
   
        
            handl=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=idnum)
            record=Entrez.read(handl)
            #print(record[0].keys())#["Entrezgene_summary"])
            to=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para

            desde=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de

            identif=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                        
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )
            #print (record[0]["TSeq_defline"])
            #print handle.read()
            
        #[u'Entrezgene_unique-keys', u'Entrezgene_track-info', u'Entrezgene_gene', u'Entrezgene_location', 
            #u'Entrezgene_xtra-properties', u'Entrezgene_locus', u'Entrezgene_homology', u'Entrezgene_source', 
         #   u'Entrezgene_type', u'Entrezgene_gene-source', u'Entrezgene_comments', u'Entrezgene_xtra-index-terms', 
          #  u'Entrezgene_properties', u'Entrezgene_prot', u'Entrezgene_summary']   
           
           
           
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(gene).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
  ### #         record = SeqIO.read("genome_escherichia coli.fasta", "fasta")
  ###  #        print record             
           # filename = "genome_"+especie+".gb"
        
        elif gene==0:
            handl=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=idn)
            record=Entrez.read(handl)
            #print(record[0].keys())#["Entrezgene_summary"])
            to=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para

            desde=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de

            identif=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                        
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )
            #print (record[0]["TSeq_defline"])
            #print handle.read()
            
        #[u'Entrezgene_unique-keys', u'Entrezgene_track-info', u'Entrezgene_gene', u'Entrezgene_location', 
            #u'Entrezgene_xtra-properties', u'Entrezgene_locus', u'Entrezgene_homology', u'Entrezgene_source', 
         #   u'Entrezgene_type', u'Entrezgene_gene-source', u'Entrezgene_comments', u'Entrezgene_xtra-index-terms', 
          #  u'Entrezgene_properties', u'Entrezgene_prot', u'Entrezgene_summary']   
           
           
           
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(idn).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
        #print return_filename
      #  return return_filename
        

        #f=open(filename,"w+")
        #f.write(handle.read())    
        #print handle.read()
        #return record["IdList"]#especie
        #rint record
        
    
        
        
       

      
      
      
      
      
      
      
from Bio import Entrez      
      
##############################################################3
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
        #print st.encurtar(1)
        #print st.print_tree()
        #print type('drt')
        #print st.getLastNodeBelow(5)
        #print st.getLastNodeBelow_main(1)
        #st.s_to_c_tree()
        #print st.print_tree()
        #print st.findPattern('TA')
        #print st.getLeafesBelow(1)
        #print st.get_genome_file('escherichia coli')
        
        
        
    test3()