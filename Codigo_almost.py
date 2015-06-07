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
                print (str(k)+"->"+str(self.nodes[k][1]))
            else:
                print (str(k)+":"+str(self.nodes[k][0]))

###########################################################
###Construcao da arvore de compactacao###
###########################################################

#adicionar um no a arvore
    
    def addNode(self,origin,symbol,leafnum=-1):
        self.num+=1
        self.nodes[origin][1][symbol]=self.num
        self.nodes[self.num]=(leafnum,{})

#adicionar um sufixo (palavra)
         
    def addSufix(self,p,sufnum):
        pos=0
        node=0
        while pos<len(p):
            if not (p[pos]) in self.nodes[node][1]:
                if pos==len(p)-1:
                    self.addNode(node,p[pos],sufnum)
                else:
                    self.addNode(node,p[pos])
            node=self.nodes[node][1][p[pos]]
            pos+=1


#apartir de uma sequencia criar uma arvore de sufixos

    def suffixTrieFromSeq(self,text):
        t=text+"$"
        for i in range (len(t)):
            self.addSufix(t[i:],i)

            
#para saber se e biforcado ou nao, dado o no respetivo
    def biforcado(self, node):
        try:
            if len(self.nodes[node][1].keys())>1:#muitas chaves, logo biforcado
                return True
            elif len(self.nodes[node][1].keys())==0:#folha
                return -1
            else:#nao biforcado
                return False
        except:
            return False
        

#funcao que retira o texto, dados os nos iniciais e o caminho para atingir o no final do ramo em questao (em forma de lista)
    
    def text_with_cut (self,init,lista):
        text=""
        keys=self.nodes[init][1].keys()#todas as chaves do primeiro no
        for key in keys:#por cada uma:
            next_node=self.nodes[init][1][key]#no seguinte
            if next_node==lista[0]:#verifica a existencia do no seguinte no caminho dado pela lista
                text+=key
                for i in lista:#confirmada a situacao anterior percorre-se a lista para retirar o texto e eliminar os nos
                    letter=self.nodes[i][1].keys()
                    letter2=list(letter)
                    if len(letter2) >=1 and not self.biforcado(i):#nao pode ser 0
                        text+=letter2[0]
                    else:
                        pass
                for i in range(len(lista)-1):#todos menos o ultimo elemento 
                    del self.nodes[lista[i]]#apaga os nos do dicionario
        return text
            
#nodos sao eliminados na funcao anterior mas temos de ligar o no inicial ao no final
            
    def encurtar(self,nod):
        dicio={}
        nos=self.getLastNodeBelow_main(nod)#vai buscar os caminhos em forma de lista (no biforcado: lista de listas porque um no biforcado tem varios caminhos)
        for stop in range(len(nos)):
            text=self.text_with_cut(nod,nos[stop])#passa a funcao o no inicial e uma lista do caminho a seguir
            dicio[text]=nos[stop][-1]#cria uma entrada no dicionario com o texto e o ultimo no como valor
        self.nodes[nod][1].clear()#limpa o dicionario do no na arvore
        for key in dicio.keys():#acrescentar as chaves do dicionario temporario criado ao dicionario do no na arvore
            self.nodes[nod][1][key]=dicio[key]
        
       
#cria uma lista do caminho, dado um no nao biforcado
       
    def getLastNodeBelow(self,node):#so aceita nos nao biforcados iniciais
        res=[]
        if self.nodes[node][0]>=0:
            res.append(node)#ultimo no (correspondente a uma folha)
        elif not self.biforcado(node):#percorre o caminho ate encontrar uma folha ou um no biforcado
            res.append(node)            
            newnode=self.nodes[node][1].values()
            list_newnode=list(newnode)
            stop=self.getLastNodeBelow(list_newnode[0])
            res.extend(stop)         
        else:#se for biforcado
            res.append(node)
        return res
        

#cria a lista do caminho mas, sendo no biforcado retorna uma matriz de listas, sendo cada uma correspondente as opcoes do no

    def getLastNodeBelow_main(self,no):#biforcado ou nao
        res=[]
        keys=self.nodes[no][1].keys()
        if len(keys)==1:#se nao for biforcado
            final_node=self.getLastNodeBelow(no)
            res.append(final_node)
        elif len(keys)>1:#biforcado
            for key in keys:
                next_node=self.nodes[no][1][key]
                final_node_bifor=self.getLastNodeBelow(next_node)
                res.append(final_node_bifor)
        return res


#transforma uma arvore de sufixos numa arvore de compactacao
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

#encontar um padrao na arvore

    def findPattern(self,pattern):
        pos=0
        node=0
        while pos<len(pattern)-1:
            for key in self.nodes[node][1].keys():
                lenght=len(key)
                if pattern[pos:pos+lenght] in self.nodes[node][1]:# a pesquisa e feita pelo tamanho de cada chave encontrada
                    node=self.nodes[node][1][key]
                    pos+=(lenght-1)#como encontrou a posicao passa a ser a pos inicial mais o tamanho da key
                else:#nao corresponde por isso passa pra chave seguinte
                    pass
        if pos==len(pattern)-1:#depois de terminado o ciclo anterior, verifica-se onde a posicao acabou
            posic=self.getLeafesBelow(node)#chama a funcao para saber o valor da folha do ultimo no
            print ('O padrao foi encontrado na sua totalidade na sequencia na posicao %s da sequencia' %posic)        
        elif pos>=10:#considera-se que se encontra de forma parcial quando mais de 10 nucleotidos sao iguais
            posic=self.getLeafesBelow(node)
            for i in posic:#pode ser encontrado em varios sitios            
                tam=int(i)+int(pos)
                print ("O padrao foi encontado, de forma parcial, da posicao %s ate a posicao %s da sequencia" %posic,tam)                   
        return posic
    
#retorna a folha ou folhas abaixo de um dado no
    
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
###Retirar as sequencias das bases de dados (padrao-gene) e sufixos (genoma)###
###############################################################################        

#dado um gi ou um nome da especie a funcao retorna a sequencia desse dado input, em formato de texto
      
    def get_seq_from_genome(self,especie=0,idn=0):
        if idn==0:
            gf=self.get_genome_file(especie)
        elif especie==0:
            gf=self.get_genome_file(0,idn)
        f=open(gf, 'r')#a funcao chamada retorna o nome do ficheiro criado, chamado aqui
        sequence = SeqIO.read(f, 'fasta')
        print ("A sequencia usada para construir a arvore sera:\n"+str(sequence.id))#para informar qual sera a sequencia usada
        return str(sequence.seq).strip(" ")
        
#dado um gi ou um nome da especia a funcao vai extrair a sequencia da base dados no ncbi e guarda-la num ficheiro, retornando o seu nome
        
    def get_genome_file(self,especie=0,idn=0):
        return_filename=""        
        if idn==0:#damos a especie
            hand=Entrez.esearch(db='nucleotide',term=especie+"[ORGN]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            handle=Entrez.efetch(db='nucleotide',rettype="fasta",retmode="text",id=idnum)
            read=SeqIO.read(handle,"fasta")
            name="genome_"+str(especie).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
        elif especie==0:#damos o id a funcao
            handle=Entrez.efetch(db='nucleotide',rettype="fasta",retmode="text",id=idn)
            read=SeqIO.read(handle,"fasta")
            name="genome_"+str(idn).strip(" ")+".fasta"
            SeqIO.write(read,name,"fasta")
            handle.close()
            return_filename+=name
        return return_filename

#dado um gi ou um nome do gene a funcao retorna a sequencia desse dado input, em formato de texto

    def get_seq_from_gene(self,gene=0,idn=0):
        if idn==0:
            gf=self.get_gene_file(gene)
        elif gene==0:
            gf=self.get_gene_file(0,idn)
            print (gf)
        f=open(gf, 'r')
        sequence = SeqIO.read(f, 'fasta')
        print ("A sequencia do gene pesquisado é:\n"+str(sequence.id))#para informar qual sera a sequencia usada
        return str(sequence.seq).strip(" ")
       
   
    def get_gene_file(self,gene=0,idn=0):#cria o ficheiro e da o nome a funcao anterior para retirar a sequencia do ficheiro
        return_filename=""        
        if idn==0:#damos o nome da especie a funcao 
            hand=Entrez.esearch(db='gene',term=gene+"[sym]",retmax=100,retype="gb",retmode="text")
            results=Entrez.read(hand)
            idnum=results["IdList"][0]#primeiro elemento da lista de resultados
            handl=Entrez.efetch(db='gene',rettype="gb",retmode="xml",id=idnum)
            record=Entrez.read(handl)
            to=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_to']#para

            desde=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_from']#de

            identif=record[0]['Entrezgene_comments'][5]['Gene-commentary_comment'][0]['Gene-commentary_comment'][0]['Gene-commentary_seqs'][0]['Seq-loc_int']['Seq-interval']['Seq-interval_id']['Seq-id']['Seq-id_gi']#gi
                        
            handle = Entrez.efetch(db="nucleotide", rettype="fasta", retmode="text", id=identif, seq_start=desde, seq_stop=to )#depois de termos o id, e os valores onde comeca e acaba a seq vamos busca-la
            read=SeqIO.read(handle,"fasta")
            name="gene_"+str(gene).strip(" ")+".fasta"
            SeqIO.write(read,name, "fasta")
            handle.close()
            return_filename+=name
        elif gene==0:#damos o gi do gene
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
        print("1 - Contrucao de arvore de sufixos apartir de uma sequencia ja conhecida")
        print("2 - Construcao da arvore de sufixos apartir de um ficheiro FASTA")
        print("3 - Extracao de uma sequencia da base de dados")
        print("4 - Pesquisa de uma sequencia no genoma")
        print("5 - Imprimir a arvore criada") 
        print("6 - Sair")
        op=input("O que pretende fazer? ")
        if op==1:
            seq=input('\nInsira a sequencia para a construcao da arvore:\n')
            st.suffixTrieFromSeq(str(seq))#dada a seq vamos construir a arvore
            st.s_to_c_tree()#e comprimi-la
            print('\Arvore criada com sucesso\n')
        elif op=="2":
            opt=input('\nPretende a pesquisa de um genoma na base de dados?(S ou N):\n').upper()
            if opt=='S':
                  opt1=input('\nPossui o GI?(S ou N):\n').upper()
                  if opt1=='S':#procura pelo GI
                      idn=input('\n Introduza o GI:\n').upper()
                      seq=st.get_seq_from_genome(0,str(idn))
                      st.suffixTrieFromSeq(seq)
                      st.s_to_c_tree()
                      print('\nArvore criada com sucesso\n')
                     
                  if opt1=='N':#procura pela especie
                      especie=input('\n Introduza o nome da especie:\n')
                      seq=st.get_seq_from_genome(especie)
                      st.suffixTrieFromSeq(seq)
                      st.s_to_c_tree()
                      print('\nArvore criada com sucesso\n')

            if opt=='N':#entao e necessario dar o ficheiro fasta da sequencia
                file_fasta=input('\nIntroduza o nome do ficheiro FASTA (sem extensao):\n')
                seq=st.get_seq_from_genome(str(file_fasta)+'.fasta')
                st.suffixTrieFromSeq(seq)
                st.s_to_c_tree()
                print('\nArvore criada com sucesso\n')
            
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
                    especie=input('\n Introduza o nome da especie:\n')
                    st.get_seq_from_genome(especie)
                    print ('\nGenoma extraido com sucesso para a diretoria de trabalho\n')
                    
        elif op=="4":
            opt=input('\nPretende a pesquisa de um gene na base de dados?(S ou N):\n').upper()
            if opt=='S':
                  opt1=input('\nPossui o GI?(S ou N):\n').upper()
                  if opt1=='S':
                      idn=input('\n Introduza o GI:\n').upper()
                      seq=st.get_seq_from_gene(0,str(idn))
                      st.findPattern(seq)
                     
                  if opt1=='N':
                      gene=input('\n Introduza o nome do gene:\n')
                      seq=st.get_seq_from_gene(gene)
                      st.suffixTrieFromSeq(seq)
                      st.findPattern(seq)

            if opt=='N':
                inp_seq=input('\nIntroduza a sequencia a pesquisar?\n')
                st.findPattern(str(inp_seq))
              
        elif op=="5":
            st.print_tree()
            
        elif op=="6":
            print("\n Adeus!")
            break
        else:
            print("\n Opcao nao valida! Tente novamente.")
            

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
        seq="TCAA"
        st=CreateTree()
        st.suffixTrieFromSeq(seq)
        #print(st.biforcado(2))##
        #print st.text_with_cut(1,13)##
        #print st.print_tree()##corta bem os ficheiros 
        #print (st.encurtar(0))
        #print (st.print_tree())
        #print type('drt')
        #print st.getLastNodeBelow()
        #print st.getLastNodeBelow_main(1)
        st.s_to_c_tree()
        #(st.print_tree())
        st.findPattern('CAA$')
        #print st.getLeafesBelow(1)
        #print st.get_genome_file('escherichia coli')
        
        
        
    #main()
    test3()