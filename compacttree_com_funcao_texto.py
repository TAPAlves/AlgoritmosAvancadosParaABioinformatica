# -*- coding: utf-8 -*-
"""
Created on Tue May 26 20:02:54 2015

@author: Tiago
"""

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
             ############################
                
                
#para saber se é biforcado ou nao, dado o no respetivo
                
    def biforcado(self, node):
        if len(self.nodes[node][1].keys())==1:
            return False
        elif len(self.nodes[node][1].keys())==0:
            return -1
        else:
            return True
            
        
#funcao que retirao texto, dados os nos iniciais e finais do ramo em questao
        
    def text_with_cut (self,init,final):#distancia entre nos inicial e final
        text=""
        letter=self.nodes[init][1].keys()[0]
        text+=letter
        destin0=self.nodes[init][1].values()
        init_cicle_node=destin0[0]
        while init_cicle_node!=final:
            letter2=self.nodes[init_cicle_node][1].keys()[0]
            text+=letter2
            destin=self.nodes[init_cicle_node][1].values()
            destin_node=destin[0]
            del self.nodes[init_cicle_node]            
            init_cicle_node=destin_node
        #if len(self.nodes[final][1].keys())==0:
         #   text+="$"
        texto=text.strip(",")#em principio nao sera necessario uma vez que foi usada a concatenacao de strings
        if len(texto)>=2:
            return texto
        else:
            return False
            
        
            
            
#nodos sao eliminados na funcao anterior mas temos de ligar o no inicial ao no final            
    def encurtar(self,node):
        dicio={}        
        #if self.biforcado(node):
        nos=self.getLastNodeBelow(node)
        for leaf in nos:
            text=self.text_with_cut(node,leaf)
            dicio[text]=leaf
        #eliminar o dicionario existente    
        self.nodes[node][1].clear()
        for key in dicio.keys():
            self.nodes[node][1][key]=dicio[key]
        #atribuicao do dicionario ao no inicial correspondente
        #else:
         #   text=self.text_with_cut(init,final)        
          #  dicio[text]=final
           # self.nodes[init][1]=dicio
        
        
            
    def getLastNodeBelow(self,node):
        res=[]
        if self.nodes[node][0]>=0:
            res.append(node)#ultimo no (correspondente a uma folha)
        else:
            for k in self.nodes[node][1].keys():
                newnode=self.nodes[node][1][k]
                leafes=self.getLastNodeBelow(newnode)
                res.extend(leafes)         
        return res
        
    
    
    

        
            
                    
#criar um metodo para encurtar as folhas ate a raiz dooss nao biforcados
#depois e so chamar aqui na funcao de cima            
            
            
            
            


########################################################
###Algoritmo de pesquisa###
########################################################





    def findPattern(self,pattern):
        pos=0
        node=0
        for pos in range(len(pattern)):
            if self.nodes[node][1].has_key(pattern[pos]):
                node=self.nodes[node][1][pattern[pos]]
                pos+=1
            else:
                return None
        return self.getLeafesBelow(node)
        
        
        
        
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
        #print st.getLastNodeBelow(6)
        print st.encurtar(6)
        
    test2()