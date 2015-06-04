# -*- coding: utf-8 -*-
"""
Created on Thu Jun 04 17:18:17 2015

@author: Tiago
"""

def text_with_cut (self,init,lista):#distancia entre nos inicial e final
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