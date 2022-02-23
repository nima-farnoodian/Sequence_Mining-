"""
AUTHOR(S): 
Nima Farnoodian <nima.farnoodian@uclouvain.be>, 
EPL, UCLouvain, Belgium 
"""

import sys
from collections import defaultdict
from itertools import combinations
import time

class Dataset_Sequence:
    """Utility class to manage a dataset stored in a external file."""

    def __init__(self, filepath):
        """reads the dataset file and initializes files"""
        """This class reads the file for sequence mining. Each line is an item 
        with its location and the transactions are separated by '\n'"""
        self.transactions = list()
        self.items = set()
        self.Freq_Items={}
        #self.minSup=minSup
        self.Vertical_Rep={}
        self.V2H=[]
        try:
            f=open(filepath,'r')
            line=[]
            idx=0
            Freq_Items={}
            visited_items=set()
            horizental=[]
            for itemLoc in f.readlines():
                if itemLoc!='\n':
                    itemLoc=itemLoc.strip()
                    item=itemLoc.split(' ')[0]
                    loc=int(itemLoc.split(' ')[1])
                    line.append(item)
                    if tuple([item]) not in self.Vertical_Rep:
                        self.Vertical_Rep[tuple([item])]={}
                    if idx not in self.Vertical_Rep[tuple([item])]:
                        self.Vertical_Rep[tuple([item])][idx]=[]
                    
                    self.Vertical_Rep[tuple([item])][idx].append(loc)
                    #self.Vertical_Rep[frozenset([item])].append((idx,loc))
                    visited_items.add(item)
                    horizental.append((item,loc))
      
                else:
                    idx+=1
                    if idx>1:
                        for item in visited_items:
                            if item not in self.Freq_Items:
                                self.Freq_Items[item]=1
                                self.items.add(frozenset([item]))
                            else:
                                self.Freq_Items[item]+=1
                        visited_items=set()
                        self.transactions.append(line)
                        self.V2H.append(horizental)
                        horizental=[]
                    line=[]
            #self.transactions.append(line)
        except IOError as e:
            print("Unable to read dataset file!\n" + e)

    def trans_num(self):
        """Returns the number of transactions in the dataset"""
        return len(self.transactions)

    def items_num(self):
        """Returns the number of different items in the dataset"""
        return len(self.items)

    def get_transaction(self, i):
        """Returns the transaction at index i as an int array"""
        return self.transactions[i]
    def get_items(self):
        return self.items
    def get_freq_items(self,minSup):
        TrNo=len(self.transactions)
        AboveFreq=[i for i in self.Freq_Items if (self.Freq_Items[i]/TrNo)>=minSup]
        return AboveFreq
    
    def get_freq_Vertical_Rep(self,minSup):
        #  it returns the vertical representation of the frequent items
        if minSup!=None:
            TrNo=self.items_num()
            High_freq=self.get_freq_items(minSup)
            Vertical_Rep_high={}
            for it in High_freq:
                Vertical_Rep_high[tuple([it])]=self.Vertical_Rep[tuple([it])]
        return Vertical_Rep_high
    def get_items_in_order(self):
        sequences=[i for i in sorted(self.Freq_Items.items(), key=lambda item: item[1],reverse=True)]
        return sequences

def vertical_finding(vr,projected_vr,right_item):
    """This function is used to find the sequence using vertical representation"""
    right_item=tuple(right_item)
    left=projected_vr
    right=vr.get(right_item,dict())
    res={}
    support=0
    for x in left:
        if x in right:
            break_Flage=False
            for loc_left in left[x]:
                if break_Flage==False:
                    for loc_right in right[x]:
                        if loc_left<loc_right:
                            if x not in res:
                                res[x]=[]
                            res[x].append(loc_right)
                            support+=1
                            break_Flage=True
                            break
                else:
                    break
    return res,support


def f2_spade(data, minFrequency,printing=True):
    """f2_SPADE Algorithm- Recursive implementation"""
                    
    vr=data.get_freq_Vertical_Rep(minFrequency) 
    total_trans = data.trans_num() 
    SequenceSet={}
    items=data.get_freq_items(minFrequency)
    explored=items.copy()
    for idx in range(len(items)):
        item=items[idx]
        ExploredNew=explored.copy()
        #ExploredNew=ExploredNew[idx:]
        supp=(len(vr[tuple([item])]))
        if supp/total_trans>=minFrequency:
            SequenceSet[tuple([item])]=supp
            
    #  Computing frequent 2 sequence
    f2_to_extend=set()
    t_id=1
    for transaction in data.V2H:
        sequences = combinations(transaction, 2)
        visited={}
        transaction=1
        for seq in sequences:
            item1=seq[0][0]
            item1_loc=seq[0][1]
            item2=seq[1][0]
            item2_loc=seq[1][1]
            f2seq=(item1,item2)
            #print(seq)
            if f2seq not in vr:
                vr[f2seq]={}
            if f2seq not in visited:
                visited[f2seq]=item1_loc
                SequenceSet[f2seq]=SequenceSet.get(f2seq,0)+1
                vr[f2seq][t_id]=[item2_loc]
                if SequenceSet[f2seq]/total_trans>=minFrequency:
                    f2_to_extend.add(f2seq)
            else:
                '''
                if item1_loc>visited[f2seq]:
                    if item2_loc>vr[f2seq][t_id][-1]:
                        vr[f2seq][t_id].append(item2_loc)
                '''
        t_id+=1
        
    for f2_item in f2_to_extend:
        depthFirstSearch(list(f2_item), minFrequency, total_trans,vr,SequenceSet,ExploredNew,vr[f2_item])
    if printing==True:
        for i in SequenceSet:
            print(str(list(i))+' '+"("+str(SequenceSet[i])+')')
    return SequenceSet

def depthFirstSearch(item,minFrequency, total_trans,vr,SequenceSet,explored,projected_vr):
    exploredLocal=explored.copy()
    #print(exploredLocal)
    if type(item)!=list:
        item=[item]
    #while (len(exploredLocal)>0):
    for item2 in exploredLocal:
        #item2=exploredLocal.pop()
        if True:
        #if item2 not in item:
            item2=[item2]
            res,support=vertical_finding(vr,projected_vr,item2)
            supp=support/total_trans
 
            if(supp > minFrequency):
                sequence=item+item2
                #print(tuple(sequence))
                if tuple(sequence) not in SequenceSet:
                    SequenceSet[tuple(sequence)]=support
                    #projected_vr[tuple(sequence)]=res
                    #print(itemsorted, '('+str(supp)+')')
                    depthFirstSearch(sequence,minFrequency, total_trans,vr,SequenceSet,explored,res)
    return None
        
#!/usr/bin/env python3

import sys


