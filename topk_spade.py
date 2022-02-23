# Author: Nima Farnoodian, EPL, UCLouvain.
import numpy as np
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
    
class k_selector:
    """Utility class to control k-top elements."""



    def __init__(self, k):
        self.score_list = list() # a list of top k (or less than k) scores that have been observed so far
        self.score_link = dict() # a dictionary where keys are the scores and the values are the sequences whose scores are equale to the key
        self.k=k # the number of unique scores
    def append(self,score,sequence):
    
        '''
        Input: 
            score e.g, total support or Wracc
            obtained sequence
        Output:
            return True if the sequence could be added to the top-k list sequences
        '''
        ret=False
        if score not in self.score_link:
            if len(self.score_list)<self.k:
                self.score_list.append(score)
                self.score_link[score]=set()
                self.score_link[score].add(sequence)
                ret=True
            else:
                minimum=np.min(self.score_list)
                if score>minimum:
                    self.score_list.pop(self.score_list.index(minimum))
                    del(self.score_link[minimum])
                    self.score_link[score]=set()
                    self.score_list.append(score)
                    self.score_link[score].add(sequence)
                    ret=True
                if score==minimum:
                    self.score_link[score].add(sequence)
                    ret=True
        else:
            self.score_link[score].add(sequence)
            ret=True
        return ret
    
def vertical_finding(vr,left_vr,right_item):
    """This function is used to find the sequence using vertical representation (Sparse Lattice suggested in Zaki's paper)"""

    '''
    Input: 
        vertical representation provided by Dataset_Sequence class
        left: vertical representation of the prefix e.g, AB
        right_item= the candidate item added to the prefix e.g, C

    Output:
        res= vertical representation of the prefix + right_item  e.g, ABC
        support= the support of the found sequence
    '''
    right_item=tuple(right_item)
    right=vr.get(right_item,dict())
    res={}
    support=0
    for x in left_vr:
        if x in right:
            break_Flage=False
            for loc_left in left_vr[x]:
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


def topk_spade(data_positive,data_negative, k,printing=True):
    """ 
    Top k SPADE Algorithm- Recursive DFS implementation
    k is the number of top sequences that should be seleceted.
    Notice that, both positive and negative supports for sequences are found in the same call using vertical finding function, thereby speeding up the algorithm.
    
    """

    '''
    Input: 
        data_positive= positive class data provided by Dataset_Sequence class
        data_negative= negative class data provided by Dataset_Sequence class
        k: k is the number of top sequences that should be seleceted.
        printing= if True, print the top k frequent sequences (e.g, [A, B, C] supp_positive supp_negative score)

    Output:
        None 
    '''
    if k==0:
        print('k==0 is not feasible')
        return None, None
    selector=k_selector(k)
    vr_positive=data_positive.get_freq_Vertical_Rep(0) 
    vr_negative=data_negative.get_freq_Vertical_Rep(0) 
    
    total_trans_positive = data_positive.trans_num() 
    total_trans_negative = data_negative.trans_num() 
    total_trans=total_trans_positive+total_trans_negative
    SequenceSet={}
    items=list(set(data_positive.get_freq_items(0)).union(set(data_negative.get_freq_items(0)))) # find all items observed in positive and negative classes
    explored=items.copy()
    ordered_item={}
    # the following block select only top k single items. It is used to prune the search as top k problem benefits from Anti-monotonicity 
    for item in items:
        supp_positive=(len(vr_positive.get(tuple([item]),[])))
        supp_negative=(len(vr_negative.get(tuple([item]),[])))
        total_supp=supp_positive+supp_negative
        selector.append(total_supp,tuple([item]))
    selected_item=[]
    scores=sorted(selector.score_link.keys(),reverse=True)
    for score in scores:
            for sequence in selector.score_link[score]:
                selected_item.append(sequence)
    selected_item=[i[0] for i in selected_item]
    # end of block 
    for item in selected_item:
        supp_positive=(len(vr_positive.get(tuple([item]),[])))
        supp_negative=(len(vr_negative.get(tuple([item]),[])))
        total_supp=supp_positive+supp_negative
        if total_supp>0:
            SequenceSet[tuple([item])]=[supp_positive,supp_negative,total_supp]
            project_positive=vr_positive.get(tuple([item]),dict())
            project_negative=vr_negative.get(tuple([item]),dict())
            topk_depthFirstSearch(item, total_trans,vr_positive,vr_negative,SequenceSet,selected_item,project_positive,project_negative,selector)
    if printing==True:
        counter=0
        scores=sorted(selector.score_link.keys(),reverse=True)
        for score in scores:
            for sequence in selector.score_link[score]:
                out='['+', '.join(list(sequence))+']' 
                print(out,SequenceSet[sequence][0],SequenceSet[sequence][1],SequenceSet[sequence][2] )
                counter+=1
                #print(counter)
    return selector.score_link

def topk_depthFirstSearch(item, total_trans,vr_positive,vr_negative,SequenceSet,exploredLocal,projected_vr_positive,projected_vr_negative,selector):
    if type(item)!=list:
        item=[item]
    for item2 in exploredLocal:
        if True:
            item2=[item2]
            total_supp=0
            res_positive,supp_positive=vertical_finding(vr_positive,projected_vr_positive,item2)
            res_negative,supp_negative=vertical_finding(vr_negative,projected_vr_negative,item2)
            total_supp=supp_positive+supp_negative
            if(total_supp/total_trans > 0.0):
                sequence=item+item2
                selected=selector.append(total_supp,tuple(sequence))
                if selected: # if the sequence was added to the top k list, so it is frequent and its supersequences may be top-k frequent sequences as well
                    if tuple(sequence) not in SequenceSet:
                        SequenceSet[tuple(sequence)]=[supp_positive,supp_negative,total_supp]
                        topk_depthFirstSearch(sequence, total_trans,vr_positive,vr_negative,SequenceSet,exploredLocal,res_positive,res_negative,selector) # discover supersequences
    return SequenceSet
import sys

'''
def main():
    pos_filepath = sys.argv[1] # filepath to positive class file
    neg_filepath = sys.argv[2] # filepath to negative class file
    k = int(sys.argv[3])
    # TODO: read the dataset files and call your miner to print the top k itemsets
    ds_positive=Dataset_Sequence(pos_filepath)
    ds_negative=Dataset_Sequence(neg_filepath)
    
    selected=topk_spade(ds_positive, ds_negative,k,printing=True)


if __name__ == "__main__":
    main()
'''
