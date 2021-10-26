#!/usr/bin/python
# ecole polytechnique - c.durr - 2009

# Kuhn-Munkres, The hungarian algorithm.  Complexity O(n^3)
# Computes a max weight perfect matching in a bipartite graph
# for min weight matching, simply negate the weights.
import copy
import unittest
import random
import numpy as np
import pickle
""" Global variables:
       n = number of vertices on each side
       U,V vertex sets
       lu,lv are the labels of U and V resp.
       the matching is encoded as 
       - a mapping Mu from U to V, 
       - and Mv from V to U.
    
    The algorithm repeatedly builds an alternating tree, rooted in a
    free vertex u0. S is the set of vertices in U covered by the tree.
    For every vertex v, T[v] is the parent in the tree and Mv[v] the
    child.

    The algorithm maintains minSlack, s.t. for every vertex v not in
    T, minSlack[v]=(val,u1), where val is the minimum slack
    lu[u]+lv[v]-w[u][v] over u in S, and u1 is the vertex that
    realizes this minimum.

    Complexity is O(n^3), because there are n iterations in
    maxWeightMatching, and each call to augment costs O(n^2). This is
    because augment() makes at most n iterations itself, and each
    updating of minSlack costs O(n).
"""

class Dynamic_Hungarian_Algorithm:
    def __init__(self,weights) :
        self.H=Hungarian(weights)
        
    def modify_row(self,i,new_row):
        for j in range(self.H.n):
            self.H.w[i][j]=new_row[j]
        try:
            v=self.H.Mu[i]
        except KeyError:
            self.H.lu[i]=max([self.H.w[i][v]-self.H.lv[v] for v in self.H.V])
            return
        del self.H.Mu[i]
        del self.H.Mv[v]
        self.H.lu[i]=max([self.H.w[i][v]-self.H.lv[v] for v in self.H.V])
    
    def modify_col(self,i,new_col):
        for j in range( self.H.n):
            self.H.w[j][i]=new_col[j]
        try:
            u= self.H.Mv[i]
        except KeyError:
            self.H.lv[i]=max([ self.H.w[u][i]- self.H.lu[u] for u in  self.H.U])
            return
        del  self.H.Mu[u]
        del  self.H.Mv[i]
        self.H.lv[i]=max([ self.H.w[u][i] - self.H.lu[u] for u in  self.H.U])
        
    def modify_val(self,i,j):
        pass
    
    def max_dynamic_cal(self):
        while len(self.H.Mu)<self.H.n:
            free = [u for u in self.H.U if u not in self.H.Mu]      # choose free vertex u0
            u0 = free[0]
            self.H.S = {u0: True}                            # grow tree from u0 on
            self.H.T = {}
            self.H.minSlack = [[self.H.slack(u0,v), u0] for v in self.H.V]
            self.H.augment()
        #                                    val. of matching is total edge weight
        val = sum(self.H.lu)+sum(self.H.lv)
    
        #return (Mu, Mv,lu,lv, val)
        return val
        
    def min_dynamic_cal(self):
        self.H.negate()
        while len(self.H.Mu)<self.H.n:
            #print("#")
            free = [u for u in self.H.U if u not in self.H.Mu]      # choose free vertex u0
            u0 = free[0]
            self.H.S = {u0: True}                            # grow tree from u0 on
            self.H.T = {}
            self.H.minSlack = [[self.H.slack(u0,v), u0] for v in self.H.V]
            self.H.augment()
        #                                    val. of matching is total edge weight
        val = sum(self.H.lu)+sum(self.H.lv)
        self.H.negate()
        #return (Mu, Mv,lu,lv, val)
        return -val
        
    def copy_Hungarian(self):
        new_H=Dynamic_Hungarian_Algorithm(self.H.w)
        new_H.H=self.H.copy_Hungarian()
        return new_H
        

class Hungarian:
    def __init__(self,weights) :
        #self.init_max(weights)
        self.w  = weights
        self.n  = len(self.w)
        self.U  = self.V = range(self.n)
    
    
    def negate(self):
        for i in range(self.n):
            for j in range(self.n):
                self.w[i][j]=-self.w[i][j]
                
    def improveLabels(self,val):
        """ change the labels, and maintain minSlack.
        """
        for u in self.S:
            self.lu[u] -= val
        for v in self.V:
            if v in self.T:
                self.lv[v] += val
            else:
                self.minSlack[v][0] -= val
            
    def improveMatching(self,v):
        """ apply the alternating path from v to the root in the tree.
        """
        u = self.T[v]
        if u in self.Mu:
            self.improveMatching(self.Mu[u])
        self.Mu[u] = v
        self.Mv[v] = u
        
    def slack(self,u,v):
        return self.lu[u]+self.lv[v]-self.w[u][v]
    
    
    
    def augment(self):
        """ augment the matching, possibly improving the lablels on the way.
        """
            #global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
        while True:
            # select edge (u,v) with u in S, v not in T and min slack
        
            ((val, u), v) = min([(self.minSlack[v], v) for v in self.V if v not in self.T])
        
            assert u in self.S
            if val>0:
                self.improveLabels(val)
        # now we are sure that (u,v) is saturated
            #print(self.lu[u],self.lv[v],self.w[u][v],self.slack(u,v))
            assert self.slack(u,v)==0
            self.T[v] = u                            # add (u,v) to the tree
            if v in self.Mv:
                u1 = self.Mv[v]                      # matched edge,
                assert not u1 in self.S
                self.S[u1] = True                    # ... add endpoint to tree
                for v in self.V:                     # maintain minSlack
                    if not v in self.T and self.minSlack[v][0] > self.slack(u1,v):
                        self.minSlack[v] = [self.slack(u1,v), u1]
            else:
                self.improveMatching(v)              # v is a free vertex
                return

    def WeightMatching(self):
        #global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
        self.lu = [ max([self.w[u][v] for v in self.V]) for u in self.U]  # start with trivial labels
        self.lv = [ 0                         for v in self.V]
        self.Mu = {}                                       # start with empty matching
        self.Mv = {}
        while len(self.Mu)<self.n:
            
            free = [u for u in self.U if u not in self.Mu]      # choose free vertex u0
            u0 = free[0]
            self.S = {u0: True}                            # grow tree from u0 on
            self.T = {}
            self.minSlack = [[self.slack(u0,v), u0] for v in self.V]
            self.augment()
        # val. of matching is total edge weight
        val = sum(self.lu)+sum(self.lv)
        #return (Mu, Mv,lu,lv, val)
        return val
        
    def minWeightMatching(self):
        self.negate()
        ans=-self.WeightMatching()
        self.negate()
        return ans
        
    def maxWeightMatching(self):
        return self.WeightMatching()
#copy.deepcopy(self.weights)
    def copy_Hungarian(self):
        #new_hungarian=Hungarian(copy.deepcopy(self.w))
        new_hungarian=Hungarian(pickle.loads(pickle.dumps(self.w, -1)))
        new_hungarian.lu=pickle.loads(pickle.dumps(self.lu, -1))
        #new_hungarian.lu=copy.deepcopy(self.lu)
        new_hungarian.lv=pickle.loads(pickle.dumps(self.lv, -1))
        #new_hungarian.lv=copy.deepcopy(self.lv)
        new_hungarian.Mu=pickle.loads(pickle.dumps(self.Mu, -1))
        #new_hungarian.Mu=copy.deepcopy(self.Mu)
        new_hungarian.Mv=pickle.loads(pickle.dumps(self.Mv, -1))
        #new_hungarian.Mv=copy.deepcopy(self.Mv)
        new_hungarian.S=pickle.loads(pickle.dumps(self.S, -1))
        #new_hungarian.S=copy.deepcopy(self.S)
        new_hungarian.T=pickle.loads(pickle.dumps(self.T, -1))
        #new_hungarian.T=copy.deepcopy(self.T)
        new_hungarian.minSlack=pickle.loads(pickle.dumps(self.minSlack, -1))
        #new_hungarian.minSlack=copy.deepcopy(self.minSlack)
        return new_hungarian
        


def improveLabels(val):
    """ change the labels, and maintain minSlack. 
    """
    for u in S:
        lu[u] -= val
    for v in V:
        if v in T:
            lv[v] += val
        else:
            minSlack[v][0] -= val

def improveMatching(v):
    """ apply the alternating path from v to the root in the tree. 
    """
    u = T[v]
    if u in Mu:
        improveMatching(Mu[u])
    Mu[u] = v
    Mv[v] = u

def slack(u,v): return lu[u]+lv[v]-w[u][v]

def augment():
    """ augment the matching, possibly improving the lablels on the way.
    """
    while True:
        # select edge (u,v) with u in S, v not in T and min slack
        
        ((val, u), v) = min([(minSlack[v], v) for v in V if v not in T])
        
        assert u in S
        if val>0:
            improveLabels(val)
        # now we are sure that (u,v) is saturated
        assert slack(u,v)==0
        T[v] = u                            # add (u,v) to the tree
        if v in Mv:
            u1 = Mv[v]                      # matched edge, 
            assert not u1 in S
            S[u1] = True                    # ... add endpoint to tree 
            for v in V:                     # maintain minSlack
                if not v in T and minSlack[v][0] > slack(u1,v):
                    minSlack[v] = [slack(u1,v), u1]
        else:
            
            improveMatching(v)              # v is a free vertex
           
            return

def maxWeightMatching(weights):
    """ given w, the weight matrix of a complete bipartite graph,
        returns the mappings Mu : U->V ,Mv : V->U encoding the matching
        as well as the value of it.
    """
    global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
    w  = weights
    n  = len(w)
    U  = V = range(n)
    lu = [ max([w[u][v] for v in V]) for u in U]  # start with trivial labels
    lv = [ 0                         for v in V]
    Mu = {}                                       # start with empty matching
    Mv = {}
    while len(Mu)<n:
        #print("1")
        free = [u for u in U if u not in Mu]      # choose free vertex u0
        u0 = free[0]
        S = {u0: True}                            # grow tree from u0 on
        T = {}
        minSlack = [[slack(u0,v), u0] for v in V]
     
        augment()
    #                                    val. of matching is total edge weight
    val = sum(lu)+sum(lv)
    
    #return (Mu, Mv,lu,lv, val)
    return val
    
    
def Dynamic_row_change(i,new_row):
    global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
    n=len(w)
    for j in range(n):
        w[i][j]=new_row[j]
    v=Mu[i]
    del Mu[i]
    del Mv[v]
    
    #lu[i]=max([w[i][v] for v in V])
   
    lu[i]=max([w[i][v]-lv[v] for v in V])
   
    
    while len(Mu)<n:
        #print("1")
        free = [u for u in U if u not in Mu]      # choose free vertex u0
        u0 = free[0]
      
        S = {u0: True}                            # grow tree from u0 on
       
        T = {}
     

        minSlack = [[slack(u0,v), u0] for v in V]
        
        augment()
    #                                    val. of matching is total edge weight
    val = sum(lu)+sum(lv)
    
    #return (Mu, Mv,lu,lv, val)
    return val
    
def Dynamic_col_change(i,new_col):
    global U,V,S,T,Mu,Mv,lu,lv, minSlack, w
    n=len(w)
    for j in range(n):
        w[j][i]=new_col[j]
    u=Mv[i]
    del Mu[u]
    del Mv[i]
    lv[i]=max([w[u][i]-lu[u] for u in U])
    while len(Mu)<n:
        #print("1")
        free = [u for u in U if u not in Mu]      # choose free vertex u0
        u0 = free[0]
        S = {u0: True}                            # grow tree from u0 on
        T = {}
        minSlack = [[slack(u0,v), u0] for v in V]
        augment()
    #                                    val. of matching is total edge weight
    val = sum(lu)+sum(lv)
    
    #return (Mu, Mv,lu,lv, val)
    return val




if __name__ == '__main__':
    #  a small example
    M=[[3,5,5,4,10],[2,2,1,2,2],[2,4,4,1,0],[0,1,1,0,0],[1,2,1,3,3]]
    M=np.asarray(M)
    new_row=[3,12,5,12,1,3]
    print("M val =",maxWeightMatching(M))
    print(Dynamic_row_change(2,new_row))
    print(Dynamic_col_change(3,new_row))

    M=[[3,5,5,4,10],[3,3,3,3,2],[2,4,4,1,0],[0,1,1,0,0],[1,2,1,3,3]]
    DH=Dynamic_Hungarian_Algorithm(M)
    print(DH.H.minWeightMatching())
    print(DH.H.Mu)
    new_row=[0,100000,100000,100000,100000,100000]

    new_col=[100000,100000,100000,100000,0,100000]
    DH.modify_row(0,new_row)
    print(DH.H.Mu)
    DH.modify_col(4,new_col)
    print("DHM=",DH.min_dynamic_cal())
