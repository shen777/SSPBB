import instance
import KTNS
import DynamicHungarianAlgorithm
import copy
import pickle
inf = 1000
class Edge :
   def __init__(self, arg_src, arg_dst, arg_weight) :
       self.src = arg_src
       self.dst = arg_dst
       self.weight = arg_weight

class Graph :
    def __init__(self, arg_num_nodes, arg_edgelist,i) :
        self.num_nodes = arg_num_nodes
        self.edgelist  = arg_edgelist
        self.parent    = []
        self.rank      = []
        # mst stores edges of the minimum spanning tree
        self.mst       = []
        self.i=i
    def FindParent (self, node) :
        # With path-compression.
        if node != self.parent[node] :
            self.parent[node] = self.FindParent(self.parent[node])
        return self.parent[node]

        # Without path compression
        # if node == self.parent[node] :
        #    return node
        # return self.FindParent(self.parent[node])

    def KruskalMST (self) :
        
        # Sort objects of an Edge class based on attribute (weight)
        self.edgelist.sort(key = lambda Edge : Edge.weight)

        self.parent = [None] * self.num_nodes
        self.rank   = [None] * self.num_nodes

        for n in range(self.num_nodes) :
            self.parent[n] = n # Every node is the parent of itself at the beginning
            self.rank[n] = 0   # Rank of every node is 0 at the beginning

        for edge in self.edgelist :
            root1 = self.FindParent(edge.src)
            root2 = self.FindParent(edge.dst)

            # Parents of the source and destination nodes are not in the same subset
            # Add the edge to the spanning tree
            if root1 != root2 :
               self.mst.append(edge)
               if self.rank[root1] < self.rank[root2] :
                  self.parent[root1] = root2
                  self.rank[root2] += 1
               else :
                  self.parent[root2] = root1
                  self.rank[root1] += 1
        cost = 0
        for edge in self.mst :
            cost += edge.weight
        
        return cost
        
def lb2(L,i,solution) :
    # start from solution i+1's spanning tree
    #O(n^2 log(n))
    # Edge(source, destination, weight)
    #assert i<len(solution)-1
    if i>=len(solution)-1:
        return 0
    
    num_nodes = len(solution)-i-1
    l=[]
    for j in range(i+1,len(solution)):
        for k in range(i+1,len(solution)):
            if j!=k:
                e=Edge(j-i-1, k-i-1, L[solution[j]][solution[k]])
                l.append(e)
    
        
    g1 = Graph(num_nodes,l,i)
    cost=g1.KruskalMST()
    if i<0:
        return cost
    closest=L[solution[i]][solution[i+1]]
    for j in range(i+1,len(solution)):
        closest=min(closest,L[solution[i]][j])
    #print("closest",closest)
    cost+=closest
    return cost

def lb1(i,solution,matrix,capacity):
    # start from solution i=jp
    
    T=set()
    for j in range(i,len(solution)):
        if j<0:
            continue
        for k in range(len(matrix[0])):
            if matrix[solution[j]][k]==1:
                T.add(k)
    return len(T)-capacity

def lb3(i,solution,DH,L):
    if i==-1:
        La=copy.deepcopy(L)
        La=augment_matrix(La,solution)
        DH=DynamicHungarianAlgorithm.Dynamic_Hungarian_Algorithm(La)
        val=DH.H.minWeightMatching()
        return val,DH
    if i==0:
        return 0,DH
    if i>=len(solution)-1:
        return 0,None

    new_DH=copy.deepcopy(DH)
    #new_DH=DH.copy_Hungarian()
    new_row=[inf]*len(solution)
    new_row[solution[i]]=0
    new_col=[inf]*len(solution)
    new_col[solution[i-1]]=0
    
    new_DH.modify_col(solution[i],new_col)
    new_DH.modify_row(solution[i-1],new_row)
    val=new_DH.min_dynamic_cal()
    return val,new_DH
        
def lb3_m(i,solution,DH,L):
    if i==-1:
        #La=copy.deepcopy(L)
        La=pickle.loads(pickle.dumps(L, -1))
        La=augment_matrix(La,solution)
        DH=DynamicHungarianAlgorithm.Dynamic_Hungarian_Algorithm(La)
        val=DH.H.minWeightMatching()
        return val,DH
    if i==0:
        return 0,DH
    if i>=len(solution)-1:
        return 0,None

    #new_DH=copy.deepcopy(DH)
    new_DH=DH.copy_Hungarian()
    new_row=[inf]*len(solution)
    new_row[solution[i]]=0
    new_col=[inf]*len(solution)
    new_col[solution[i-1]]=0
    
    new_DH.modify_col(solution[i],new_col)
    new_DH.modify_row(solution[i-1],new_row)
    val=new_DH.min_dynamic_cal()
    return val,new_DH
    
def Greedy_lb3(i,solution,DH,L):
    if i==0:
        #La=copy.deepcopy(L)
        La=pickle.loads(pickle.dumps(L, -1))
        La=augment_matrix(La,solution)
        DH=DynamicHungarianAlgorithm.Dynamic_Hungarian_Algorithm(La)
        val=DH.H.minWeightMatching()
        return val,DH
    if i>=len(solution)-1:
        return 0,None

    #new_DH=copy.deepcopy(DH)
    new_DH=DH.copy_Hungarian()
    new_row=[inf]*len(solution)
    new_row[solution[i]]=0
    new_col=[inf]*len(solution)
    new_col[solution[i-1]]=0
    
    new_DH.modify_col(solution[i],new_col)
    new_DH.modify_row(solution[i-1],new_row)
    val=new_DH.min_dynamic_cal()
    return val,new_DH

def lb4(i,solution):
    if 1 in solution[:i] and 0 in solution[i:]:
        return True
    return False
 
def lij(matrix,capacity):
    L=[]
   
    for i in range(len(matrix)):
        l=[0]*len(matrix)
        for j in range(len(matrix)):
            count=0
            for k in range(len(matrix[0])):
                if matrix[i][k]==1 or matrix[j][k]:
                    count+=1
            l[j]=max(0,count-capacity)
        l[i]=inf
        L.append(l)
    #print(L)
    return L
    
def augment_matrix(matrix,solution):
    for i in range(len(matrix)):
        matrix[i].append(0)
    matrix.append([0]*(len(matrix)+1))
    matrix[len(matrix)-1][len(matrix)-1]=inf
    for i in range(len(matrix)-1):
        matrix[i][solution[0]]=inf
    return matrix
    

if __name__ == '__main__':
    inf = 1000
    k=instance.load_data("dat1")
    solution=list(range(len(k.matrix)))
    W=KTNS.ktns(solution,k.matrix,k.capacity)
    L=lij(k.matrix,k.capacity)
    print("solution=",solution)
    L=augment_matrix(L,solution)
    solution.append(len(solution))
    DH=DynamicHungarianAlgorithm.Dynamic_Hungarian_Algorithm(L)
    print("######start",DH.H.minWeightMatching())
    print(DH.H.Mu)
    index=0
    for j in range(index,len(solution)-1):
        solution[index],solution[j]=solution[j],solution[index]
        print(solution)
        val,DH=lb3(index-1,solution,DH,L)
        print(DH.H.Mu)
        print(val)
        #self.branch(index+1,DH)
        solution[index],solution[j]=solution[j],solution[index]

    
    # dummy node need to define better
