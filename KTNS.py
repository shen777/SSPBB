#Keep tool needest soonest algoritms
import instance
import copy
from operator import itemgetter

def ktns(solution,matrix,capacity):
    #O(MN)
    assert len(solution)>0
    N=len(solution)
    M=len(matrix[0])
    L=pre(solution,matrix)
    #step 1
    J=[0]*M
    Li0=[]
    #print(L)
    for i in range(M):
        #print(i)
        temp=[L[i][0],i]
        Li0.append(temp)
    Li0=sorted(Li0,key=itemgetter(0))
    for i in range(capacity):
        J[Li0[i][1]]=1
    n=1
    W=[]
    #step 2
    while(n<N):
        W.append(copy.deepcopy(J))
        flag=1
        #step 3
        for i in range(M):
            if L[i][n]==n and J[i]==0:
                flag=0
        if flag:
            n=n+1
        else:
            #step 4
            for i in range(M):
                if L[i][n]==n and J[i]==0:
                    J[i]=1
            Lik=[]
            #step 5
            count=0
            for i in range(M):
                if J[i]==1:
                    count+=1
                    temp=[L[i][n],i]
                    Lik.append(temp)
            Lik=sorted(Lik,key=itemgetter(0), reverse=True)
            for i in range(count-capacity):
                J[Lik[i][1]]=0
            n=n+1
    W.append(copy.deepcopy(J))
    #print("######",W)
    return W


def pre(solution,matrix):
    m=len(matrix[0])
    n=len(solution)
    T = []
    L=[]
    for i in range(m):
        set={n}
        T.append(copy.deepcopy(set))
        temp=[0]*n
        L.append(temp)
        
        
    for i in reversed(range(n)):
        for j in range(m):
            if matrix[solution[i]][j]==1:
                T[j].add(i)
            L[j][i]=min(T[j])
    return L

def KTNS_Cost(solution,k):
    if len(solution)<=0:
        return 0
    
    W=ktns(solution,k.matrix,k.capacity)
   
    cost=0
    for i in range(len(W[0])):
        #print(W[0][i])
        cost+=W[0][i]
    for i in range(1,len(solution)):
        for j in range(len(W[0])):
            if W[i-1][j]==0 and W[i][j]==1:
                cost+=1
    return cost

    
if __name__ == '__main__':
    solution=[0,1,2,3,4,5,6,7,8,9]
    k=instance.load_data("datA1")
    #print(k.matrix)
   
    #print(KTNS_Cost(solution,k))
    """
    for i in range(1,11):
        print("solution=",range(i))
        W=ktns(range(i),k.matrix,k.capacity)
        #print("ans")
        print(W)
        print(KTNS_Cost(range(i),k))
    """
