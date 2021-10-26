import instance
import KTNS
import lb
from itertools import permutations
import time
import copy
import generate_data
import pandas as pd
import cProfile
import pickle
import numpy as np
import os
inf=1000

class Branch_and_Bound:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.node=0
    def start(self):
        ans=self.branch(0)
        return self.upper_bound
    
    def branch(self,index):
        #ip=index-1
        self.node+=1
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
            
        #bound
        
        lower_bound=KTNS.KTNS_Cost(self.solution[:max(0,index-1)],self.ins)+max(lb.lb1(index-1,self.solution,self.ins.matrix,self.ins.capacity),lb.lb2(self.L,index-1,self.solution))
        
        if lower_bound>=self.upper_bound:
            return lower_bound
            
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.branch(index+1)
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
        return lower_bound
        
class Greedy_Branch_and_Bound:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.node=0
    def start(self):
        ans=self.branch(0)
        return self.upper_bound
    
    def branch(self,index):
        #ip=index-1
        self.node+=1
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        
        s=[]
            
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            lower_bound=KTNS.KTNS_Cost(self.solution[:max(0,index)],self.ins)+max(lb.lb1(index,self.solution,self.ins.matrix,self.ins.capacity),lb.lb2(self.L,index,self.solution))
            if lower_bound<self.upper_bound:
                s.append((lower_bound,j))
            else:
                self.node+=1
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]

        s.sort()
        for j in range(len(s)):
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
            self.branch(index+1)
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
        
        return lower_bound
        
class Greedy_Symmetric_Branch_and_Bound:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.node=0
    def start(self):
        ans=self.branch(0)
        return self.upper_bound
    
    def branch(self,index):
        #ip=index-1
        self.node+=1
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        if lb.lb4(index,self.solution):
            return self.upper_bound
        s=[]
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            lower_bound=KTNS.KTNS_Cost(self.solution[:max(0,index)],self.ins)+max(lb.lb1(index,self.solution,self.ins.matrix,self.ins.capacity),lb.lb2(self.L,index,self.solution))
            if lower_bound<self.upper_bound:
                s.append((lower_bound,j))
            else:
                self.node+=1
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]

        s.sort()
        for j in range(len(s)):
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
            self.branch(index+1)
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
        
        return lower_bound

class easy_BB:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.node=0
    def start(self):
        ans=self.branch(0)
        return self.upper_bound
    
    def branch(self,index):
        self.node+=1
        #print(self.solution[:index],"     ",self.upper_bound)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
            
        #bound
        #lower_bound=KTNS.KTNS_Cost(self.solution[:index],self.ins)+max(lb.lb1(index-1,self.solution,self.ins.matrix,self.ins.capacity),lb.lb2(self.L,index-1,self.solution))
      
        lower_bound=KTNS.KTNS_Cost(self.solution[:index],self.ins)
        #print(lower_bound)
        if lower_bound>=self.upper_bound:
            return lower_bound
            
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.branch(index+1)
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
        return lower_bound

class DHungarian_BB:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.solution_a=copy.deepcopy(self.solution)
        self.solution_a.append(len(self.solution_a))
        self.node=0
        self.lb1_dominate=self.lb2_dominate=self.lb3_dominate=0
       
    def start(self):
        ans=self.branch(0,None)
        return self.upper_bound
    
    def branch(self,index,DH):
        self.node+=1
        #print(self.node)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        #bound
        
        val,new_DH=lb.lb3(index-1,self.solution_a,DH,self.L)
        lower_bound=KTNS.KTNS_Cost(self.solution[:max(index-1,0)],self.ins)+val
        
        #print(solution_a)
        #print(new_DH.H.Mu)
        #print(index,val)
       
        #lower_bound=KTNS.KTNS_Cost(self.solution[:index],self.ins)+val
        
        if lower_bound>=self.upper_bound:
            return lower_bound
        
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            self.branch(index+1,new_DH)
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
        return lower_bound

class Hungarian_BB:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.solution_a=copy.deepcopy(self.solution)
        self.solution_a.append(len(self.solution_a))
        self.node=0
        self.lb1_dominate=self.lb2_dominate=self.lb3_dominate=0
       
    def start(self):
        ans=self.branch(0,None)
        return self.upper_bound
    
    def branch(self,index,DH):
        self.node+=1
        #print(self.node)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        #bound
        
        val,new_DH=lb.lb3(index-1,self.solution_a,DH,self.L)
        lower_bound=KTNS.KTNS_Cost(self.solution[:max(index-1,0)],self.ins)+val
        
        #print(solution_a)
        #print(new_DH.H.Mu)
        #print(index,val)
       
        #lower_bound=KTNS.KTNS_Cost(self.solution[:index],self.ins)+val
        
        if lower_bound>=self.upper_bound:
            return lower_bound
        
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            self.branch(index+1,new_DH)
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
        return lower_bound

class BB_count:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        #self.solution_a=copy.deepcopy(self.solution)
        self.solution_a=pickle.loads(pickle.dumps(self.solution, -1))
        self.solution_a=np.append(self.solution_a,len(self.solution_a))
        self.node=0
        self.lb1_dominate=self.lb2_dominate=self.lb3_dominate=0
        self.lb1_lb2_tie=self.lb2_lb3_tie=self.lb1_lb3_tie=0
        self.tie=0
    def start(self):
        ans=self.branch(0,None)
        return self.upper_bound
    
    def branch(self,index,DH):
        self.node+=1
        #print(self.node)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        #bound
        kc=KTNS.KTNS_Cost(self.solution[:max(0,index-1)],self.ins)
        val,new_DH=lb.lb3_m(index-1,self.solution_a,DH,self.L)
        
        val3=kc+val
        val1=kc+lb.lb1(index-1,self.solution,self.ins.matrix,self.ins.capacity)
        val2=kc+lb.lb2(self.L,index-1,self.solution)
        
        
        if val1>val2 and val1>val3:
            self.lb1_dominate+=1
        elif val2>val1 and val2>val3:
            self.lb2_dominate+=1
        elif val3>val1 and val3>val2:
            self.lb3_dominate+=1
            
        elif val2==val1 and val2==val3:
            self.tie+=1
            
        elif val2==val1 and val2>val3:
            self.lb1_lb2_tie+=1
        elif val2==val3 and val2>val1:
            self.lb2_lb3_tie+=1
        elif val3==val1 and val1>val2:
            self.lb1_lb3_tie+=1
        else:
            print(val1,val2,val3)
        
        lower_bound=max(val1,val2,val3)
        
        #print(solution_a)
        #print(new_DH.H.Mu)
        #print(index,val)
       
        #lower_bound=KTNS.KTNS_Cost(self.solution[:index],self.ins)+val
        
        if lower_bound>=self.upper_bound:
            return lower_bound
        
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            self.branch(index+1,new_DH)
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
        return lower_bound
        
class BB_mix:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.solution_a=copy.deepcopy(self.solution)
        self.solution_a=np.append(self.solution_a,len(self.solution_a))
        self.node=0
    def start(self):
        ans=self.branch(0,None)
        return self.upper_bound
    
    def branch(self,index,DH):
        self.node+=1
        #print(self.node)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        #bound
        kc=KTNS.KTNS_Cost(self.solution[:max(0,index-1)],self.ins)
        val,new_DH=lb.lb3(index-1,self.solution_a,DH,self.L)
        
        val3=kc+val
        val1=kc+lb.lb1(index-1,self.solution,self.ins.matrix,self.ins.capacity)
        val2=kc+lb.lb2(self.L,index-1,self.solution)
        
        lower_bound=max(val1,val2,val3)
        if lower_bound>=self.upper_bound:
            return lower_bound
        
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            self.branch(index+1,new_DH)
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
        return lower_bound
        
        
class Greedy_BB_mix:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.solution_a=copy.deepcopy(self.solution)
        self.solution_a=np.append(self.solution_a,len(self.solution_a))
        self.node=0
    def start(self):
        ans=self.branch(0,None)
        return self.upper_bound
    
    def branch(self,index,DH):
        self.node+=1
        #print(self.node)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        #bound
        s=[]
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            kc=KTNS.KTNS_Cost(self.solution[:max(0,index)],self.ins)
            val,new_DH=lb.Greedy_lb3(index,self.solution_a,DH,self.L)
        
            val3=kc+val
            val1=kc+lb.lb1(index,self.solution,self.ins.matrix,self.ins.capacity)
            val2=kc+lb.lb2(self.L,index,self.solution)
        
            lower_bound=max(val1,val2,val3)
            if lower_bound<self.upper_bound:
                s.append((lower_bound,j,new_DH))
            else:
                self.node+=1
            #
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            
            
        s.sort()
        for j in range(len(s)):
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
            self.solution_a[index],self.solution_a[s[j][1]]=self.solution_a[s[j][1]],self.solution_a[index]
            self.branch(index+1,s[j][2])
            
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
            self.solution_a[index],self.solution_a[s[j][1]]=self.solution_a[s[j][1]],self.solution_a[index]
        return lower_bound
        
        
class Greedy_Symmetric_BB_mix:
    def __init__(self,solution,k):
        self.solution=solution
        self.index=0
        self.ins=k
        self.upper_bound=KTNS.KTNS_Cost(self.solution,self.ins)
        self.L=lb.lij(self.ins.matrix,self.ins.capacity)
        self.solution_a=copy.deepcopy(self.solution)
        self.solution_a=np.append(self.solution_a,len(self.solution_a))
        self.node=0
    def start(self):
        ans=self.branch(0,None)
        return self.upper_bound
    
    def branch(self,index,DH):
        self.node+=1
        #print(self.node)
        if index==len(self.solution):
            self.upper_bound=min(self.upper_bound,KTNS.KTNS_Cost(self.solution,self.ins))
            return KTNS.KTNS_Cost(self.solution,self.ins)
        #bound
        if lb.lb4(index,self.solution):
            return self.upper_bound
        s=[]
        for j in range(index,len(self.solution)):
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            kc=KTNS.KTNS_Cost(self.solution[:max(0,index)],self.ins)
            val,new_DH=lb.Greedy_lb3(index,self.solution_a,DH,self.L)
        
            val3=kc+val
            val1=kc+lb.lb1(index,self.solution,self.ins.matrix,self.ins.capacity)
            val2=kc+lb.lb2(self.L,index,self.solution)
        
            lower_bound=max(val1,val2,val3)
            if lower_bound<self.upper_bound:
                s.append((lower_bound,j,new_DH))
            else:
                self.node+=1
            #
            self.solution[index],self.solution[j]=self.solution[j],self.solution[index]
            self.solution_a[index],self.solution_a[j]=self.solution_a[j],self.solution_a[index]
            
            
        s.sort()
        for j in range(len(s)):
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
            self.solution_a[index],self.solution_a[s[j][1]]=self.solution_a[s[j][1]],self.solution_a[index]
            self.branch(index+1,s[j][2])
            
            self.solution[index],self.solution[s[j][1]]=self.solution[s[j][1]],self.solution[index]
            self.solution_a[index],self.solution_a[s[j][1]]=self.solution_a[s[j][1]],self.solution_a[index]
        return lower_bound
        

def experient(job_num,tool_num,capacity,capacity_range,filename):
    #generate_data.generate(job_num,tool_num,capacity,filename)
    generate_data.generate_non_full(job_num,tool_num,capacity,capacity-capacity_range,filename)
    #k=instance.load_data("./data/test/T1")
    k=instance.load_data(filename)
    print(k.job_num,k.tool_num,k.capacity,capacity_range)
    solution=list(range(k.job_num))
    
    print("start as",KTNS.KTNS_Cost(solution,k))
    start = time.time()
    BB=Branch_and_Bound(solution,k)
    ans1=BB.start()
    end = time.time()
    print("ans1=",ans1,"node=",BB.node)
    print("time=",end - start)
    time1=end-start
    """
    start = time.time()
    BB3=BB_mix(solution,k)
    ans3=BB3.start()
    end = time.time()
    print("ans3=",ans3,"node=",BB3.node)
    print("time=",end - start)
    time3=end-start
    assert ans1==ans3
    """
    
    start = time.time()
    BB4=BB_count(solution,k)
    ans4=BB4.start()
    

    end = time.time()
    print("ans4=",ans4,"node=",BB4.node)
    #print(BB4.lb1_dominate,BB4.lb2_dominate,BB4.lb3_dominate)
    #print(BB4.lb1_lb2_tie,BB4.lb2_lb3_tie,BB4.lb1_lb3_tie,BB4.tie)
    #L=[BB4.lb1_dominate,BB4.lb2_dominate,BB4.lb3_dominate,BB4.lb1_lb2_tie,BB4.lb2_lb3_tie,BB4.lb1_lb3_tie,BB4.tie]
    L=[BB4.lb1_dominate,BB4.lb2_dominate,BB4.lb3_dominate,BB4.lb1_lb2_tie,BB4.lb2_lb3_tie,BB4.lb1_lb3_tie,BB4.tie]
    print("time=",end - start)
    time4=end-start
    #line=[int(time1),int(time3),int(time4),BB.node,BB3.node,BB4.node]
    line=[int(time1),int(time4),BB.node,BB4.node]
    line=line+L
    return line
    
def time_experient(job_num,tool_num,capacity,capacity_range,filename):
    generate_data.generate_non_full(job_num,tool_num,capacity,capacity-capacity_range,filename)
    k=instance.load_data(filename)
    print(k.job_num,k.tool_num,k.capacity,capacity_range)
    solution=list(range(k.job_num))

    print("start as",KTNS.KTNS_Cost(solution,k))
    start = time.time()
    BB1=Greedy_Branch_and_Bound(solution,k)
    ans1=BB1.start()
    end = time.time()
    print("ans1=",ans1,"node=",BB1.node)
    print("time=",end - start)
    time1=end-start

    start = time.time()
    BB2=Greedy_BB_mix(solution,k)
    ans2=BB2.start()
    end = time.time()
    print("ans2=",ans2,"node=",BB2.node)
    print("time=",end - start)
    time2=end-start

    start = time.time()
    BB3=Greedy_Symmetric_Branch_and_Bound(solution,k)
    ans3=BB3.start()
    end = time.time()
    print("ans3=",ans3,"node=",BB3.node)
    print("time=",end - start)
    time3=end-start

    start = time.time()
    BB4=Greedy_Symmetric_BB_mix(solution,k)
    ans4=BB4.start()
    end = time.time()
    print("ans4=",ans4,"node=",BB4.node)
    print("time=",end - start)
    time4=end-start
    
    
    #line=[int(time1),int(time3),int(time4),BB.node,BB3.node,BB4.node]
    line=[int(time1),int(time2),int(time3),int(time4),BB1.node,BB2.node,BB3.node,BB4.node]
    #df = pd.DataFrame(D, columns = ['Time1','Time2','Node1','Node2'])
    #print(df)
    #df.to_csv('output10105.csv')
    #df.to_csv(str(jobnum)+"_"+str(toolnum)+"_"+str(capacity)+"_"+str(capacity_range)+"_"+'out.csv')
    return line

def test1():
    generate_data.generate(11,11,5,"Dat1")
    k=instance.load_data("Dat1")
    print(k.job_num,k.tool_num,k.capacity)
    solution=np.asarray(list(range(k.job_num)))
    start = time.time()
    BB3=Greedy_Symmetric_Branch_and_Bound(solution,k)
    ans3=BB3.start()

    end = time.time()
    print("ans=",ans3,"node=",BB3.node)
    print("time=",end - start)
    time3=end-start
    
def test2():
    #generate_data.generate(12,12,6,"Dat1")
    k=instance.load_data("Dat1")
    print(k.job_num,k.tool_num,k.capacity)
    solution=np.asarray(list(range(k.job_num)))
    start = time.time()
    BB3=Greedy_Branch_and_Bound(solution,k)
    ans3=BB3.start()

    end = time.time()
    print("ans=",ans3,"node=",BB3.node)
    print("time=",end - start)
    time3=end-start

def test3():
    #generate_data.generate(10,10,4,"Dat1")
    k=instance.load_data("Dat1")
    print(k.job_num,k.tool_num,k.capacity)
    solution=np.asarray(list(range(k.job_num)))
    start = time.time()
    BB3=Greedy_Symmetric_BB_mix(solution,k)
    ans3=BB3.start()

    end = time.time()
    print("ans=",ans3,"node=",BB3.node)
    print("time=",end - start)
    time3=end-start

def test4():
    #generate_data.generate(10,10,4,"Dat1")
    k=instance.load_data("Dat1")
    print(k.job_num,k.tool_num,k.capacity)
    solution=np.asarray(list(range(k.job_num)))
    start = time.time()
    BB3=Greedy_BB_mix(solution,k)
    ans3=BB3.start()

    end = time.time()
    print("ans=",ans3,"node=",BB3.node)
    print("time=",end - start)
    time3=end-start
    
    
if __name__ == '__main__':
    # add lb4:symmetric property and use greedy in branching
    #experient(10,10,5,"Dat1")
    #fp = open("experiment_data.txt", "w")
    #fp.write("10,10,5\n")
    #fp.write("iter Time1 Time2 Time3 Node1 Node2 Node3 lb1_dominate lb2_dominate lb3_dominate lb1_lb2_tie lb2_lb3_tie lb1_lb3_tie tie\n")
    #fp.close()
    #fp = open("experiment_data.txt", "a")
    #generate_data.generate(10,10,5,"Dat1")
    #cProfile.run('test1()')
    #cProfile.run('test2()')
    #cProfile.run('test3()')
    #cProfile.run('test4()')
    """
    jobnum=10
    toolnum=10
    capacity=5
    capacity_range=0
    D=[]
    for i in range(10):
        L=experient(jobnum,toolnum,capacity,capacity_range,"Dat1")
        D.append(L)
    #df = pd.DataFrame(D, columns = ['Time1','Time2','Time3','Node1','Node2','Node3','lb1_dominate','lb2_dominate','lb3_dominate','lb1_lb2_tie','lb2_lb3_tie','lb1_lb3_tie','tie'])
    df = pd.DataFrame(D, columns = ['Time1','Time3','Node1','Node3','lb1_dominate','lb2_dominate','lb3_dominate','lb1_lb2_tie','lb2_lb3_tie','lb1_lb3_tie','tie'])
    print(df)
    #df.to_csv('output10105.csv')
    df.to_csv(str(jobnum)+"_"+str(toolnum)+"_"+str(capacity)+"_"+str(capacity_range)+"_"+'out.csv')
    """
    
    jobnum=15
    toolnum=15
    capacity=5
    capacity_range=0
    D=[]
    for i in range(10):
        L=time_experient(jobnum,toolnum,capacity,capacity_range,"Dat1")
        D.append(L)
    #df = pd.DataFrame(D, columns = ['Time1','Time2','Time3','Node1','Node2','Node3','lb1_dominate','lb2_dominate','lb3_dominate','lb1_lb2_tie','lb2_lb3_tie','lb1_lb3_tie','tie'])
    df = pd.DataFrame(D, columns = ['Greedy+BB time','Greedy+MCPM+BB time','Greedy+Symmetric+BB time','Greedy+Symmetric+MCPM+BB','Node1','Node2','Node3','Node4'])
    print(df)
    #df.to_csv('output10105.csv')
    file_num=0
    filename="output_data/"+str(jobnum)+"_"+str(toolnum)+"_"+str(capacity)+"_"+str(capacity_range)+"_"+'greedy_symmetric_{}.csv'.format(file_num)
    while os.path.exists(filename):
        file_num+=1
        filename="output_data/"+str(jobnum)+"_"+str(toolnum)+"_"+str(capacity)+"_"+str(capacity_range)+"_"+'greedy_symmetric_{}.csv'.format(file_num)
    df.to_csv(filename)
    
