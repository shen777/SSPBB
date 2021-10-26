# jobs and tools
import numpy as np
class instance:
    def __init__(self,name,job_num,tool_num,capacity,matrix):
        self.name=name
        self.capacity=capacity
        self.job_num=job_num
        self.tool_num=tool_num
        self.matrix=np.asarray(matrix)
        
def load_data(name):
    f = open(name, "r")
    k=f.readlines()
    job_num=int(k[0])
    tool_num=int(k[1])
    capacity=int(k[2])
    #print(job_num,tool_num,capacity)
    matrix=[]
    for i in range(job_num):
        col=[0]*tool_num
        for j in range(tool_num):
            col[j]=int(k[j+3][i*2])
        matrix.append(col)
    loadong_ins=instance(name,job_num,tool_num,capacity,matrix)
    return loadong_ins
        
        
if __name__ == '__main__':
    k=load_data("datA1")
    
