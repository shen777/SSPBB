import random
def generate(job_num,tool_num,capacity,filename):
    L=[str(job_num)+"\n",str(tool_num)+"\n",str(capacity)+"\n"]
    M=[]
    for i in range(tool_num):
        M.append([])
    for i in range(job_num):
        choose=random.sample(list(range(tool_num)),capacity)
        for j in range(tool_num):
            if j in choose:
                M[j].append(1)
            else:
                M[j].append(0)
    #print(M)
    #filename="dat1"
    fp = open(filename, "w")
    fp.writelines(L)
    for i in range(tool_num):
        L=""
        for j in range(job_num):
            L+=str(M[i][j])+" "
        L+="\n"
        fp.writelines(L)
    fp.close()

def generate_non_full(job_num,tool_num,capacity,min_need,filename):
    L=[str(job_num)+"\n",str(tool_num)+"\n",str(capacity)+"\n"]
    M=[]
    for i in range(tool_num):
        M.append([])
    for i in range(job_num):
        need=random.randrange(min_need, capacity+1)
        choose=random.sample(list(range(tool_num)),need)
        for j in range(tool_num):
            if j in choose:
                M[j].append(1)
            else:
                M[j].append(0)
    #print(M)
    #filename="dat1"
    fp = open(filename, "w")
    fp.writelines(L)
    for i in range(tool_num):
        L=""
        for j in range(job_num):
            L+=str(M[i][j])+" "
        L+="\n"
        fp.writelines(L)
    fp.close()


if __name__ == '__main__':
    job_num=10
    tool_num=10
    capacity=5
    filename="dat1"
    generate(job_num,tool_num,capacity,filename)
    
