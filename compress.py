#!/usr/bin/python3
import math
import os
import sys
import time
from queue import PriorityQueue
import platform
plt = platform.system()
Result = []
ans=[]
graph = []
vis =[]
wieghts={}
BITSINNUM=12
BITSIID = 12
BITPERBASE = 3
PER_CHARACTER_BIT=0#
SKIPLIMIT = 1
PER_CHARACTER_BIT2=0#
final_final_string=''#
final_final_string2=''


def check_mates(dict,r,t,file,reference_index,non_reference_index):
    final_final_string=''
    final_final_string=str(reference_index)+","+str(non_reference_index)#cuto
    for key in dict:#key holo kon value
        final_final_string=final_final_string+key+'#'#cuto
     
        prev=0
        for j in dict[key].split(','):
            if(prev==0):
                prev=int(j)
                num=prev
            else:
                num=int(j)-prev
                prev=int(j)
            final_final_string=final_final_string+str(num)+','#cuto
        final_final_string=final_final_string[:len(final_final_string)-1]
        final_final_string+='#'
    final_final_string+='\n'
    file.write(final_final_string)

def loadData(fileName):
    f = open(fileName, 'r')
    Data = []
    Name = []
    a=0#
    ind=0#
    filetype='p'
    data=f.read().split("\n")
    if(">" in data[0]):
    	filetype='f'
    if(filetype=='f'): 
        f = open(fileName, 'r')
        data=f.read().split(">")  
        for i in data[1:]:
            j=i.split("\n")
            s=[]
            Name.append(j[0])
            for k in j[1:]:
                s.extend(k)
            Data.append(''.join(s).strip())
            ind+=1
   
       
        len_check = [ len(i) for i in Data]
        if min(len_check)!= max(len_check):
            print("FATAL ERROR > LENGTHS NOT EQUAL ")
            raise Exception("FATAL ERROR > LENGTHS NOT EQUAL ")
        f.close()
        return len(Data), Data,Name#
    elif(filetype=='p'):

        lenn=len(data)
        first=len(data[0].split())
        i=len(data[1].split())
        if(i!=2 or first!=2):
            print("Not a fasta or phylip file ")
            raise Exception("FATAL ERROR > LENGTHS NOT EQUAL ")
        l1=(data[0].split())[0]
        l2=(data[0].split())[1]
        if(l1.isdigit()==False or l2.isdigit()==False):
            print(" Not a fasta or phylip file ")
            raise Exception("FATAL ERROR > NOT A PHYLIP FILE")
        for i in data[1:]:
            if(i=='\n' or i==''):
                break
            species_name=i.split()[0]
            species_data=i.split()[1]
            Data.append(species_data)
            Name.append(species_name)
            #print("data load",number)
            #number+=1
     
        len_check = [ len(i) for i in Data]
        print(len(len_check))
        if min(len_check)!= max(len_check):
            print("FATAL ERROR > LENGTHS NOT EQUAL ")
            raise Exception("FATAL ERROR > LENGTHS NOT EQUAL ")
        f.close()
        return len(Data), Data,Name#
    raise Exception('please print appropriate file type')

def getProbability(Data, interval):
    length = len(Data[0])
    probabilities = [dict() for i in range(length)]
    num = len(Data)
    actual_lenth = len(Data[0])
    for i in range(actual_lenth):
        for j in range(len(Data)):
            if i+interval<=len(Data[0]):
                c = Data[j][i:i+interval]
            else :
                c = Data[j][i:] 
            if c in probabilities[i]:
                probabilities[i][c]+=1
            else:
                probabilities[i][c]=1
        for c in probabilities[i]:
            probabilities[i][c]=probabilities[i][c]/num
    return probabilities

def getExpectations(Data,probabilities,interval):
    expectation = []
    for i in Data:
        expec = 0
        #print(i)
        for j in range(len(i)) :
            #print(i[j])
            if j+interval<=len(i):
                c=i[j:j+interval]
            else:
                c=i[j:]
            expec += math.log(probabilities[j][c])
        
        expectation.append(expec)
    return expectation


def getExpectations_2(Data,probabilities,interval):
    l=[]
    distance_2=[]
    
    
    for i in probabilities:
        l.append( max(i,key=i.get))
    
    k=0
    for i in Data:
        #dis=0
        dis2=0
        for j in range(len(i)):
            if( l[j]!=i[j]):
                #dis+=1
                dis2+=1-probabilities[j][i[j]]
        #distance.append(dis)
        distance_2.append(dis2)
        k+=1
    orderedpairs_3 = [ (i,j) for  j,i in enumerate(distance_2) ]
    orderedpairs_3.sort()
    order_3 = [ i for j,i in orderedpairs_3]
    return order_3


def getDiff(r, t):
    Map = {}
    dict={}#
    j=-1
    for i in range(len(r)):
        if i<=j:
            continue
        if r[i] != t[i]:
            j = i
            while i < len(r) and (r[i] != t[i] or ((i+1)<len(r) and r[i:i+SKIPLIMIT]!=t[i:i+SKIPLIMIT])) :
                if( r[i]==t[i]):
                    if(t[j:i] in Map.keys()): break
                i += 1
            subs = t[j:i]
            index=j#
            j=i
            if subs in Map.keys():
                Map[subs] += 1
                dict[subs]+=","+str(index)#
            else:
                Map[subs] = 1
                dict[subs]=str(index)#
    cost = 0
    for j in Map.keys():
       # print(j)
        cost += len(j) * BITPERBASE + (Map[j]) * (BITSINNUM)
    return cost,dict#

def getDiff_2(r, t):
    Map = {}
    dict={}#
    j=-1 
    for i in range(len(r)):
        if i<=j:
            continue
        if r[i] != t[i]:
            j = i
            while i < len(r) and r[i] != t[i]:
                i += 1
            subs = t[j:i]
            index=j#
            j=i
            if subs in Map.keys():
                Map[subs] += 1
                dict[subs]+=","+str(index)#
            else:
                Map[subs] = 1
                dict[subs]=str(index)#
    cost = 0
    for j in Map.keys():
       # print(j)
        cost += len(j) * BITPERBASE + (Map[j]) * (BITSINNUM)
    return cost,dict#

def getReduceVal(s):
    global BITSINNUM
    cost=0
    count=0
    for c in s:
        if c!='-':
            if count>=5:
                count=0
                cost+=BITSINNUM
            cost+=3
        else:
            count+=1
    if count>=5:
        count=0
        cost+=BITSINNUM+6
    else:
        cost+=count*3
    return cost/2


def getmatrix(Data, num, length, overlap,order1):
    weights = {}
    dict = {}
    #for i in range(num):
    	#for j in range(num):
    		#weights[i,j]=-1

    for i in range(0,num-length,length-overlap):
        if(i+length>num):
            m=order1[i:]
            #m2=order2[i:]
        else:
            m=order1[i:i+length]
            #m2=order1[i:i+length]
        for i in m:
            for j in m:
                
                if(i!=j and (i,j) not in weights):
                    dict[i,j]={}#
                    weights[i,j],dict[i,j]=getDiff(Data[i], Data[j])#
    
    for i in range(num):
        if i!=num-1:
            dict[num-1,i]={}#
            weights[num-1,i],dict[num-1,i]=getDiff(Data[num-1],Data[i])#
    return weights,dict#



def compressDataExpectation(fileName, window, overlap):
    global BITSINNUM
    global BITSIID
    global final_final_string
    global final_final_string2
    start=time.time()
    Name=[]
    num, Data,Name = loadData(fileName)#
    BITSIID = int(math.log10(num))*4
    BITSINNUM = int(math.log10(len(Data[0])))*4
  
    prob = getProbability(Data,1)
   

    order3 = getExpectations_2(Data,prob,1)
    
    own = []
    for d in Data:
        own.append(getReduceVal(d))

    weights1,dict = getmatrix(Data,num,window,overlap,order3)
    del order3[:]
   
    Edges = []
    for i in range(num):
        for j in range(num):
            
            if(i==j):
                Edges.append((0,i+1,own[i],0,i+1))
            elif (i,j) not in weights1:
                continue
            else:
                Edges.append((i+1,j+1,weights1[i,j],i+1,j+1))
    solution, refmap2 =dmst(num+1,Edges,0)
    del own[:]
    
    del weights1
   
    import os
    if plt == "Windows":    
        FolderName=fileName+'.chapao\\'
    elif plt == "Linux":
        FolderName=fileName+'.chapao/'
    
    os.makedirs(os.path.dirname(FolderName), exist_ok=True)
    fileName=FolderName+"ref.txt"
    fileName2= FolderName+"metadata.txt"

    file=open(fileName,'w')
    file2=open(fileName2,'w')

    
    total_ref=0
    total_non_ref=0
    for i in refmap2:
        #print( i[0],i[1])
        j=i[0]
        k=i[1]
     
        
        if(i[0]==0):
            total_ref+=1
           
            file.write(str(i[1]-1)+','+Data[i[1]-1]+"\n")
        else :
            total_non_ref+=1
            check_mates(dict[i[0]-1,i[1]-1],Data[i[0]-1],Data[i[1]-1],file2,i[0]-1,i[1]-1)#,Data[i[3
            
    fff=str(num)
    for i in Name:
        fff=fff+"|"+i
    fff+="\n"
    file.write(fff)
  
    file2.write(final_final_string)
    file2.flush()
   

    file2.close()
    file.close()

    
    #print("ref",total_ref,"non",total_non_ref)
    import bz2
    compressionLevel=9
    tarbz2contents1 = bz2.compress(open(fileName, 'rb').read(), compressionLevel)
    tarbz2contents2 = bz2.compress(open(fileName2, 'rb').read(), compressionLevel)

    file=open(fileName,'wb')#
    file2=open(fileName2,'wb')#
    file.write(tarbz2contents1)
    file2.write(tarbz2contents2)


    file.flush()
    file2.flush()

    file.close()
    file2.close()


    size1_7 = sys.getsizeof(tarbz2contents1)
    size2_2 = sys.getsizeof(tarbz2contents2)
    size3_72 = size1_7 + size2_2
    print('compressed size '+str(size3_72))
    #print("ref",total_ref,"non",total_non_ref)


def dfs(s,s1):
    q=PriorityQueue()
    q.put((s,0,0))
    while not q.empty():
        w,v,u=q.get()
        if vis[v]==1:
            continue
        vis[v]=1
        if(v!=u):
            ans.append((u,v))
        for i in graph[v]:
                q.put((w+wieghts[v,i],i,v,))
                
    






def addEdge(e):
    #print('addedge ', e)
    '''for i in Result:
        if(i[4] == e[4]):
            #print("remove",i[0], i[1], i[3], i[4])
            Result.remove(i)'''

    Result.append(e)


def dmst(vertices, test, root=0):
    oo = int(1e9)
    #print(test)
    cost = {}
    global graph
    global vis
    global wieghts
    #bigmap={}
    back = {}
    label = []
    bio = []
    edge = {}
    ret = 0
    nodes = vertices
    #for e in test:
        #bigmap[(e[3],e[4])]=[]
    #wieghts=[[100000000000 for i in range(vertices)] for j in range(vertices)]

    for j in test:
        #wieghts[j[0]][j[1]]=j[2]
        wieghts[j[0],j[1]]=j[2]
    while(True):
        for i in range(vertices):
            cost[i] = oo

        for e in test:
            if(e[0] == e[1]):
                continue
            if(e[2] < cost[e[1]]):
                #print("found", e[1])
                cost[e[1]] = e[2]
                back[e[1]] = e[0]
                edge[e[1]] = e

        cost[root] = 0

        for i in range(vertices):
            if(cost[i] == oo):
                print('problem at ', i)
                return -1


        for i in range(vertices):
            ret += cost[i]
            if(i!=root):
                #print("add   ",edge[i][0], edge[i][1],edge[i][3], edge[i][4])
                addEdge(edge[i])

        K = 0
        label=[-1 for i in range(vertices)]
        bio = [-1 for i in range(vertices)]
        '''for i in range(vertices):
            label[i] = -1
            bio[i] = -1'''

        for i in range(vertices):
            x = i
            while(True):
                if(x != root and bio[x] == -1):
                    bio[x] = i
                    x = back[x]
                else:
                    break

            if(x != root and bio[x] == i):
                while(True):
                    if(label[x] == -1):
                        label[x] = K
                        x = back[x]
                    else:
                        break
                K = K + 1

        if(K==0):
            break

        for i in range(vertices):
            if(label[i] == -1):
                label[i] = K
                K = K + 1


        for e in range(len(test)):
            xx = label[test[e][0]]
            yy = label[test[e][1]]
            if(xx != yy):
                zz = test[e][2] - cost[test[e][1]]
                #m=edge[test[e][1]]
                #bigmap[(test[e][3],test[e][4])].append((m[3],m[4]))
            else:
                zz = test[e][2]
                #test[e] = (test[e][0],test[e][1], zz, test[e][3],test[e][4])

            test[e] = (xx, yy,zz,test[e][3],test[e][4])
            #test[e] = (test[e][0], yy, test[e][2],test[e][3],test[e][4])

        root = label[root]
        vertices = K
    
    
    result= [(i[3], i[4]) for i in Result]
    Result.clear()
    answer=[]
    
    graph = [set() for i in range(nodes)]
    vis = [ 0 for i in range(nodes)]
    #print( len(result))
    for i in result:
        #print(i[0],i[1])
        graph[i[0]].add(i[1])
    dfs(0,0)
    answer = ans[:]
    ans.clear()
    #print( "answe length ", len(answer))

    return ret, answer
if __name__ == "__main__":
    filename=sys.argv[1]
    window = int(sys.argv[2])
    overlap = int(sys.argv[3]) 
    start_time=time.time()	
    compressDataExpectation(filename,window,overlap)
    end_time=time.time()
    print('time '+str(end_time-start_time))
  
