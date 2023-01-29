#team members
#Mohammed Mohammed Aladdin Ahmed
#Manar Ashraf Fouad
#Eslam Ayman Labib Abdelnaby
#Mennatullah Hassanin Fares
#Ibrahim Ahmed Abdel Ghaffar Mohamed
#Osama Usry Mohamed Al-Adham
#Ahmed Abdelhamed Mohamed Ali

#single read

#read file and kmers length
def readsinglefile(filename):
    with open(filename) as file:
        lines = [[line.rstrip()] for line in file]         #loop in file to read each line and strip it
        k = lines[0]                                       #specify kmers
        k = ''.join(k)                                     #convert to string
        k =int(k)                                          # convert to integer
        sequence = []
        for i in range(1, len(lines)):                     # start from 1 to skip first line as it is kmer number
            sequence.append(lines[i][0])                   #[row][col]---->[all rows except first line][col=0 ] as we have many rows but one column

        #print("single reads:", sequence)
        #print("kmer number:", k)
        return sequence, k


#create DE Bruijn Graph
def singleDBG(sequence,k):
    combination = []
    graph = []
    for seq in sequence:                                  #take k-1 from right and from left and append in combination list
        combination.append(seq[0:k-1])
        combination.append(seq[1:k])

    for i in range (0,len(combination)-1,2):               #then loop inside combination list to make graph
                                                           #that line is for ex: CTT ----> TTA then append in graph list
        graph.append([combination[i], combination[i+1]])   #to append read and its combination together as explained in the previous comment
    #print("graph:",graph)
    return graph

#help us to find start point in graph
# start point is the point that found in list 2 and not found in list 1
def single_findStartPoint(graph):
    tmp_tofind_startpoint = []
    for i in range (len(graph)):
        tmp_tofind_startpoint.append(graph[i][1])            #to find start point in second list [ as unique value]-->[all rows][list 2]
    return tmp_tofind_startpoint

#create path from tracing graph and begin with start point
def singlePath(graph,tmp_tofind_startpoint):
    path = []
    for i in range(len(graph) - 1):
        if (graph[i][0] not in tmp_tofind_startpoint):     #check that graph of list 1 not in list 2
            path.insert(0, graph[i][0])                    #that part for making the path
            path.insert(1, graph[i][1])                    # by moving from list 1 to list 2
            z = i                                          # begin with start point and saving position to move to other list
    for l in graph:                                        # and so on until we find end point or end of 2 lists
        for j in range(len(graph)):                        # for ex: GGC-->GCT-->CTT -------
            if (graph[z][1] == graph[j][0]):
                path.append(graph[j][1])
                z = j
    #print("path:",path)
    return path

#start create genome from path
def createsinglegenome(path,k):
    genome = []                         #create empty list to append genome to it
    genome.insert(0,path[0])            #take first node in path
    for i in range (1,len(path)):       #start from 1 to skip first node
        genome.append(path[i][k-2])     #take last character from each read after first node
    wholegenome = ''.join(genome)       #convert list to string
    return wholegenome



def SingleRead():       #recall all previous functions
    filename="D:\Downloads\SingleReadsInput.txt"
    sequence, k =readsinglefile(filename)
    graph=singleDBG(sequence,k)
    startpoint=single_findStartPoint(graph)
    path=singlePath(graph,startpoint)
    wholegenome=createsinglegenome(path,k)
    print("wholegenome:",wholegenome)
    return wholegenome



#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------


#pair read

#read file , kmers and gap length
def readpairfile(filename):
    with open(filename) as file:
        lines = [(line.strip()).split() for line in file]                  #loop in file to read each line, strip and split it
    k = lines[0][0]                                                        #specify kmers
    k = ''.join(k)                                                         #convert to string
    k = int(k)                                                             # convert to integer
    gap = lines[0][1]                                                      #specify gap
    gap = ''.join(gap)                                                      #convert to string
    gap = int(gap)                                                          # convert to integer
    sequence = []
    for i in range(1, len(lines)):                                           #start from 1 to skip first line
        sequence.append(lines[i][0])                                         #[row][col]---->[all rows except first line][col=0 ] as we have many rows but one column
    for i in range(0, len(lines) - 1):
        sequence[i] = "".join(ch for ch in sequence[i] if ch.isalnum())      #remove '|' from the list
    #print("pair reads:",sequence)
    #print("kmer number:",k)
    #print("gap number:",gap)
    return sequence,k,gap


#create DE Bruijn Graph
def pairDBG(sequence,k):
    combination = []

    for seq in sequence:
        combination.append(seq[0:k-1])                                 #take k-1 from right and from left and append in combination list but for both reads
        combination.append(seq[k:k+k-1])
        combination.append(seq[1:k])
        combination.append(seq[k+1:k+k])
    #print(combination)
    graph = [list(combination[i: i + 2]) for i in range(0, len(combination), 2)]   #that line is for ex: GAC ----> ACC then append in graph list
                                                                                   #                     GCG ----> CGC
    Graph = [list(graph[i: i + 2]) for i in range(0, len(graph), 2)]               #to append read and its combination together as explained in the previous comment
    #print("Graph:",Graph)
    return Graph


#help us to find start point in graph
# start point is the point that found in list 2 and not found in list 1
def pair_findStartPoint(graph):
    tmp_tofind_startpoint = []
    for i in range (len(graph)):
        tmp_tofind_startpoint.append(graph[i][1])                     #to find start point in second list [ as unique value]-->[all rows][list 2]
    return tmp_tofind_startpoint
    #print(tmp_tofind_startpoint)


#create path from tracing graph and begin with start point
def pairPath(graph,tmp_tofind_startpoint):
    path = []

    for i in range(len(graph) - 1):
        if (graph[i][0] not in tmp_tofind_startpoint):               #check that graph of list 1 not in list 2
            path.insert(0, graph[i][0])                              # that part for making the path
            path.insert(1, graph[i][1])                              # by moving from list 1 to list 2
            z = i                                                    # begin with start point and saving position to move to other list
                                                                     # and so on until we find end point or end of 2 lists
                                                                     # for ex: GAC-->ACC-->CCG -------
    for l in graph:                                                  #         GCG-->CGC-->GCC -------
        for j in range(len(graph)):
            if (graph[z][1] == graph[j][0]):
                path.append(graph[j][1])
                z = j

    #print("path:",path)
    return path


#start create genome from path
def createpairgenome(path,k,gap):
    prefix=[]
    suffix =[]

    for i in range (0,len(path)-1):
        prefix.append(path[i][0][0])               #take first character of the initial kmers
    prefix.append(path[len(path)-1][0])            #and the last one from the initial
    prefix = ''.join(prefix)                       #convert to string

    for i in range (0,len(path)-1):
        suffix.append(path[i][1][0])               #take first character of the terminal kmers
    suffix.append(path[len(path)-1][1])            #and the last one from the terminal
    suffix = ''.join(suffix)                       #convert to string

    suffixlastindex = k+gap                             #specify k+gap
    wholegenome = prefix + suffix[-suffixlastindex:]    #take whole prefix and last k+gap of suffix
    #print(prefix)
    #print(suffix)
    return wholegenome


def PairRead():                              #recall all previous functions
    filename="D:\Downloads\ReadPairsInput.txt"
    sequence,k,gap=readpairfile(filename)
    graph=pairDBG(sequence,k)
    startpoint=pair_findStartPoint(graph)
    path=pairPath(graph,startpoint)
    wholegenome=createpairgenome(path,k,gap)
    print("wholegenome:", wholegenome)
    return wholegenome


#---------------------------------------------------------------------------------------
#---------------------------------------------------------------------------------------


#user interface
Read_input=input('write the type of read (single/pair):')
if Read_input=='single':
   wholegenome=SingleRead()
elif Read_input=='pair':
   wholegenome=PairRead()








