'''
This function takes in a sequence of indices and removes duplicates next to each other
Input: rSeq (list of indices)
Output: red (reduced list of indices)
'''
def rCancel(rSeq):
    red = rSeq
    j = 0 # Begins at the first element of the sequence
    while j < (len(red)-1):
        if red[j] == red[j+1]: # Checks if the current element is equal to the next element
            red1 = red[0:j]
            red2 = []
            if (j+2 <= len(red)):
                red2.extend(red[j+2:])
            red = red1 + red2 # Removes the two equal neighboring elements
            j = 0 # Begins again
        else:
            j += 1 # Advances one element
    return red

'''
This function takes in a matrix and outputs its image under matrix mutation
Indexing of the matrix begins at 1 not 0, as we are accustomed to
Input: M (matrix), w (mutation sequence len(w) >= 0)
Output: mut (mutated Matrix)
'''
def mMutation(M, w):
    length = len(w)
    r = M.nrows()
    c = M.ncols()
    mut = copy(M) # Builds a a copy of the original matrix
    if length == 0: # In case no mutation happens
        return mut
    else:
        for i in w:
            if i < 1 or i > r or i > c: # Checks to see if we are mutating at non-vertices
                print("Invalid mutation at vertex: ", i)
                return mut
    if length == 1: # Only one mutation happens
         k = w[0]
         for i in range(r):
             for j in range(c): # Standard matrix mutation here
                 if i == k-1 or j == k-1: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = -M[i,j]
                 elif M[i,k-1]*M[k-1,j] > 0: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = M[i,j] + M[i,k-1]*M[k-1,j]*(M[i,k-1]/abs(M[i,k-1])) # Has the -1's here to account for sage indexing at 0
                 else:
                     mut[i,j] = M[i,j]
    else: # Multiple mutations happen
        mut = mMutation(mMutation(M,w[0:length-1]), [w[length-1]]) # Recursively "goes up" the mutation sequence
    return mut

'''
This function takes in a matrix and a sequence representing a total order
It then returns a GIM
Indexed at 1 not 0
Input: M (skew-symmatrizable matrix), l (sequence representing linear order len(l) = M.nrows())
Output: A (GIM)
'''
def corrGIM(M, l):
    n = M.nrows()
    length = len(l)
    A = zero_matrix(n,n) # Starts A as a zero matrix
    if n != length:
        print("Invalid order")
    else:
        for i in range(n):
            for j in range(n): # Standard rules for constructing GIM
                if i == j:
                    A[i,j] = 2
                else:
                    li = l.index(i+1) # Finds the index of the first and only instance of i+1 in the order
                    lj = l.index(j+1) # Ditto ^
                    if li < lj: # Means that i+1 < j+1 in the order
                        A[i,j] = M[i,j]
                    else:
                        A[i,j] = -M[i,j]
    return A

'''
This function takes in an index i, the initial exchange matrix, the initial c-vector matrix, and a mutation sequence w
It then returns a sequence representing the indices of the r_i's
Indexing of the matrices begins at 1 not 0, as we are accustomed to
Both matrices should be square
Inputs: i (index), bM (initial exchange matrix), cM(initial c-vector matrix), w (mutation sequence len >= 0)
Output: rW (index sequence representing r_i^w)
'''
def rMutation(i, bM, cM, w):
     wPrime = []
     length = len(w)-1
     rW = []
     if i < 1 or i > bM.nrows(): # Checks to see if this is a valid index we are picking
         print ("Invalid mutation at index: ", i)
         return w
     if length == -1: # When w is the empty mutation
         rW = [i]
     elif length == 0: # This is when w is a single mutation
         k = w[0]
         if (i == k): # Identity we showed
             rW = [i]
         else:
             cK = bM[i-1,k-1] * cM[k-1, :] # Has the -1's here to account for sage indexing at 0
             if cK >= 0: # Definition of r_i^[k] from paper
                 rW = [k, i, k]
             else:
                 rW = [i]
     else:
         wPrime = w[0:length]
         k = w[length] # Hence w = wPrime[k]
         if (i == k): # Identity we showed
             rW.extend(rMutation(i, bM, cM, wPrime)) # Recursive step here
         else:
             M = mMutation(block_matrix([[bM,cM]]), wPrime) # Puts the matrix in [ B C ] form so we can perform matrix mutation on it
             # The cluster algebra mutation provided with sage does column c-vectors instead, and it was messing up my calculations for some reason
             r = bM.nrows()
             cK = M[i-1,k-1] * M[k-1, r:] # Has the -1's here to account for sage indexing at 0
             rW1 = rMutation(i,bM,cM,wPrime)
             rW2 = [] # Corresponds to the case when r_i^(w'[k]) = r_i^(w')
             if cK >= 0:
                       rW2.extend(rMutation(k,bM,cM,wPrime)) # Corresponds to the case when r_i^(w'[k]) = r_k^(w')r_i^(w')r_k^(w')
             rW = rW2 + rW1 + rW2 # Cocatenates the lists
     rW = rCancel(rW) # Automatically removes some of the r_i's when possible
     return rW

'''
This function takes in an index sequence representing r_i^w and returns a matrix representing its action under pi
Indexing begins at 1 not 0
Returns a square Matrix
A acts on elements of Gamma on the left
Inputs: iSeq (index sequence len > 0), aM (GIM from linear ordering. Don't know how to represent linear orderings yet)
Output: A (square matrix representing pi(r_i^w))
'''
def rAction(iSeq, aM):
    length = len(iSeq)
    i = iSeq[length-1] # As the r_i's act on the left, we need to evalute the last one first
    r = aM.nrows()
    c = aM.ncols()
    A = matrix(r,c,0) # Gives a fresh matrix to represent the action
    if length == 1:
        for j in range(r):
            for k in range(c): # This follows the rules set for the action of r_i given in paper
                if j == k:
                    if k == i-1: #-1 here to account for indexing
                        A[j,k] = -1
                    else:
                        A[j,k] = 1
                elif j == i-1:
                    A[j,k] = -aM[k,j]
                else:
                    A[j,k] = 0
    else:
        A = rAction(iSeq[0:length-1], aM)*rAction([i], aM) # Recursively multiplies on the left
    return A

'''
This function checks if the actions of two sequences representing r_i^w and r_j^v are equal
Indexing begins at 1 not 0
Returns true or false
Input: iSeq, jSeq, aM (similar to rAction)
Output: result (true or false)
'''
def rEqual(iSeq, jSeq, aM):
    A1 = rAction(iSeq, aM) # Gets the two actions to compare
    A2 = rAction(jSeq, aM)
    result = (A1 == A2) # Compares them
    return result

'''
This function checks if r_i^w = r_i^v for all i
Indexing begins at 1 not 0
Returns true or false
Input: bM, cM, (exchange, coefficient matrices), w, v (the two mutation sequences in question), l (linear ordering of the indices of bM)
Output: result (true or false)
'''
def allREqual(bM, cM, w, v, l):
    aM = corrGIM(bM, l)
    n = aM.nrows() # Gets the number of indices
    result = true # default return is true
    for i in [1..n]:
        iSeqW = rMutation(i, bM, cM, w) # Calculates the two sequences of r_i's
        iSeqV = rMutation(i, bM, cM, v)
        if not rEqual(iSeqW, iSeqV, aM): # compares the two. If false they are not equal
            print('Linear order is not suitable')
            result = false
            break
    return result

'''
This function determines if the C-vectors are equal under matrix mMutation
Indexing begins at 1 not 0
Returns true or false
Input: bM, cM (exchange and c-vector matrices), w, v(mutation sequences len(w), len(v) >= 0)
Output: result (true or false)
'''
def cEqual(bM, cM, w, v):
    n = bM.nrows()
    M = block_matrix([(bM, cM)]) # Mutates the two together
    cW = mMutation(M, w)[:,n:] # Gets the C matrices
    cV = mMutation(M, v)[:,n:]
    result = (cW == cV) #compares the two
    return result

'''
Finally, this function checks to see if the conjecture is satisfied for two given mutation sequences
Indexing begins at 1 not 0
Returns true or false
Input: bM, cM, (exchange, coefficient, matrices), w, v (the two mutation sequences in question), l (total ordering of indices)
Output: result (true or false. Default is true)
'''
def conjSatisfied(bM, cM, w, v, l):
    n = bM.nrows()
    result = true
    if cEqual(bM, cM, w, v):
        result = allREqual(bM, cM, w, v, l)
    else:
        print("C^w does not equal C^v. Does not satisfy assumptions of conjecture")
        result = false
    return result




'''
It's important to remember that we only care about based loops in the exchange graph. Now we need a way to generate loops in the exchange graph. This is probably not an efficent way to do this. adj is an adjacency matrix, i is that starting seed/vertex. Returns a mutation sequence that gives a loop in the exchange graph. Indexing starts at 1 not 0.
'''

def randomWalk(oB):
    if oB.nrows() == oB.ncols():
        n = oB.ncols()
        B = oB.set_immutable() #immutable version of B if we need it
        G = fullExGraph(oB,Graph()) #creates the associhedron for B
        MarkovB = matrix(QQ,G.adjacency_matrix()) #A markov chain for the associhedron
        u = sum(posB.columns()) #this gives a vector of the row sums of our adjacency matrix
        for i in [1..n]:
            for j in [1..n]:
                MarkovB[i-1,j-1] = MarkovB[i-1,j-1]/u[i-1] #turns B into a markov chain matrix with uniform probabilities
        walk = [] #empty walk vector
        walk = rCancel(walk) #this removes any backtracking/extraneous steps
    else:
        print("B matrices must be square.")
    return walk




'''
##############################################################################################
USING B INSTEAD OF M GIVES THE DESIRED EXCHANGE GRAPH.
##############################################################################################
'''

'''
Adds the vertices and edges connected to a matrix in the "extended" exchange graph
Input: oM (n x 2n matrix), gO (graph that is being mutated)
Output: G (Graph with mutations added to it)
'''
def exGraphMutation(oM, gO):
    G = copy(gO)
    o = G.order()
    M = copy(oM)
    M.set_immutable()
    n = M.nrows()
    if o == 0:
        G.add_vertex(M)
    else:
        for i in [1..n]:
            mut = mMutation(M, [i])
            mut.set_immutable()
            edge = {M, mut}
            str = 'mu %i'%i
            if G.has_vertex(mut):
                if not G.has_edge(edge):
                    G.add_edge(edge, label=str)
                    #print(str)
            else: 
                G.add_vertex(mut)
                G.add_edge(edge, label=str)
                #print(str)
    return G

'''
Adds the vertices and edges connected to a exchange graph for all vertices currently in the graph
Input: M (n x 2n matrix), gO (graph)
Output: G (graph with mutations added)
'''

def exGraph(M, gO):
    G = copy(gO)
    if G.order() == 0:
        G = exGraphMutation(M,G)
    vertices = G.vertices()
    n = M.nrows()
    for vert in vertices:
        edges = G.neighbors(vert)
        if len(edges) != n:
            newG = exGraphMutation(vert, G)
            G = G.union(newG)
    return G

'''
Creates the full "extended" exchange graph for a given matrix and graph
Input: M( n x 2n matrix), gO (graph)
Output: G (extended exchange graph)
'''

def fullExGraph(M, g0):
    G = copy(g0)
    newG = exGraph(M, G)
    while G != newG:
        G = copy(newG)
        newG = exGraph(M, G)
    return G
    
'''
Add G.plot() in the Jupyter notebook to make it plot it.
Be warned, these are not the normal exchange graphs, but extended ones (Have not figured out how to make the plot bigger)
The standard orientation produces 84 vertices
By which I mean they are the graph formed by mutating [B C] until it pans out
'''
'''
Trying to find the full extended exchange graph of A_10 does not load
MAkes it not computationally feasible
But does the normal exchange graph give us what we want
'''



