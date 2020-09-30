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
Input: M (matrix), w (mutation sequence)
Output: mut (mutated Matrix)
'''
def mMutation(M, w):
    r = M.nrows()
    c = M.ncols()
    for i in w:
        if i < 1 or i > r or i > c: # Checks to see if we are mutating at non-vertices
            print("Invalid mutation at vertex: ", i)
            return M
    mut = matrix(r,c,0) # Builds a 0 matrix because letting mut = M runs into problems
    k = w[0]
    if len(w) == 1:
         for i in range(r):
             for j in range(c): # Standard matrix mutation here
                 if i == k-1 or j == k-1: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = -M[i,j]
                 elif M[i,k-1]*M[k-1,j] > 0: # Has the -1's here to account for sage indexing at 0
                     mut[i,j] = M[i,j] + M[i,k-1]*M[k-1,j]*(M[i,k-1]/abs(M[i,k-1])) # Has the -1's here to account for sage indexing at 0
                 else:
                     mut[i,j] = M[i,j]
    else:
        mut = mMutation(mMutation(M,w[0:len(w)-1]), [w[len(w)-1]]) # Recursively "goes up" the mutation sequence
    return mut

'''
This function takes in an index i, the initial exchange matrix, the initial c-vector matrix, and a mutation sequence w
It then returns a sequence representing the indices of the r_i's
Indexing of the matrices begins at 1 not 0, as we are accustomed to
Both matrices should be square
Inputs: i (index), bM (initial exchange matrix), cM(initial c-vector matrix), w (mutation sequence)
Output: rW (index sequence representing r_i^w)
'''
def rMutation(i, bM, cM, w):
     wPrime = []
     length = len(w)-1
     rW = []
     if i < 1 or i > bM.nrows(): # Checks to see if this is a valid index we are picking
         print ("Invalid mutation at index: ", i)
         return w
     if length == 0: # This is when w is a single mutation
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
