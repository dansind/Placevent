'''
Creates a pickle containing shell indices up to a limit
Doesn't need to be run if you already made a pickle
'''
import pickle
def main():
    maxindex=50
    from math import sqrt
    shellindices=[[[0,0,0]]]
    for index in range(1,maxindex):
        #range[0]
        indicesinthisshell=[]
        for i in range(-index,index+1):
            for j in range(-index,index+1):
                for k in range(-index,index+1):
                    #print "math.sqrt(i*i + j*j + k*k))=",math.sqrt(i*i + j*j + k*k),"index = ",index
                    if int(sqrt(i*i + j*j + k*k)) == index : # I think this will miss some
                        indicesinthisshell.append([i,j,k])
        shellindices.append(indicesinthisshell)
    f=open("shells.pickle","wb")
    pickle.dump(shellindices,f)
    f.close()

def readshellindices():
    f=open("shells.pickle","rb")
    shellindices=pickle.load(f)
    return(shellindices)

if __name__ == '__main__' :
    main()
