'''Daniel J. Sindhikara
    sindhikara@gmail.com
    Does stuff related to dx files
automatically detects and reads gzipped dx files
    '''

def readdx(filenames):
    '''
    Reads one or more dx files into memory
    returns grid data and single 4D array containing dx data [index][xindex][yindex][zindex]
    '''
    import gzip
    import numpy as np
    def opendxfile(filename):
        if "gz" in filename:
            dxfile = gzip.open(filename,"rb")
        else:
            dxfile = open(filename,"r")
        return(dxfile)
        
    dxfile = opendxfile(filenames[0])
    dxlines = []
    for i in range(10): # only need the first few lines to get grid data
        dxlines.append(dxfile.readline())
    dxfile.close()
    gridcount = []
    origin = []
    deltas = [0,0,0]
    startline=0
    for i,line in enumerate(dxlines) :
        splitline = line.split()
        if len(splitline) > 2 :
            if(splitline[1]=="1") :
                gridcount.append(int(splitline[5]))
                gridcount.append(int(splitline[6]))
                gridcount.append(int(splitline[7]))
                #print "# gridcounts ",gridcount
            if(splitline[0]=="origin") :
                origin.append(float(splitline[1]))
                origin.append(float(splitline[2]))
                origin.append(float(splitline[3]))
            if(splitline[0]=="delta") :
                if(float(splitline[1])>0) :
                    deltas[0]=float(splitline[1])
                if(float(splitline[2])>0) :
                    deltas[1]=float(splitline[2])
                if(float(splitline[3])>0) :
                    deltas[2]=float(splitline[3])
            if(splitline[1]=="3") :
                numpoints=int(splitline[9])
                #print "# Total number of gridpoints is ",numpoints
                startline=i+1
        if(startline>1) :
            break
    distributions = [[[[0.0 for x in range(gridcount[2])] for y in range(gridcount[1])] for z in range(gridcount[0])]for w in range(len(filenames))]
    gridvolume=deltas[0]*deltas[1]*deltas[2]
    #print "# I have to read",len(filenames)*gridcount[2]*gridcount[1]*gridcount[0],"values."
    for i, dxfilename in enumerate(filenames): 
        dxfile = opendxfile(dxfilename)
        dxtext = dxfile.read()
        dxfile.close()
        splittext=dxtext.split()
        del splittext[0:splittext.index("follows")+1] # get rid of header text, last word is "follows"
        del splittext[-4:-1]
        del splittext[-1] # dont know why i cant just use a 0 instead of -1 above
        floats = [float(element) for element in splittext]
        index=0
        for x in range(gridcount[0]):
            for y in range(gridcount[1]):
                for z in range(gridcount[2]):
                    distributions[i][x][y][z]=floats[index]
                    index+=1 # there's probably a more pythonic way to do this
    #print "# corr[0][0][0][0] = ",distributions[0][0][0][0],"[0][-1][-1][-1]=",distributions[0][-1][-1][-1]
    #distributions=np.array(distributions)
    return(distributions,origin,deltas,gridcount)    


def printdxfrom1dzfast(values,origin,delta,gridcounts,filename):
    ''' Print a dx file'''
    f = open(filename,"w")
    f.write("#DX file from Dan's program\n")
    f.write("object 1 class gridpositions counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("origin {0} {1} {2}\n".format(origin[0],origin[1],origin[2]))
    f.write("delta {0} 0 0\n".format(delta[0]))
    f.write("delta 0 {0} 0\n".format(delta[1]))
    f.write("delta 0 0 {0}\n".format(delta[2]))
    f.write("object 2 class gridconnections counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("object 3 class array type double rank 0 items {0} data follows\n".format(gridcounts[0]*gridcounts[1]*gridcounts[2]))
    for value in values:
        f.write("{0}\n".format(value))
    f.write("object {0} class field\n".format(filename))
    f.close()



def printdx(values,indices,origin,delta,gridcounts,filename):
    ''' Print a dx file'''
    '''first construct full grid as 1d list (z fast, y med, z slow)
    note values do not span entire grid, so just replce them if they exist
    values = 1D array containing all values
    indices[value index][dimension] = index in the dimension specified (as opposed to 1D index)
    origin = length 3 float list of origin coordinate
    delta = length 3 float list of grid spacing
    gridcounts = length 3 integer list of number of grid points in each dimension
    filename = string of output file name
    '''
    gridvalues=[0.0]*gridcounts[0]*gridcounts[1]*gridcounts[2]
    if len(values) != len(gridvalues):
        exit("error! len(gridvalues) %s != len(values) %s" % (len(gridvalues),len(values)))
    for i,value in enumerate(values):  # warning, this doubles the memory usage
        gridindex=indices[i][2]+indices[i][1]*gridcounts[2]+indices[i][0]*gridcounts[2]*gridcounts[1]
        gridvalues[gridindex]=value

    f = open(filename,"w")
    f.write("#DX file from Dan's program\n")
    f.write("object 1 class gridpositions counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("origin {0} {1} {2}\n".format(origin[0],origin[1],origin[2]))
    f.write("delta {0} 0 0\n".format(delta[0]))
    f.write("delta 0 {0} 0\n".format(delta[1]))
    f.write("delta 0 0 {0}\n".format(delta[2]))
    f.write("object 2 class gridconnections counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("object 3 class array type double rank 0 items {0} data follows\n".format(gridcounts[0]*gridcounts[1]*gridcounts[2]))
    for gridvalue in gridvalues:
        f.write("{0}\n".format(gridvalue))
    f.write("object {0} class field\n".format(filename))
    f.close()

def printdxfrom3d(distribution,origin,delta,gridcounts,filename):
    ''' print a dx file given a 3d list'''

    f = open(filename,"w")
    f.write("#DX file from Dan's program\n")
    f.write("object 1 class gridpositions counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("origin {0} {1} {2}\n".format(origin[0],origin[1],origin[2]))
    f.write("delta {0} 0 0\n".format(delta[0]))
    f.write("delta 0 {0} 0\n".format(delta[1]))
    f.write("delta 0 0 {0}\n".format(delta[2]))
    f.write("object 2 class gridconnections counts {0} {1} {2}\n".format(gridcounts[0],gridcounts[1],gridcounts[2]))
    f.write("object 3 class array type double rank 0 items {0} data follows\n".format(gridcounts[0]*gridcounts[1]*gridcounts[2]))
    for i in range(gridcounts[0]):
        for j in range(gridcounts[1]):
            for k in range(gridcounts[2]):
                f.write("{0}\n".format(distribution[i][j][k]))
    f.write("object {0} class field\n".format(filename))
    f.close()


def getcoordfromindices(indices,origin,deltas):
    '''Returns coordinates as a length 3 list of floats 
    '''
    coords=[]
    for i in range(3) :
        coords.append(float(indices[i])*deltas[i]+origin[i])
    return(coords)

def getindicesfromcoord(coord,origin,deltas):
    indices=[]
    for i in range(3):
        indices.append(int((coord[i]-origin[i])/deltas[i]+0.5))
    return indices


def precomputeshellindices(maxindex):
    '''return a list of 3d lists containing the indices in successive search shells
i.e. 0 = [0,0,0]
     1 = [[-1,-1,-1],[-1,0,-1],..
    essentially how much to shift i,j,k indices from center to find grid point on shell
at index radius
    This will make evacuation phase faster
    '''
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
    return(shellindices)

def linearinterpolatevalue(distribution,origin,deltas,coord):
    '''given a 3d coordinate, using a linear interpolation from the 8 nearest gridpoints,
estimate the value at that coordinate
    '''
    ccrd=[]# coordinates of corners
    cindex=[]# indices of corners
    cdist=[]# distances to corner
    #below store indices and coordinates of 8 nearby corners
    cindex.append([int((coord[0]-origin[0])/deltas[0]),int((coord[1]-origin[1])/deltas[1]),int((coord[2]-origin[2])/deltas[0])])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas)) 
    cindex.append([cindex[0][0]+1,cindex[0][1],cindex[0][2]])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    cindex.append([cindex[0][0],cindex[0][1]+1,cindex[0][2]])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    cindex.append([cindex[0][0],cindex[0][1],cindex[0][2]+1])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    cindex.append([cindex[0][0]+1,cindex[0][1]+1,cindex[0][2]])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    cindex.append([cindex[0][0],cindex[0][1],cindex[0][2]+1])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    cindex.append([cindex[0][0]+1,cindex[0][1],cindex[0][2]+1])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    cindex.append([cindex[0][0]+1,cindex[0][1]+1,cindex[0][2]+1])
    ccrd.append(getcoordfromindices(cindex[-1],origin,deltas))
    totalweight=0.0
    weights=[]
    for crd in ccrd: # coordinates of corners
        myweight=((coord[0]-crd[0])**2+(coord[1]-crd[1])**2+(coord[2]-crd[2])**2)**(-0.5)
        weights.append(myweight)
        totalweight+=myweight
    value=0.0
    for i,mycindex in enumerate(cindex):
        try:
            value+=distribution[mycindex[0]][mycindex[1]][mycindex[2]]*weights[i]/totalweight
        except:
            print "Failed to find gridpoint at",mycindex[0],mycindex[1],mycindex[2]
            print "coordinate=",coord
    return(value)


def calcrdf(distribution,origin,deltas,coord,deltar=0.1,maxradius=20.0,numsumgrids=20):
    '''
Calculates the radial distribution function about a point using the 3d distribution
    '''
    def shellintegral(radius,delta,center):
        sum=0.0
        count=0.0
        for i in range(int(2.0*radius/delta)):
            for j in range(int(2.0*radius/delta)):
                for k in range(int(2.0*radius/delta)):
                    x=float(i)*delta-radius
                    y=float(j)*delta-radius
                    z=float(k)*delta-radius
                    thisrad=(x**2+y**2+z**2)**(0.5)
                    if thisrad > radius - delta/2.0 and thisrad < radius + delta/2.0 : 
                        mycoord=[center[0]+x,center[1]+y,center[2]+z]
                        sum+=linearinterpolatevalue(distribution,origin,deltas,mycoord)
                        count+=1.0
        sum=sum/count
        return(sum) 
    radii=[(float(i)+0.5)*deltar for i in range(int(maxradius/deltar))]
    gr=[]
    for i,rad in enumerate(radii) :
        subdelta=deltar*(float(i)+0.5)/float(numsumgrids)
        gr.append(shellintegral(rad,subdelta,coord))

    return(radii,gr) 
   





