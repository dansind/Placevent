#!/usr/bin/env python
'''Created by Daniel Sindhikara, sindhikara@gmail.com
Placevent
This program is designed to automatically place explicit solvent atoms/ions based
on 3D-RISM data. The 3D-RISM correlations should be in a DX file

The details of the algorithm are described in: 
Placevent: An algorithm for prediction of explicit solvent atom distribution -- Awpplication to HIV-1 protease and F-ATP synthase 
Daniel J. Sindhikara, Norio Yoshida, Fumio Hirata
http://dansindhikara.com/Software/Entries/2012/6/22_Placevent_New.html

This package requires Numpy.

There is a degree of uncertainty in the placement due to the conversion of continuous distribution to
a discrete one. The atoms are printed in PDB format with the beta term reflecting the uncertainty.

printing g(r) in beta spot

//    Copyright (C) 2012 Daniel J. Sindhikara
//
//    This program is free software: you can redistribute it and/or modify
//    it under the terms of the GNU General Public License as published by
//    the Free Software Foundation, either version 3 of the License, or
//    (at your option) any later version.
//
//    This program is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//    GNU General Public License for more details.
//
//    You should have received a copy of the GNU General Public License
//    along with this program.  If not, see <http://www.gnu.org/licenses/>.



TODO: 


'''
from math import *
import sys
import numpy as np

import modules.djsgrid as grid
import modules.djspdb as pdb 


def converttopop(distribution,delta,conc):
    ''' Convert distribution function to population function
    output numpy array, total population as float
    '''
    xlen=len(distribution)
    ylen=len(distribution[0])
    zlen=len(distribution[0][0])
    popzero = [[[0.0 for z in range(zlen)] for y in range(ylen)] for x in range(xlen)]
    gridvolume=delta[0]*delta[1]*delta[2]
    print "# volume of each grid point is ",gridvolume," cubic angstroms"
    # z fast, y med, x slow
    #     #convert to population array
    totalpop=0;
    for i in range(xlen) :
        for j in range(ylen) :
            for k in range(zlen) :
                 popzero[i][j][k]=distribution[i][j][k]*conc*gridvolume
                 totalpop+=popzero[i][j][k]
    print "# There are approximately ",totalpop," sites within the 3D-RISM grid"
    return(popzero,totalpop,gridvolume)   

def doplacement(popzero,conc,gridvolume,origin,delta,shellindices,grcutoff):
    ''' Does actual placement, returns array of "Atom" objects to be printed later 
    '''
    popi=np.array(popzero)# popi will change, popzero will stay constant
    max = 0
    topindices = [[[]]]
    #print "# printing g(r) in beta column"
    #print "# printing g(r)_i in volume column"
    placedcenters=[]
    finished=0
    print "# Doing placement..."
    while finished==0 : 
        mymax=popi.max()
        maxi,maxj,maxk=np.argwhere(popi==mymax)[0]
        topindices.append([maxi,maxj,maxk])
        # need to subtract out one populations worth
        # if using shell index list, replace below
        remainder = 1 # population of one
        indexradius = 0 # first shell radius
        while remainder > 0 :
            availablepopulation=0
            try: 
                for indices in shellindices[indexradius]: # search only through indices lying on shell
                    try :
                        availablepopulation+=popi[maxi+indices[0],maxj+indices[1],maxk+indices[2]]
                    except :
                        print "# Stopping! Population search went outside grid. You may want to try again with a higher concentration"
                        remainder = -1 # cue quitting search   
                        finished=1
                        break
            except: 
                 print "# Stopping! Population search went outside maximum shell size. This usually happens with a very dilute concentration.  You may want to try again with a higher concentration"
                 remainder = -1 # cue quitting search   
                 finished=1
                 break
            if availablepopulation>remainder and remainder > -1 :
                #There is enough in this shell to bring the remainder to zero
                fraction=remainder/availablepopulation
                remainder=0
                for indices in shellindices[indexradius]:
                    popi[maxi+indices[0]][maxj+indices[1]][maxk+indices[2]]=popi[maxi+indices[0]][maxj+indices[1]][maxk+indices[2]]*fraction
            elif remainder > -1:
                # not enough in this shell, set it all to zero
                for indices in shellindices[indexradius]:
                    popi[maxi+indices[0]][maxj+indices[1]][maxk+indices[2]]=0
                remainder-=availablepopulation
                indexradius+=1
                if indexradius > len(shellindices) : # too big
                    print "# Stopping! Population search went outside grid, you may want to try again with a higher concentration"
                    exit()

        myg0=popzero[maxi][maxj][maxk]/(conc*gridvolume) #actual original 
        mygi=mymax/(conc*gridvolume) #remaining correlation
        #betaterm=float(indexradius)*float(delta[0])# Really hope all deltas are the same
        #betaterm=float(fiftypercentradius)*float(delta[0])
        placedcenters.append(pdb.Atom(serial=len(placedcenters),x=maxi*delta[0]+origin[0],y=maxj*delta[1]+origin[1],z=maxk*delta[2]+origin[2],occ=myg0,tfac=mygi))
        #above can be modded to use grid
        if(myg0<grcutoff) :
            finished=1 # It's dilute enough already if the highest is less than g(r)=1.5
    return(placedcenters)

def main():
    print "# input format placevent.py <dxfile> <concentration M (molar)> [cutoff g(r) (default 1.5)]"
    print "# concentration (#/A^3) ~= molarity*6.022E-4"
    if(len(sys.argv)<3) :
    	print "Insufficient arguments need 2 : ",len(sys.argv)-1
    	exit()
    if len(sys.argv)==4 :
        grcutoff=float(sys.argv[3]) #setting cutoff
    else :
        grcutoff=1.5
    dxfilename=[sys.argv[1]]
    molar=float(sys.argv[2])
    print "# your dx file is",dxfilename,"molarity is",molar,"M"
    conc=molar*6.0221415E-4
    distributions,origin,delta,gridcount=grid.readdx(dxfilename)
    if molar > 10 :  # it's unlikely we'll need more than 20 shells
        numshells = 20
    else :
        numshells = 50
    print "# precalculating indices lying on",numshells," concentric shells"
    shellindices=grid.precomputeshellindices(numshells) 
    popzero,totalpop,gridvolume=converttopop(distributions[0],delta,conc)
    placedcenters=doplacement(popzero,conc,gridvolume,origin,delta,shellindices,grcutoff)
    print len(placedcenters)
    for i,center in enumerate(placedcenters):
        print center.makestring(), 
        #mystring=makepdbstring(i+1,"A","BCD"," ",1,center.x,center.y,center.z,center.gi,center.g0,"A",0.0)
        #print mystring[:-1]

if __name__ == '__main__' :
    main()
