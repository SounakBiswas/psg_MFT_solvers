import sys
from psg_hamiltonian_class import psg_lattice
#from class_backup import psg_lattice
import numpy as np
   


     

def make_psg_lattice(L,ji,jo,jint) :
    ht=np.sqrt(3)/2.0
    #bravais = [[1., 0.], [0.0, 1]]
    bravais = [[1., 0.], [0.5, ht]]
    bravais=np.array(bravais)
    sh=0.5;
    unit_cell=[[0.,0.]]*6
    print(unit_cell)
    
    
    bikag = psg_lattice(bravais,unit_cell, [L, L])
    bonds=[]
    #black kagome
    bonds= bonds+[[0,1,0,0],[1,2,0,0],[2,0,0,0]]
    bonds= bonds+[[2,1,1,0],[2,0,0,1],[1,0,-1,1]]
    #red kagome
    bonds= bonds+[[3,4,0,0],[4,5,0,0],[5,3,0,0]]
    bonds= bonds+[[4,3,1,0],[4,5,1,-1],[5,3,0,1]]
    #bonds= bonds+[[3,4,1,0],[4,5,1,-1],[3,5,0,1]]
    #inter 
    bonds= bonds+[[0,3,0,0],[0,4,0,0],[4,2,0,0]]
    bonds= bonds+[[2,5,0,0],[5,1,0,0],[1,3,0,0]]
    bonds=np.array(bonds)

    ep=0.000001
    for i in range(bonds.shape[0]) :
        dist= bonds[i,2]*bravais[0]+bonds[i,3]*bravais[1]
        bikag.add_bond(int(bonds[i,0]+ep),int(bonds[i,1]+ep),dist)
    
    #j1=0.4; j2 =-0.15; j3=1; j4=5.4; j5=2.54
    jlist=np.zeros(18)

    jlist[0:3]=ji
    jlist[3:6]=jo
    jlist[6:9]=ji
    jlist[9:12]=jo
    jlist[12:18]=jint
    
    bikag.fixed_couplings(jlist)

    return bikag
   


     

