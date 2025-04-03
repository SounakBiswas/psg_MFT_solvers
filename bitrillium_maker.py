import sys
from psg_hamiltonian_class import psg_lattice
from largeN_class import largeN_lattice
import numpy as np
def make_largeN_lattice(L) :
    threedim_bravais = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    u=0.3354;
    v=0.5952
    unit_cell=[[u,u,u], [0.5+u,0.5-u,1-u], [1-u,0.5+u,0.5-u], [0.5-u,1-u,0.5+u],
              [v,v,v], [0.5+v,0.5-v,1-v], [1-v,0.5+v,0.5-v], [0.5-v,1-v,0.5+v]]
    unit_cell=np.array(unit_cell)
    bitrillium = largeN_lattice(threedim_bravais,unit_cell, [L, L, L])

    bonds=np.loadtxt("bitrillium_bonds5nn");
    print("test")
    ep=0.000001
    for i in range(bonds.shape[0]) :
        bitrillium.add_bond(int(bonds[i,0]+ep),int(bonds[i,1]+ep),[bonds[i,2], bonds[i,3],bonds[i,4]])
    jlist=np.zeros(52);
    #j1=0.4; j2 =-0.15; j3=1; j4=5.4; j5=2.54
    #j1=0.364; j2 =-0.144; j3=0.798; j4=5.545; j5=2.657
    #j1=0.0; j2 =0.00; j3=1; j4=0.0; j5=2.54
    j1=0.0; j2 =0.00; j3=0; j4=1.0; j5=1

    jlist[0:4]=j1
    jlist[4:16]=j2
    jlist[16:28]=j3
    jlist[28:40]=j4
    jlist[40:52]=j5
    
    bitrillium.fixed_couplings(jlist)

    return bitrillium
   


     

def make_psg_lattice(L) :
    threedim_bravais = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    u=0.3354;
    v=0.5952
    unit_cell=[[u,u,u], [0.5+u,0.5-u,1-u], [1-u,0.5+u,0.5-u], [0.5-u,1-u,0.5+u],
              [v,v,v], [0.5+v,0.5-v,1-v], [1-v,0.5+v,0.5-v], [0.5-v,1-v,0.5+v]]
    unit_cell=np.array(unit_cell)
    bitrillium = psg_lattice(threedim_bravais,unit_cell, [L, L, L])

    bonds=np.loadtxt("bitrillium_bonds5nn");
    ep=0.000001
    for i in range(bonds.shape[0]) :
        bitrillium.add_bond(int(bonds[i,0]+ep),int(bonds[i,1]+ep),[bonds[i,2], bonds[i,3],bonds[i,4]])
    jlist=np.zeros(52);
    #j1=0.4; j2 =-0.15; j3=1; j4=5.4; j5=2.54
    j1=0.0; j2 =0.00; j3=1; j4=5.4; j5=2.54

    jlist[0:4]=j1
    jlist[4:16]=j2
    jlist[16:28]=j3
    jlist[28:40]=j4
    jlist[40:52]=j5
    
    bitrillium.fixed_couplings(jlist)

    return bitrillium
   


     

