import sys
from psg_hamiltonian_class import psg_lattice
from largeN_class import largeN_lattice
import numpy as np
def make_largeN_lattice(L) :
    threedim_bravais = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    u=0.138;
    unit_cell=[[u,u,u], [0.5+u,0.5-u,1-u], [1-u,0.5+u,0.5-u], [0.5-u,1-u,0.5+u]]
    
    unit_cell=np.array(unit_cell)
    trillium = largeN_lattice(threedim_bravais,unit_cell, [L, L, L])
    print(trillium.norbital)

    bonds=np.loadtxt("trill_bondsnn");
    ep=0.000001
    print(bonds.shape)
    for i in range(bonds.shape[0]) :
        trillium.add_bond(int(bonds[i,0]+ep),int(bonds[i,1]+ep),[bonds[i,2], bonds[i,3],bonds[i,4]])
    jlist=np.ones(12);
    j1=1.0


    trillium.fixed_couplings(jlist)

    return trillium
   


     

def make_psg_lattice(L) :
    threedim_bravais = [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]]
    u=0.25
    unit_cell=[[u,u,u], [0.5+u,0.5-u,1-u], [1-u,0.5+u,0.5-u], [0.5-u,1-u,0.5+u]]
    unit_cell=np.array(unit_cell)
    trillium = psg_lattice(threedim_bravais,unit_cell, [L, L, L])

    bonds=np.loadtxt("trill_bondsnn");
    ep=0.000001
    print(bonds.shape)
    for i in range(bonds.shape[0]) :
        trillium.add_bond(int(bonds[i,0]+ep),int(bonds[i,1]+ep),[bonds[i,2], bonds[i,3],bonds[i,4]])
    jlist=np.ones(12);
    
    trillium.fixed_couplings(jlist)

    return trillium
   


     

