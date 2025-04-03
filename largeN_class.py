import numpy as np

import fillH_largeN as fillH
import sys
import time
imag = complex(0, 1)
#norbitals is the number of sites in each unit cell


class largeN_lattice:
    def __init__(self, prim_vec, orbitals, n_of_sites):
        self.prim_vec = np.array(prim_vec)
        self.ansatz=0;
        self.orbitals = np.array(orbitals)
        self.n_of_sites = np.array(n_of_sites)
        self.bvectors = []  # the primitive vector of the reciprocal lattice
        self.norbital = self.orbitals.shape[0]
        self.nbonds = 0 # keeps track the number of added bonds
        self.bondlist = [] # keeps track of the info of added bonds
        self.ml=[] # don't know yet
        self.nl=[] # don't know yet
        self.phaselist = [] # keeps track of phase factor related to the added bonds
        self.nsites = np.prod(n_of_sites) # the number of unit cells in the system
        self.jlist = [] # the heisenberg couplings
        self.mulist = [] # ux,uy->these become our chis and etas
        self.temperature=1.0
        if self.n_of_sites.size == 1:
            n_1 = self.n_of_sites[0]
            a_1_vec = self.prim_vec[0]
            b_1_vec = 2. * np.pi / a_1_vec
            self.bvectors.append(b_1_vec)
            unit_k_1_vec = b_1_vec / n_1
            lambda_1_list=np.array(range(0,n1))
            self.k_list = np.empty((n_1, 1))
            for i in range(n_1):
                self.k_list[i][0] = unit_k_1_vec[0] * lambda_1_list[i]
        elif self.n_of_sites.size == 2:
            n_1 = self.n_of_sites[0]
            n_2 = self.n_of_sites[1]
            a_1_vec = self.prim_vec[0]
            a_2_vec = self.prim_vec[1]
            eta = 2 * np.pi / (np.cross(a_1_vec, a_2_vec))
            b_1_vec = np.array([eta * a_2_vec[1], -eta * a_2_vec[0]])
            b_2_vec = np.array([-eta * a_1_vec[1], eta * a_1_vec[0]])
            self.bvectors.append(b_1_vec)
            self.bvectors.append(b_2_vec)
            unit_k_1_vec = b_1_vec / n_1
            unit_k_2_vec = b_2_vec / n_2
            lambda_1_list=np.array(range(0,n_1))
            lambda_2_list=np.array(range(0,n_2))
            self.k_list = np.empty((n_1, n_2, 2))
            for i in range(n_1):
                for j in range(n_2):
                    self.k_list[i][j] = lambda_1_list[i] * unit_k_1_vec + lambda_2_list[j] * unit_k_2_vec
            self.k_list = np.reshape(self.k_list, (n_1 * n_2, 2))
        elif self.n_of_sites.size == 3:
            n_1 = self.n_of_sites[0]
            n_2 = self.n_of_sites[1]
            n_3 = self.n_of_sites[2]
            a_1_vec = self.prim_vec[0]
            a_2_vec = self.prim_vec[1]
            a_3_vec = self.prim_vec[2]
            b_1_vec = 2 * np.pi * np.cross(a_2_vec, a_3_vec) / (np.dot(a_1_vec, np.cross(a_2_vec, a_3_vec)))
            b_2_vec = 2 * np.pi * np.cross(a_3_vec, a_1_vec) / (np.dot(a_1_vec, np.cross(a_2_vec, a_3_vec)))
            b_3_vec = 2 * np.pi * np.cross(a_1_vec, a_2_vec) / (np.dot(a_1_vec, np.cross(a_2_vec, a_3_vec)))
            self.bvectors.append(b_1_vec)
            self.bvectors.append(b_2_vec)
            self.bvectors.append(b_3_vec)
            unit_k_1_vec = b_1_vec / n_1
            unit_k_2_vec = b_2_vec / n_2
            unit_k_3_vec = b_3_vec / n_3
            lambda_1_list=np.array(range(0,n_1))
            lambda_2_list=np.array(range(0,n_2))
            lambda_3_list=np.array(range(0,n_3))
            self.k_list = np.empty((n_1, n_2, n_3, 3))
            for i in range(n_1):
                for j in range(n_2):
                    for l in range(n_3):
                        self.k_list[i][j][l] = (lambda_1_list[i] * unit_k_1_vec + lambda_2_list[j] * unit_k_2_vec +
                                                lambda_3_list[l] * unit_k_3_vec)
            self.k_list = np.reshape(self.k_list, (n_1 * n_2 * n_3, 3))
        else:
            print('Feature not yet supported.')
        self.bvectors = np.array(self.bvectors)
    def add_ansatz(self, n_icouplings, transform_matrix) :
        self.ansatz=transform_matrix

    def add_bond(self, m, n, distance):
        self.bondlist.append([[m, n], distance])
        self.ml.append(m)
        self.nl.append(n)
        self.phaselist.append(np.exp(imag * np.dot(self.k_list, distance)))
        self.nbonds += 1

    def initH(self): # organise the bond information
        self.phaselist = np.array(self.phaselist, dtype=np.cdouble)
        self.repeated_phase = np.repeat(self.phaselist, self.norbital, axis=1)
        self.ml=np.array(self.ml,dtype=np.intc)
        self.nl=np.array(self.nl,dtype=np.intc)
        self.H=np.zeros([self.nsites,self.norbital,self.norbital],dtype=np.cdouble)

    def fixed_couplings(self, jlist,filling=1): # the fixed input of the system
        if (len(jlist) == self.nbonds):
            self.jlist = np.array(jlist)
            self.filling = filling
        else:
            # raise ValueError('Inconsistent Input')
            print(jlist)
            print("Invalid input.")



    def make_hamiltonian(self, coupling): 
        #print(coupling.shape,self.nbonds+self.norbital)
        if (coupling.shape[0] ==self.norbital):
            #print("here",self.nbonds,coupling.shape)
            coupling = np.array(coupling, dtype=np.double)
            self.mulist=np.array(coupling)
            #print("fc=",self.H.flags.f_contiguous)
            #print("cc=",self.H.flags.c_contiguous)
            fillH.fillH(self.norbital, self.H, 1/self.temperature,self.jlist,self.mulist, self.nbonds,
                                self.ml, self.nl, self.phaselist, self.nsites)
        else:
            print('Invalid input, incorrect number of bond/site amplitudes')


    def calc_exps(self, icoupling):
        transform_matrix = self.ansatz        
        couplings = np.dot(transform_matrix, icoupling)
        self.make_hamiltonian(couplings)
        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        energies = energies.flatten()
        #print("emin=",np.min(energies),"np.max=",np.max(energies));
        eigenvec = eigenvec.swapaxes(1, 2).reshape(self.nsites*self.norbital, self.norbital).T
        inve=1/energies;

        occs=np.zeros([self.norbital],dtype=complex)
        for i in range(self.norbital) :
            occs[i]=np.dot(eigenvec[i].conj() *inve ,eigenvec[i].T)/self.nsites

        return occs
    def calc_ev(self, icoupling):
        transform_matrix = self.ansatz        
        couplings = np.dot(transform_matrix, icoupling)
        self.make_hamiltonian(couplings)
        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        return energies,eigenvec



    def print_corrfuncs(self, icoupling):
        cfuncs=np.zeros([self.nsites,self.norbital,self.norbital],dtype=complex)

        transform_matrix = self.ansatz
        couplings = np.dot(transform_matrix, icoupling)
        #couplings is [delta,bonds+sites]
        self.make_hamiltonian(couplings)
        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        inve=1/energies;
        for i in range(self.norbital) :
            for j in range(self.norbital) :
                for mom in range(self.nsites) :
                    ek=energies[mom,:]
                    inve=1.0/ek
                    cfuncs[mom,i,j]= np.dot( np.conjugate(eigenvec[mom, i, :] *inve), eigenvec[mom,j,:])


        np.save('corrfuncs%d.npy'%(self.n_of_sites[0]),cfuncs);




