import numpy as np
import fillH
import sys
import time
imag = complex(0, 1)
#norbitals is the number of sites in each unit cell


class psg_lattice:
    def __init__(self, prim_vec, orbitals, n_of_sites):
        self. current_ansatz=-1
        self.ansatzen = []
        self.n_ansatzen=0
        self.prim_vec = np.array(prim_vec)
        self.orbitals = np.array(orbitals)
        self.n_of_sites = np.array(n_of_sites)
        self.bvectors = []  # the primitive vector of the reciprocal lattice
        self.norbital = self.orbitals.shape[0]
        self.nband = 2 * self.norbital
        self.nbonds = 0 # keeps track the number of added bonds
        self.bondlist = [] # keeps track of the info of added bonds
        self.ml=[] # don't know yet
        self.nl=[] # don't know yet
        self.phaselist = [] # keeps track of phase factor related to the added bonds
        self.nsites = np.prod(n_of_sites) # the number of unit cells in the system
        self.jlist = [] # the heisenberg couplings
        self.ulist = [] # ux,uy->these become our chis and etas
        self.Llist= []# lx,ly
        self.filling=1
        #print(self.nbonds,self. norbital)
        self.n_couplings =  2*self.norbital
        self.temperature=0.00001
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
    def draw_BZ(self):  # paint the BZ given the bravais lattice structure
        pass

    def draw_unit_cell(self):  # paint the unit cell of the lattice
        pass

    def add_bond(self, m, n, distance):
        self.bondlist.append([[m, n], distance])
        self.ml.append(m)
        self.nl.append(n)
        #self.phaselist.append(np.exp(imag * np.dot(self.k_list, distance+self.orbitals[n]-self.orbitals[m])))
        self.phaselist.append(np.exp(imag * np.dot(self.k_list, distance)))
        self.nbonds += 1
        self.n_couplings+=2

    def add_ansatz(self, n_icouplings, transform_matrix) :
        print(self.n_couplings,n_icouplings)
        #print(transform_matrix.shape)
        #if(transform_matrix.shape == (self.n_couplings, n_icouplings)) :
        (self.ansatzen).append ([n_icouplings, transform_matrix])
            #print(self.ansatzen)
        #else :
        #    print("wrong parameters.")
    def initH(self): # organise the bond information
        self.phaselist = np.array(self.phaselist, dtype=np.cdouble)
        self.repeated_phase = np.repeat(self.phaselist, self.nband, axis=1)
        self.ml=np.array(self.ml,dtype=np.intc)
        self.nl=np.array(self.nl,dtype=np.intc)
        self.H=np.zeros([self.nsites,self.nband,self.nband],dtype=np.cdouble)
        #print("nbonds=",self.nbonds)

    def fixed_couplings(self, jlist,filling=1): # the fixed input of the system
        if (len(jlist) == self.nbonds):
            self.jlist = np.array(jlist)
            self.filling = filling
            print(jlist)
        else:
            # raise ValueError('Inconsistent Input')
            print(jlist,)
            print("Invalid input.")



    def make_hamiltonian(self, coupling): 
        #print(coupling.shape,self.nbonds+self.norbital)
        if (coupling.shape[0] == 4 and coupling.shape[1]==(self.nbonds + self.norbital)):
            #print("here",self.nbonds,coupling.shape)
            coupling = np.array(coupling, dtype=np.double)
            u0list= coupling[0,0:self.nbonds]
            uxlist= coupling[1,0:self.nbonds]
            uylist= coupling[2,0:self.nbonds]
            uzlist= coupling[3,0:self.nbonds]
            muxlist= coupling[1,self.nbonds:self.nbonds+self.norbital]
            muylist= coupling[2,self.nbonds: self.nbonds+self.norbital]
            muzlist= coupling[3,self.nbonds: self.nbonds+self.norbital]

            #fillH.fillH(self.H, u0list,  uxlist, uylist, uzlist, muxlist, muylist, muzlist, self.nbonds,
            #                    self.ml, self.nl, self.phaselist, self.nsites, self.nband,self.jlist)
            fillH.fillH(self.H, u0list,  uxlist, uylist, uzlist, muxlist, muylist, muzlist, self.nbonds,
                                self.ml, self.nl, self.phaselist, self.nsites, self.nband)
        else:
            print('Invalid input, incorrect number of bond/site amplitudes')

    def fermi_function(self, energy, fermi_energy):
        if self.temperature < 0.001:
            return np.heaviside(fermi_energy-energy, 0.5)
        else:
            return 1. / (np.exp((energy - fermi_energy) / self.temperature) + 1.)
    def calc_energy(self, icoupling):
        if(self.current_ansatz==-1) :
            raise Exception("set the ansatz, you cunt")
        #print(self.ansatzen)
        cost_arr =np.zeros(7,dtype=complex)

        transform_matrix = self.ansatzen[self. current_ansatz][1]
        n_icouplings = self.ansatzen[self. current_ansatz][0]
        couplings = np.dot(transform_matrix, icoupling)
        #couplings is [delta,bonds+sites]
        self.make_hamiltonian(couplings)
        #print("norm=",np.linalg.norm(self.H[30]-np.conjugate(np.transpose(self.H[30]))))

        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        flattened_energies = energies.flatten()
        flattened_eigenvec = eigenvec.swapaxes(1, 2).reshape(self.nsites*self.nband, self.nband).T
        fermi_energy = np.sort(energies, axis = None)[(int(self.filling+self.norbital)//2)*self.nsites - 1]
        fermi_energy=0
        fermi_dist = self.fermi_function(flattened_energies, fermi_energy)
        return np.dot(fermi_dist,flattened_energies)/(self.nsites)

    def calc_exps(self, icoupling):
        if(self.current_ansatz==-1) :
            raise Exception("set the ansatz, you cunt")
        #print(self.ansatzen)
        cost_arr =np.zeros(7,dtype=complex)

        transform_matrix = self.ansatzen[self. current_ansatz][1]
        n_icouplings = self.ansatzen[self. current_ansatz][0]
        couplings = np.dot(transform_matrix, icoupling)
        #couplings is [delta,bonds+sites]
        self.make_hamiltonian(couplings)
        #print("norm=",np.linalg.norm(self.H[30]-np.conjugate(np.transpose(self.H[30]))))

        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        flattened_energies = energies.flatten()
        flattened_eigenvec = eigenvec.swapaxes(1, 2).reshape(self.nsites*self.nband, self.nband).T
        fermi_energy = np.sort(energies, axis = None)[(int(self.filling+self.norbital)//2)*self.nsites - 1]
        #print("actual fermi energy=",fermi_energy)
        fermi_energy=0
        fermi_dist = self.fermi_function(flattened_energies, fermi_energy)


        #for general su2 symmetric case, we might need all of tx,ty,tz. We have tz appearing in TR, so
        # we have two real couplings per bond and site
        bond_exps=np.zeros([self.nbonds,4],dtype=complex) # four elements of the matrix
        site_exps=np.zeros([self.nsites,4],dtype=complex)
        for i in range(self.nbonds):

            m_u = self.ml[i] #m_up
            n_u = self.nl[i] #n_up
            m_dd = self.norbital+self.ml[i] #m_down_dagger
            n_dd = self.norbital+self.nl[i] #n_down_dagger

            hop1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_u].T))/self.nsites
            hop2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_dd].T))/self.nsites


            pair1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_dd].T))/self.nsites


            pair2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_u].T))/self.nsites
            bond_exps[i,0]=1j*(hop1-hop2)
            bond_exps[i,3]=hop1+hop2
            bond_exps[i,1]=pair1-pair2
            bond_exps[i,2]=-1j*(pair1+pair2)


            #cost funcs later




        tot=0;
        for m in range(self.norbital):
            m_u=m
            m_dd=self.norbital+m

            hop1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_u].T))/self.nsites

            hop2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_dd].T))/self.nsites


            pair1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_dd].T))/self.nsites

            pair2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_u].T))/self.nsites
            #print("m=",m,"print exps",hop1,hop2)
            #print(fermi_dist[self.nsites*self.nband//2],fermi_dist[self.nsites*self.nband//2+1],fermi_dist[self.nsites*self.nband//2-1])
            #print(flattened_energies[self.nsites*self.nband//2],flattened_energies[self.nsites*self.nband//2+1],flattened_energies[self.nsites*self.nband//2-1])
            #sys.stdin.read(1)

            site_exps[m,0]=(hop1-hop2)
            site_exps[m,3]=hop1+hop2
            site_exps[m,1]=pair1-pair2
            site_exps[m,2]=1j*(pair1+pair2)
            tot+=(hop1-hop2);



        return bond_exps,site_exps





    def print_expvals(self, icoupling):
        if(self.current_ansatz==-1) :
            raise Exception("set the ansatz, you cunt")
        #print(self.ansatzen)

        transform_matrix = self.ansatzen[self. current_ansatz][1]
        n_icouplings = self.ansatzen[self. current_ansatz][0]
        couplings = np.dot(transform_matrix, icoupling)
        #couplings is [delta,bonds+sites]
        self.make_hamiltonian(couplings)

        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        flattened_energies = energies.flatten()
        flattened_eigenvec = eigenvec.swapaxes(1, 2).reshape(self.nsites*self.nband, self.nband).T
        fermi_energy = np.sort(energies, axis = None)[(int(self.filling+self.norbital)//2)*self.nsites - 1]
        fermi_energy=0
        np.savetxt('flat.txt',np.sort(flattened_energies))
        fermi_dist = self.fermi_function(flattened_energies, fermi_energy)
        fe=flattened_energies
        te=0
        for i in range(fe.shape[0]) : 
            if(fe[i]<=0) :
                te+=fe[i]
        print("gs energy=",np.sum(flattened_energies*fermi_dist)/(self.nsites),te/(self.nsites))
        exp_cost = np.zeros([8],dtype=complex)

        #for general su2 symmetric case, we might need all of tx,ty,tz. We have tz appearing in TR, so
        # we have two real couplings per bond and site
        bond_exps=np.zeros([self.nbonds,4],dtype=complex) # four elements of the matrix
        site_exps=np.zeros([self.nsites,4],dtype=complex)
        for i in range(self.nbonds):

            m_u = self.ml[i] #m_up
            n_u = self.nl[i] #n_up
            m_dd = self.norbital+self.ml[i] #m_down_dagger
            n_dd = self.norbital+self.nl[i] #n_down_dagger

            hop1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_u].T))/self.nsites
            hop2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_dd].T))/self.nsites


            pair1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_dd].T))/self.nsites


            pair2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist
                                            *self.repeated_phase[i],
                                        flattened_eigenvec[n_u].T))/self.nsites
            bond_exps[i,3]=hop1+hop2
            bond_exps[i,0]=(hop1-hop2)
            bond_exps[i,1]=(pair1-pair2)
            bond_exps[i,2]=1j*(pair1+pair2)
            print("bond=",i)
            testmat=np.array([hop1,pair1,-pair2,-hop2])
            print(testmat.reshape([2,2]))
            #print(bond_exps[i,:])

            #cost funcs later



        for m in range(self.norbital):
            m_u=m
            m_dd=self.norbital+m

            hop1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_u].T))/self.nsites
            hop2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_dd].T))/self.nsites


            pair1= (np.dot(flattened_eigenvec[m_u].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_dd].T))/self.nsites

            pair2= -(np.dot(flattened_eigenvec[m_dd].conj()
                                            *fermi_dist,
                                        flattened_eigenvec[m_u].T))/self.nsites
            site_exps[m,3]=hop1+hop2
            site_exps[m,0]=(hop1-hop2)
            site_exps[m,1]=pair1-pair2
            site_exps[m,2]=1j*(pair1+pair2)

            testmat=np.array([hop1,pair1,-pair2,-hop2])
            print("site=",m)
            print(testmat.reshape([2,2]))
            #print(site_exps[m,:])
        for n in range(n_icouplings) :
            cost=0
            c1=c2=0
            for b in range(self.nbonds):
                for delta in range(4) :
         #           print("b,delta:",b,delta)
                    c1 += transform_matrix[delta,b,n]*bond_exps[b,delta]

         #           print("exp cost:",cost)
                    #second term
                    for m in range(n_icouplings) :
                        c2 -= transform_matrix[delta,b,n]*transform_matrix[delta,b,m]*icoupling[m]
                    if(b<12 and delta==1 or delta==2) :
                        #print("n=%d,b=%d,delata=%d"%(n,b,delta))
                        temp=0
                        temp2=0
                        for m in range(n_icouplings) :
                            temp += transform_matrix[delta,b,n]*transform_matrix[delta,b,m]*icoupling[m]
                            temp2 += transform_matrix[delta,b,m]*icoupling[m]
                        #print(transform_matrix[delta,b,n]*bond_exps[b,delta].real,temp.real,c1.real,c2.real,transform_matrix[delta,b,n])
                        #print(bond_exps[b,delta].real,temp2.real)

                    #print("pcost:",-c2)


            #print("cprint:",n,c1,c2)
            cost=c1+c2
            for s in range(self.norbital):
                for delta in range(4) :
                    cost += transform_matrix[delta,self.nbonds+s,n]*site_exps[s,delta]
                    #for m in range(n_icouplings) :
                    #    cost -= transform_matrix[delta,self.nbonds+s,n]*transform_matrix[delta,self.nbonds+s,m]*icoupling[m]

            print("n=,",n,icoupling[n],cost)

    def print_corrfuncs(self, icoupling,L,fname=None):
        if(self.current_ansatz==-1) :
            raise Exception("set the ansatz, you cunt")
        #print(self.ansatzen)
        print(self.nsites)
        cfuncs=np.zeros([self.nsites,self.nband,self.nband],dtype=complex)

        transform_matrix = self.ansatzen[self. current_ansatz][1]
        n_icouplings = self.ansatzen[self. current_ansatz][0]
        couplings = np.dot(transform_matrix, icoupling)
        #couplings is [delta,bonds+sites]
        self.make_hamiltonian(couplings)

        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        
        fermi_energy = np.sort(energies, axis = None)[self.nsites*self.norbital-1]
        print("fermi energy=",fermi_energy);
        np.savetxt('flat.txt',np.sort(energies.flatten()))

        for m in range(self.nband):
            for n in range(self.nband) :
                for mom in range(self.nsites) :
                    entemp= energies[mom,:]
                    mask_array= np.where(entemp<=fermi_energy,1.0,0.0);
                    cfuncs[mom,m,n]= np.dot( np.conjugate(eigenvec[mom, m, :] *mask_array), eigenvec[mom,n,:])
        if(fname is not None) :
            np.save(fname,cfuncs);

        else :
            np.save('corrfuncs%d.npy'%(L),cfuncs);

    def write_bandstrucs(self, icoupling,fname):
        if(self.current_ansatz==-1) :
            raise Exception("set the ansatz, you cunt")
        #print(self.ansatzen)
        print(self.nsites)
        cfuncs=np.zeros([self.nsites,self.nband,self.nband],dtype=complex)
        n1= self.n_of_sites[0]
        n2= self.n_of_sites[1]
        n3= self.n_of_sites[1]

        transform_matrix = self.ansatzen[self. current_ansatz][1]
        n_icouplings = self.ansatzen[self. current_ansatz][0]
        couplings = np.dot(transform_matrix, icoupling)
        #couplings is [delta,bonds+sites]
        self.make_hamiltonian(couplings)

        sol = np.linalg.eigh(self.H)
        energies = sol[0]
        eigenvec = sol[1]
        
        fermi_energy = np.sort(energies, axis = None)[self.nsites*self.norbital-1]
        print("fermi energy=",fermi_energy);
        klist=[];
        k2list=[];
        zero= n2//2
        half=n1-1
        half=n2-1
        half=n3-1
        #gamma X
        for k in range(zero,half+1,1) :
            vec= zero*n2*n3+ k*n3+  zero
            vec2= zero*n2*n3+ (zero-(k-zero))*n3+  zero
            klist.append(vec)
            k2list.append(vec2)
        #X M
        for k in range(zero,half+1,1) :
            vec= k*n2*n3+ half*n3+  zero
            klist.append(vec)
        for k in range(zero,half+1,1) :
            vec= k*n2*n3+ half*n3+  zero
            klist.append(vec)
        #M \Gamma
        for k in range(half,zero-1,-1) :
            vec= k*n2*n3+ k*n3+  zero
            klist.append(vec)
        for k in range(zero,half+1,1) :
            vec= k*n2*n3+ k*n3+  k
            klist.append(vec)
        for k in range(half,zero-1,-1) :
            vec= k*n2*n3+ half*n3+  k
            klist.append(vec)



        bstruc=np.zeros([len(klist),self.nband])
        for k in range(len(klist)) :
                    print(k,klist[k],n1*n2*n3)
                    bstruc[k,:]= energies[klist[k],:]
                    np.set_printoptions(precision=3)
                    print(self.H[klist[k]].real)
                    print(self.H[klist[k]].imag)
                    print(energies[klist[k]])
                    print(self.H[k2list[k]].real)
                    print(self.H[k2list[k]].imag)
                    print(energies[k2list[k]])
                    sys.stdin.read(1)
        np.save(fname,bstruc);

