import numpy as np
#from psg_hamiltonian_class import psg_lattice
from class_backup import psg_lattice
import scipy.optimize as opt
import scipy
import sys
import bikagome_maker
np.set_printoptions(precision=16)
def get_mf_params(psg,x) :
    b,s,_= psg.calc_exps(x,u1=True)
    b1=(np.mean(b[0:3,3])+np.mean(b[6:9,3]))/2
    b2=(np.mean(b[3:6,3])+np.mean(b[9:12,3]))/2
    b3=np.mean(b[12:18,3])
    s=np.mean(s[:,3])
    return b1,b2,b3,s

def cost_func_u1(psg,x) : 
    b1,b2,b3,_= get_mf_params(psg,x)
    return [b1-x[0],b2-x[1],b3-x[2]]



    
#        
def cost_func_u1_energy(psg,x) : 
    _,_,e= psg.calc_exps(x,u1=True)
    return e
        
def callback_func(x,f) :
    print("callback params:,",x.real)
    print("callback costs,",cost_func_u1_root_full(bikag,x))
    x=0
    #y=x
    #y[1]=x[0]
    #y[0]=x[1]
    #print("reversed")
    #print(cost_func_u1_root(bikag,y));
    #sys.stdin.read(1)
def callback_func_min(x) :
    print("callback params:,",x.real)
    print("callback costs,",cost_func_u1_root_full(bikag,x.real))
    #print("energies:",12*x[0]*(2*f[0]+x[0]),12*x[0]**2,trillium_psg.calc_energy(x))
    #sysen=trillium.calc_energy(x)
    #print("energies:",en,sysen/2,sysen+en)
    return 0

L=20

ji=1.0
jo=1.5
jint=1.1
bikag = bikagome_maker.make_psg_lattice(L,ji,jo,jint);
nbonds=bikag.nbonds
norb=bikag.norbital
print (nbonds,norb)

#psg IV
#notation : ux, uy, mu_x, mu_y
n_icouplings=3
tmat= np.zeros([4,nbonds+norb,n_icouplings])
for delta in [0,1,2,3] :
    tmat[delta,:,:]=0

tmat[3,0:3,0]=1
tmat[3,6:9,0]=1
tmat[3,3:6,1]=1
tmat[3,9:12,1]=1
tmat[3,12:18,2]=1



bikag.add_ansatz(n_icouplings,tmat)
bikag.initH()
bikag.current_ansatz=0
mine=0
guess=[0.66,0,0]

guess =[ 0.1673601237692682, 0.49256576294081855, 0.29845018703527015]
#guess =[ 0.00292100, 1.0000185, 0.001120]
b,s,e=bikag.calc_exps(guess,u1=True)
b1,b2,b3,s=get_mf_params(bikag,guess)
cost=cost_func_u1(bikag,guess)
print("cost:",cost)
print("params:",b1,b2,b3,s)
print("energy:",cost_func_u1_energy(bikag,guess)/(4*6))





