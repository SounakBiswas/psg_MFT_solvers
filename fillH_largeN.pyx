import numpy as np
import sys

def fillH(int norbital, double complex [:,:,::1]temp_ham, double beta, double [:]jlist, double [:]mulist,
        int nbonds, int[:] ml,int[:]nl, double complex [:,:] phaselist,
          int nsites):
    cdef int point
    cdef int i,  m,  n
    cdef int m1,n1,m2,n2;
    cdef double u0,ux,uy,uz;
    for point in range(nsites) :
        for m in range(norbital) :
           for n in range(norbital) :
               temp_ham[point,m,n]=0
    for point in range(nsites): #sum over momenta
        for i in range(nbonds):
            m = ml[i] #m_up
            n = nl[i]
            temp_ham[point, m, n] = temp_ham[point, m, n] + beta*jlist[i]* phaselist[i, point]
            temp_ham[point, n, m] = temp_ham[point, n, m] + beta*jlist[i]* phaselist[i, point].conjugate()


        for i in range(norbital):
            mu = mulist[i]
            temp_ham[point, i, i] = temp_ham[point, i, i]+ mu
