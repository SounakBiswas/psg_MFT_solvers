import numpy as np
import sys

def fillH(double complex [:,:,:]temp_ham, double [:]u0list, double [:]uxlist, double [:]uylist,double [:]uzlist, double [:]muxlist, double [:]muylist, double [:]muzlist,
        int nbonds, int[:] ml,int[:]nl, double complex [:,:] phaselist,
        int nsites, int nband):
    cdef int norbital = nband // 2
    cdef int point
    cdef int i,  m,  n
    cdef int m1,n1,m2,n2;
    cdef double u0,ux,uy,uz;
    cdef double mux, muy, muz;
    for point in range(nsites) :
        for m in range(nband) :
           for n in range(nband) :
               temp_ham[point,m,n]=0
    for point in range(nsites): #sum over momenta
        #x1=np.zeros([4,4],dtype=complex)
        #x2=np.zeros([4,4],dtype=complex)
        for i in range(nbonds):
            m1 = ml[i] #m_up
            n1 = nl[i]
            m2 = norbital+ml[i] #m_down_dagger
            n2 = norbital+nl[i]
            u0 = u0list[i]
            ux = uxlist[i]
            uy = uylist[i]
            uz = uzlist[i]
            temp_ham[point, m1, n2] = temp_ham[point, m1, n2] - (ux-uy*1j) * phaselist[i, point]
            temp_ham[point, m2, n1] = temp_ham[point, m2, n1] - (ux+uy*1j) * phaselist[i, point]

            temp_ham[point, n2, m1] = temp_ham[point, n2, m1] - (ux+uy*1j) * phaselist[i, point].conjugate()
            temp_ham[point, n1, m2] = temp_ham[point, n1, m2] - (ux-uy*1j) * phaselist[i, point].conjugate()

            temp_ham[point, m1, n1] = temp_ham[point, m1, n1] -  (1j*u0+uz) * phaselist[i,point]
            temp_ham[point, m2, n2] = temp_ham[point, m2, n2] - (1j*u0-uz)  *phaselist[i,point]

            temp_ham[point, n1, m1] = temp_ham[point, n1, m1] -  (-1j*u0+uz) * phaselist[i,point].conjugate()
            temp_ham[point, n2, m2] = temp_ham[point, n2, m2] - (-1j*u0-uz)  *phaselist[i,point].conjugate()

            #x1[ m1, n1] = x1[ m1, n1] -  (1j*u0+uz) * phaselist[i,point]
            #x1[ n1, m1] = x1[ n1, m1] -  (-1j*u0+uz) * phaselist[i,point].conjugate()

            #x2[ m1, n1] = x2[ m1, n1] - (1j*u0-uz)  *phaselist[i,point]
            #x2[ n1, m1] = x2[ n1, m1] - (-1j*u0-uz)  *phaselist[i,point].conjugate()
        #print(point)
        #np.set_printoptions(precision=3,suppress=True)
        #print("real")
        #print(np.real(x1))
        #print(np.real(x2))
        #print("imag")
        #print(np.imag(x1))
        #print(np.imag(x2))
        #sys.stdin.read(1)
        for i in range(norbital):
            mux = muxlist[i]
            muy = muylist[i]
            muz = muzlist[i]
            m1 = i #m_up
            m2 = norbital+i #m_ddd
            temp_ham[point, m1, m2] = temp_ham[point, m1, m2]+ (mux-1j*muy)
            temp_ham[point, m2, m1] = temp_ham[point, m2, m1]+ (mux+1j*muy)
            temp_ham[point, m1, m1] = temp_ham[point, m1, m1]+ (muz)
            temp_ham[point, m2, m2] = temp_ham[point, m2, m2]+ (-muz)
