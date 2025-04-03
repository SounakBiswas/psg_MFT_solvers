import numpy as np
import sys
cdef extern from "complex.h" :
    double complex cexp(double complex)


def expand(double complex [:,:,:,:,:]sfacKab, int L, int norb,  double [:,:]psgc, int bz_fac ) :
    cdef int m,n;
    cdef Py_ssize_t kx,ky,kz;
    cdef Py_ssize_t idx,idy,idz;
    cdef Py_ssize_t fKx,fKy,fKz;
    cdef int fL;
    cdef double pi = np.pi
    cdef double phase;
    
    fL=bz_fac*L
    sfac_full=np.zeros([fL,fL,fL],dtype=complex)
    cdef complex [:,:,:] sfacO= sfac_full

    for idx in range(fL) :
        for idy in range(fL) :
            for idz in range(fL) :
                if(fL%2==0) :
                    kx=(idx-fL/2)
                    ky=(idy-fL/2)
                    kz=(idz-fL/2)
                else :
                    kx=(idx-(fL-1)/2)
                    ky=(idy-(fL-1)/2)
                    kz=(idz-(fL-1)/2)
                fKx=(fL+kx)%L
                fKy=(fL+ky)%L
                fKz=(fL+kz)%L
                for m in range(norb) :
                    for n in range(norb) :
                        phase = (2*pi/L)*(psgc[m,0]-psgc[n,0])*kx;
                        phase += (2*pi/L)*(psgc[m,1]-psgc[n,1])*ky;
                        phase += (2*pi/L)*(psgc[m,2]-psgc[n,2])*kz;
                        sfacO[idx,idy,idz]=sfacO[idx,idy,idz]+sfacKab[fKx,fKy,fKz,m,n]*cexp(1j*phase)
    return sfac_full
    
