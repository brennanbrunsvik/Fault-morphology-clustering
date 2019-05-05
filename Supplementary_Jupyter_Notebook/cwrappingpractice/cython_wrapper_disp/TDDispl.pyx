import numpy as np
cimport numpy as cnp

cnp.import_array()

cdef extern from "examples.h":
    void DispInHalfSpace(double Displ[3], double X,double Y,double Z,double P1[3],double P2[3],double P3[3],double Ss,double Ds,double Ts,double nu)

def py_DispInHalfSpace(cnp.ndarray[cnp.float64_t,ndim=1] Displ,double X,double Y,double Z,cnp.ndarray[cnp.float64_t,ndim=1] P1,cnp.ndarray[cnp.float64_t,ndim=1] P2,cnp.ndarray[cnp.float64_t,ndim=1] P3,double Ss,double Ds,double Ts,double nu) -> None:
    DispInHalfSpace(&Displ[0], X, Y, Z, &P1[0] ,&P2[0], &P3[0],Ss, Ds, Ts, nu)   

def py_DispInHalfSpaceLoop(cnp.ndarray[cnp.float64_t,ndim=2] Displ,cnp.ndarray[cnp.float64_t,ndim=1] X,cnp.ndarray[cnp.float64_t,ndim=1] Y,cnp.ndarray[cnp.float64_t,ndim=1] Z,cnp.ndarray[cnp.float64_t,ndim=1] P1,cnp.ndarray[cnp.float64_t,ndim=1] P2,cnp.ndarray[cnp.float64_t,ndim=1] P3,double Ss,double Ds,double Ts,double nu) -> None:
    # it appears that the c function itself takes so much time that it is not pointful to make a cython optimized loop.
    # a python loop is quick enough to not noticeably slow the code. 
    cdef int i
    cdef cnp.ndarray[cnp.float64_t,ndim=1] tempDisp
    
    tempDisp = np.array([0, 0, 0], dtype = np.float64)
    for i in range(X.size):
        DispInHalfSpace(&tempDisp[0], X[i], Y[i], Z[i], &P1[0] ,&P2[0], &P3[0],Ss, Ds, Ts, nu)  
        #Displ[i, :] = tempDisp[:]
        Displ[i] = tempDisp



