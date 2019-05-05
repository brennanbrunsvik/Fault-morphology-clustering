/*
cdef extern from "examples.h":
    void hello(const char *name)

def py_hello(name: bytes) -> None:
    hello(name)
*/

cdef extern from "DispInHalfSpace.c":
    void DispInHalfSpace(double Displ[3], double X,double Y,double Z,double P1[3],double P2[3],double P3[3],double Ss,double Ds,double Ts,double nu)

def pyDisp(*inputs):
    DispInHalfSpace(*inputs)
