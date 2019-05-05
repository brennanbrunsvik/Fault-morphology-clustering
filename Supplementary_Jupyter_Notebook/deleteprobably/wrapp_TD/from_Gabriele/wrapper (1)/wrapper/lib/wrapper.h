#ifndef WRAPPER_H
#define WRAPPER_H

void        TDdispFS_inDispHS(double Displ[3],    double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double nu);   
void  TDdisp_HarFunc_inDispHS(double Displ[3],    double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double nu);   
void      CoordTrans_inDispHS(double newVal[3],  double x, double y, double z, double RotMat[3][3]); 
void   TriModeFind_inDispHS(int    TriMode[1],  double x, double y, double z, double p1a, double p1b, double p2a, double p2b, double p3a, double p3b);
void        TDSetupD_inDispHS(double DispVect[3], double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double TriVertex[3], double SideVec[3]);
void     AngSetupFSC_inDispHS(double DispVect[3], double x, double y, double z, double bx, double by, double bz, double PA[3], double PB[3], double nu);
void      AngDisDisp_inDispHS(double DispVect[3], double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu); 
void   AngDisDispFSC_inDispHS(double DispVect[3], double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu, double a); 
void DispInHalfSpace(double Displ[3], double X,double Y,double Z,double P1[3],double P2[3],double P3[3],double Ss,double Ds,double Ts,double nu);


#endif
