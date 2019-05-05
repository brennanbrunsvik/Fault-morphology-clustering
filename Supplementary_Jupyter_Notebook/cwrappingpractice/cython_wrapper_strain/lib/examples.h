#ifndef EXAMPLES_H
#define EXAMPLES_H

void StrainInHalfSpace(double Stress[6], double Strain[6], double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double Ss, double Ds, double Ts, double mu, double lambda);
void TDstress_HarFunc_inStrnHS(double StsFSC[6],double StrFSC[6], double X, double Y, double Z, double P1[3] ,double P2[3], double P3[3], double Ss, double Ds, double Ts, double mu, double lambda);
void TDstressFS_inStrnHS(double StsMS[6], double StrMS[6],  double X, double Y, double Z, double P1[3], double P2[3], double P3[3], double by, double bz, double bx, double mu, double lambda);
void       CoordTrans_inStrnHS(double newVal[3], double x_shift, double y_shift, double z_shift, double RotMat[3][3]); 
void      TriModeFind_inStrnHS(int TrimMode[1], double x,double y,double z,double p1_a,double p1_b, double p2_a,double p2_b,double p3_a,double p3_b);
void         TDSetupS_inStrnHS(double x,double y,double z,double alpha,double bx,double by,double bz,double nu, double TriVertex[3],double SideVec[3],double e[6]);
void     AngDisStrain_inStrnHS(double x, double y, double z, double alpha, double bx, double by, double bz, double nu, double e[6]);
void        TensTrans_inStrnHS(double e_in[6], double e_out[6], double B[3][3]);
void    AngSetupFSC_S_inStrnHS(double Stress1[6],double Strain1[6], double X,double Y,double Z,double bX,double bY,double bZ,double Pt1[3], double Pt2[3], double mu,double lambda);
void  AngDisStrainFSC_inStrnHS(double y1, double y2, double y3, double beta, double b1, double b2, double b3, double nu, double a, double Strain[6]);

#endif
