#!/usr/bin/env python
# Sailor.py
#
# Diagnostics of the Sailor diagram for two-dimensional vector data
# (Saenz et al., 2020, Geosci. Model Dev. 13, 3221-3240, and corrigendum)
#
# Copyright (c) 2024-2026 Jon Saenz and Ganix Esnaola
# Licensed under the MIT License (see LICENSE-MIT).
#
import numpy as np
import math 
import numpy.linalg as linalg
import scipy.linalg as spla


# SVD based PC computations of centered Uc (in columns)
class SVDPCs:
    def __init__(self,Uanom):
        self.dN=Uanom.shape[0]*1.
        # SVD decomposition of anomalies
        Au,Lhu,Eu=linalg.svd(Uanom,full_matrices=False)
        # Extract information from SVD and build the PC structures
        # These are now the eigenvalues (variances) from singular values
        self.lambdasU=Lhu*Lhu/self.dN
        # fraction of variance
        self.fracVarU=self.lambdasU/self.lambdasU.sum()
        # These are the standard deviations of PCs
        self.sdUs=np.sqrt(self.lambdasU)
        # These are the rotation matrices (originally transposed from SVD)
        self.Eu=np.transpose(Eu)
        # SVD does not fix the sign of the EOFs. Make EOF1 point east
        # (non-negative zonal component; north if exactly meridional)
        # and EOF2 follow it counterclockwise (det=+1, as Eq. (9) assumes)
        if self.Eu[0,0]<0. or (self.Eu[0,0]==0. and self.Eu[1,0]<0.):
            self.Eu[:,0]=-self.Eu[:,0]
            Au[:,0]=-Au[:,0]
        if linalg.det(self.Eu)<0.:
            self.Eu[:,1]=-self.Eu[:,1]
            Au[:,1]=-Au[:,1]
        self.Eu1=self.Eu[:,0]
        self.Eu2=self.Eu[:,1]
        # And these are (in columns) the PCs
        # First, scale them as variances
        self.PCu1=Au[:,0]*Lhu[0]
        self.PCu2=Au[:,1]*Lhu[1]
        # Total Variance U
        self.TotVarU=self.PCu1.var()+self.PCu2.var()
    def getLambdas(self):
        # get eigenvalues
        return self.lambdasU
    def getVarFrac(self):
        # get fractions of variance
        return self.fracVarU
    def getSDPCs(self):
        # get standard deviations of PCs
        return self.sdUs
    def getE(self):
        # get rotation matrix (EOFs by columns)
        return self.Eu
    def getEOFs(self):
        # get EOFs as individual vectors
        return (self.Eu1,self.Eu2)
    def getPCs(self):
        # get both principal components as arrays
        return (self.PCu1,self.PCu2)
    def getTotalVariance(self):
        # get total variance of the field
        return self.TotVarU
        
# U (reference) and V (model) must be paired: row t of both
# corresponds to the same time. Time alignment is the caller's job.
class UVError:
    def __init__(self,U,V):
        if U.ndim!=2 or U.shape[1]!=2 or U.shape!=V.shape:
            raise ValueError("U and V must be paired Nx2 arrays of the same shape, got "
                             +str(U.shape)+" and "+str(V.shape))
        # Get the number of rows (sample size) as double precission
        self.dN=float(U.shape[0])
        # Means (columnwise)
        self.meanU=U.mean(0)
        self.meanV=V.mean(0)
        # Get anomalies (remove mean by columns)
        self.Uanom=U-self.meanU[np.newaxis,:]
        self.Vanom=V-self.meanV[np.newaxis,:]
        # standard deviations of U/V in zonal/mean X/Y directions
        self.Ustd=self.Uanom.std(0)
        self.Vstd=self.Vanom.std(0)
        # EOFs of U
        # SVD decomposition of anomalies
        self.svdu=SVDPCs(self.Uanom)
        self.lambdasU=self.svdu.getLambdas()
        # fraction of variance
        self.fracVarU=self.svdu.getVarFrac()
        # These are the standard deviations of PCs
        self.sdUs=self.svdu.sdUs
        # These are the rotation matrices
        self.Eu=self.svdu.getE()
        self.Eu1, self.Eu2=self.svdu.getEOFs()
        # And these are the PCs
        self.PCu1,self.PCu2=self.svdu.getPCs()
        # Total Variance U
        self.TotVarU=self.svdu.getTotalVariance()
        # EOFs for V too
        # SVD decomposition of V anomalies (same as above)
        self.svdv=SVDPCs(self.Vanom)
        self.lambdasV=self.svdv.getLambdas()
        # fraction of variance
        self.fracVarV=self.svdv.getVarFrac()
        # These are the standard deviations of PCs
        self.sdVs=self.svdv.sdUs
        # These are the rotation matrices
        self.Ev=self.svdv.getE()
        self.Ev1, self.Ev2=self.svdv.getEOFs()
        # And these are the PCs
        self.PCv1,self.PCv2=self.svdv.getPCs()
        # Total Variance V
        self.TotVarV=self.svdv.getTotalVariance()
        # Eccentricity can be computed 
        self.epsilonU=math.sqrt(1-self.sdUs[1]*self.sdUs[1]/ \
                                (self.sdUs[0]*self.sdUs[0]))
        self.epsilonV=math.sqrt(1-self.sdVs[1]*self.sdVs[1]/ \
                                (self.sdVs[0]*self.sdVs[0]))
        self.congruence1=abs((self.Eu1*self.Ev1).sum())
        ######BUG CORRECTION#####
        # Rotations. EOF1 has a non-negative zonal component, so these
        # are in [-pi/2,pi/2]
        self.thetau=np.arctan2(self.Eu1[1],self.Eu1[0])
        self.thetav=np.arctan2(self.Ev1[1],self.Ev1[0])
        self.Rvu=np.matmul(self.Ev,np.transpose(self.Eu))
        # Relative rotation of the axes, reduced to (-pi/2,pi/2]
        self.thetavu=np.arctan2(self.Rvu[1,0],self.Rvu[0,0])
        if self.thetavu>math.pi/2.:
            self.thetavu=self.thetavu-math.pi
        elif self.thetavu<=-math.pi/2.:
            self.thetavu=self.thetavu+math.pi
        ######BUG CORRECTION#####
        # RMSE
        self.DeltaUV=np.matmul(np.transpose(U-V),(U-V))/self.dN
        self.error2=linalg.norm(self.DeltaUV,"fro")
        self.RMSE=math.sqrt(self.error2)
        ##CCA-based R2
        Sxx=np.matmul(np.transpose(self.Uanom),self.Uanom)/self.dN
        Syy=np.matmul(np.transpose(self.Vanom),self.Vanom)/self.dN
        Sxy=np.matmul(np.transpose(self.Uanom),self.Vanom)/self.dN
        Syx=np.transpose(Sxy)
        lhs=np.matmul(np.matmul(spla.fractional_matrix_power(Sxx,-0.5),\
                                Sxy),\
                  spla.fractional_matrix_power(Syy,-0.5))
        ucca,dcca,vccat=linalg.svd(lhs,full_matrices=False)
        # Canonical correlation coefficients
        self.ccacoefs=dcca
        # Based on them, R2 in 2D
        self.R2=(dcca*dcca).sum()
    def getMeanU(self):
        # Mean of dataset U
        return self.meanU
    def getMeanV(self):
        # mean of dataset V
        return self.meanV
    def getTotalVarianceU(self):
        # variance of U
        return self.TotVarU
    def getTotalVarianceV(self):
        # Variance of V
        return self.TotVarV
    def getEu(self):
        # EOF matrix of U (columns are EOFs)
        return self.Eu
    def getEv(self):
        # EOF matrix of V (columns are EOFs)
        return self.Ev
    def getRvu(self):
        # Rotation matrix
        return self.Rvu
    def getSigmaU(self):
        # standard deviations of principal axes 1 and 2 (U)
        return self.sdUs
    def getSigmaV(self):
        # standard deviations of PCs 1 and 2 (V)
        return self.sdVs
    def getSDU(self):
        # standard deviations of U dataset in zonal/meridional directions
        return self.Ustd
    def getSDV(self):
        # standard deviations of V dataset in zonal/meridional directions
        return self.Vstd
    def getThetaU(self):
        # Get the angle of E1 (U) with zonal axis
        return self.thetau
    def getThetaV(self):
        # Get the angle of E1 (V) with zonal axis
        return self.thetav
    def getThetaVU(self):
        # Get relative rotation of E1v as measured from E1u
        return self.thetavu
    def getEccentricityU(self):
        # Eccentricity of U dataset
        return self.epsilonU
    def getEccentricityV(self):
        # Eccentricity of V dataset
        return self.epsilonV
    def getRMSE(self):
        # RMSE between U and V
        return self.RMSE
    def getCCAcoefs(self):
        # CCA coefficients
        return self.ccacoefs
    def getR2vec(self):
        # 2D correlation coefficient between U and V (CCA-based)
        return self.R2
    def getCongruence(self):
        # Congruence coefficient between E1u and E1v
        return self.congruence1
    def getAll(self):
        theValues={
            "meanU":self.meanU,"meanV":self.meanV,
            "biasMag":np.sqrt(((self.meanU-self.meanV)*(self.meanU-self.meanV)).sum()),
            "TotVarU":self.TotVarU,"TotVarV":self.TotVarV,
            "Eu":self.Eu,"Ev":self.Ev,"Rvu":self.Rvu,
            "Sigmau":self.sdUs,"Sigmav":self.sdVs,
            "sdUx":self.Ustd[0],"sdUy":self.Ustd[1],
            "sdVx":self.Vstd[0],"sdVy":self.Vstd[1],
            "thetau":self.thetau,"thetav":self.thetav,"thetavu":self.thetavu,
            "EccentricityU":self.epsilonU,"EccentricityV":self.epsilonV,
            "congruenceEOF":self.congruence1,
            "RMSE":self.RMSE,
            "R2vec":self.R2,
            "CCAcoefs":self.ccacoefs
            }
        return theValues
         

