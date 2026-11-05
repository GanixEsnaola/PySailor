# Properties that UVError must satisfy for any paired input
import math
import numpy as np
import pytest
from Sailor import UVError

def rotation(angle):
    return np.array([[math.cos(angle),-math.sin(angle)],[math.sin(angle),math.cos(angle)]])

def randomPairs(seed,n=300):
    rng=np.random.default_rng(seed)
    U=rng.normal(size=(n,2))*rng.uniform(0.2,5.,2)
    V=rng.normal(size=(n,2))*rng.uniform(0.2,5.,2)
    U=np.matmul(U,np.transpose(rotation(rng.uniform(-math.pi,math.pi))))
    V=np.matmul(V,np.transpose(rotation(rng.uniform(-math.pi,math.pi))))+U
    return U,V

@pytest.mark.parametrize("seed",range(200))
def test_angles_and_eofs(seed):
    uv=UVError(*randomPairs(seed))
    # EOF bases are rotations with EOF1 pointing east
    for E in (uv.getEu(),uv.getEv()):
        assert abs(np.linalg.det(E)-1.)<1e-12
        assert E[0,0]>=0.
    # Angle ranges
    assert -math.pi/2.<=uv.getThetaU()<=math.pi/2.
    assert -math.pi/2.<=uv.getThetaV()<=math.pi/2.
    assert -math.pi/2.<uv.getThetaVU()<=math.pi/2.
    # thetavu is thetav-thetau up to pi, and its cosine is the congruence g11
    d=(uv.getThetaV()-uv.getThetaU()-uv.getThetaVU())/math.pi
    assert abs(d-round(d))<1e-9
    assert abs(math.cos(uv.getThetaVU())-uv.getCongruence())<1e-12
    # PCs and EOFs rebuild the anomalies
    P=np.column_stack(uv.svdu.getPCs())
    assert np.allclose(np.matmul(P,np.transpose(uv.getEu())),uv.Uanom)

def test_same_axes_on_opposite_sides_of_north():
    # Major axes at 85 and 95 degrees: relative rotation is 10 degrees
    rng=np.random.default_rng(0)
    B=rng.normal(size=(2000,2))*[3.,1.]
    uv=UVError(np.matmul(B,np.transpose(rotation(math.radians(85.)))),
               np.matmul(B,np.transpose(rotation(math.radians(95.)))))
    assert abs(uv.getThetaVU()-math.radians(10.))<1e-12

def test_perfect_and_linear_models():
    U,V=randomPairs(1)
    uv=UVError(U,U)
    assert uv.getRMSE()==0. and abs(uv.getR2vec()-2.)<1e-12 and abs(uv.getThetaVU())<1e-12
    # Any invertible linear transformation plus a shift keeps R2=2
    uv=UVError(U,np.matmul(U,[[1.5,0.3],[-0.2,0.7]])+[1.,2.])
    assert abs(uv.getR2vec()-2.)<1e-12

def test_unpaired_input_is_rejected():
    U,V=randomPairs(2)
    for W in (V[:-1],np.column_stack([V,V[:,0]])):
        with pytest.raises(ValueError):
            UVError(U,W)
