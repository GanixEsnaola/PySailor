#!/usr/bin/env python
# SailorIndices.py
#
# Sailor diagram indices for several models against a reference, and
# summary table (Saenz et al., 2020, Geosci. Model Dev. 13, 3221-3240,
# and corrigendum)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
import numpy as np
import pandas as pd
from Sailor import UVError

# ref and mod are DataFrames with columns mod (name), U (zonal) and
# V (meridional). Each model in mod must be paired row by row with ref.
# If pooled is True, a model k times longer than ref (e.g. k ensemble
# members stacked under the same name) is compared with ref repeated
# k times, giving one set of indices for the whole ensemble.
# Returns a dict {model name: UVError}, in order of appearance in mod.
def sailorIndices(ref,mod,pooled=False):
    U=ref[["U","V"]].to_numpy(dtype=float)
    nU=U.shape[0]
    UV={}
    for name in pd.unique(mod["mod"]):
        V=mod.loc[mod["mod"]==name,["U","V"]].to_numpy(dtype=float)
        nV=V.shape[0]
        if nV==nU:
            UV[name]=UVError(U,V)
        elif pooled and nV%nU==0:
            UV[name]=UVError(np.tile(U,(nV//nU,1)),V)
        else:
            raise ValueError("Model "+str(name)+" has "+str(nV)+" rows and the reference "
                             +str(nU)+(", not an integer multiple" if pooled else
                                       "; use pooled=True for stacked ensemble members"))
    return UV

# Summary table of the indices from sailorIndices, one row per model plus
# a first row for the reference. Values rounded to roundDigits decimals.
def sailorTable(UV,roundDigits=2):
    first=next(iter(UV.values()))
    # Reference row: values of U, which are the same for every model
    rows=[{"modelName":"Ref",
           "sdUx":first.getSDU()[0],"sdUy":first.getSDU()[1],
           "sdVx":np.nan,"sdVy":np.nan,
           "Sigmax":first.getSigmaU()[0],"Sigmay":first.getSigmaU()[1],
           "thetau":first.getThetaU(),"thetav":np.nan,"thetavu":np.nan,
           "R2vec":2.,"biasMag":0.,"RMSE":0.,
           "Eccentricity":first.getEccentricityU(),"congruenceEOF1":1.}]
    # One row per model
    for name,uv in UV.items():
        rows.append({"modelName":name,
                     "sdUx":np.nan,"sdUy":np.nan,
                     "sdVx":uv.getSDV()[0],"sdVy":uv.getSDV()[1],
                     "Sigmax":uv.getSigmaV()[0],"Sigmay":uv.getSigmaV()[1],
                     "thetau":np.nan,"thetav":uv.getThetaV(),"thetavu":uv.getThetaVU(),
                     "R2vec":uv.getR2vec(),
                     "biasMag":np.linalg.norm(uv.getMeanU()-uv.getMeanV()),
                     "RMSE":uv.getRMSE(),
                     "Eccentricity":uv.getEccentricityV(),
                     "congruenceEOF1":uv.getCongruence()})
    table=pd.DataFrame(rows)
    # Round (adding 0. turns -0. into 0.)
    numeric=table.columns[1:]
    table[numeric]=table[numeric].round(roundDigits)+0.
    return table
