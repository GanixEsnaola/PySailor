#!/usr/bin/env python
# WRF.py
#
# Vertically integrated water vapour transport from two WRF runs and ERA-Interim
# against soundings (Saenz et al., 2020, Fig. 6 right, Table 2)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
# Run from any directory: python examples/WRF.py [outputDir]
#
import os
import sys
here=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(here,".."))
import pandas as pd
from SailorIndices import sailorIndices, sailorTable
from SailorPlot import sailorPlot

data=os.path.join(here,"..","data")
out=sys.argv[1] if len(sys.argv)>1 else os.path.join(here,"output")
os.makedirs(out,exist_ok=True)

ref=pd.read_csv(os.path.join(data,"WRF_ref.csv"))
mod=pd.read_csv(os.path.join(data,"WRF_mod.csv"))
# Keep only the times with a valid sounding, in the reference and in every model
isOK=ref[["U","V"]].notna().all(axis=1).to_numpy()
ref=ref[isOK]
mod=pd.concat([mod[mod["mod"]==name][isOK] for name in pd.unique(mod["mod"])])

# Indices and summary table (Table 2)
table=sailorTable(sailorIndices(ref,mod),roundDigits=2)
print(table.to_string())
table.to_csv(os.path.join(out,"WRF_table.csv"),index=False)

# Fig. 6 (right)
fig,ax=sailorPlot(ref,mod,sfactor=0.1,docenter=True,
                  xlim=(60,110),ylim=(0,50),xlab="Qx (kg/m/s)",ylab="Qy (kg/m/s)",
                  plotmain="Water vapour transport",
                  rmseUnits=" kg/m/s",rmseRounding=0,referenceName="Sounding")
fig.savefig(os.path.join(out,"WRF_1.png"),dpi=100)

# The same with a different line type per model
fig,ax=sailorPlot(ref,mod,sfactor=0.1,docenter=True,
                  xlim=(60,110),ylim=(0,50),xlab="Qx (kg/m/s)",ylab="Qy (kg/m/s)",
                  plotmain="Water vapour transport",plotScaleLegend=False,
                  rmseUnits=" kg/m/s",rmseRounding=0,referenceName="Sounding",
                  lineStyles=["solid","dotted","twodash"])
fig.savefig(os.path.join(out,"WRF_2.png"),dpi=100)
