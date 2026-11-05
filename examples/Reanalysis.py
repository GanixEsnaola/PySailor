#!/usr/bin/env python
# Reanalysis.py
#
# Climatological 10 m wind from four reanalyses against ERA5 on a global grid
# (Saenz et al., 2020, Fig. 7)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
# Run from any directory: python examples/Reanalysis.py [outputDir]
#
import os
import sys
here=os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0,os.path.join(here,".."))
import numpy as np
import pandas as pd
from SailorIndices import sailorIndices, sailorTable
from SailorPlot import sailorPlot

data=os.path.join(here,"..","data")
out=sys.argv[1] if len(sys.argv)>1 else os.path.join(here,"output")
os.makedirs(out,exist_ok=True)

ref=pd.read_csv(os.path.join(data,"Reanalysis_ref.csv"))
mod=pd.read_csv(os.path.join(data,"Reanalysis_mod.csv"))
# Area weighting of the grid points: both reference and models are
# multiplied by sqrt(cos(latitude))
for df in (ref,mod):
    weights=np.sqrt(np.cos(np.radians(df["lat"])))
    df["U"]=df["U"]*weights
    df["V"]=df["V"]*weights

# Indices and summary table
table=sailorTable(sailorIndices(ref,mod),roundDigits=3)
print(table.to_string())
table.to_csv(os.path.join(out,"Reanalysis_table.csv"),index=False)

# Fig. 7: uncentred (left) and centred at the reference mean (right)
for i,docenter in ((1,False),(2,True)):
    fig,ax=sailorPlot(ref,mod,sfactor=0.15,docenter=docenter,
                      xlim=(-0.5,1),ylim=(-1,0.5),xlab="Ux (m/s)",ylab="Uy (m/s)",
                      plotmain="Reanalyses",
                      rmseUnits=" m/s",rmseRounding=1,referenceName="ERA5")
    fig.savefig(os.path.join(out,"Reanalysis_"+str(i)+".png"),dpi=100)
