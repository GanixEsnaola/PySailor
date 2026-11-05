#!/usr/bin/env python
# Dragonera.py
#
# 10 m wind from two WRF runs, ERA-Interim and CCMP satellite winds against
# the Dragonera buoy (Saenz et al., 2020, Fig. 5)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
# Run from any directory: python examples/Dragonera.py [outputDir]
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

ref=pd.read_csv(os.path.join(data,"Dragonera_ref.csv"))
mod=pd.read_csv(os.path.join(data,"Dragonera_mod.csv"))

# Indices and summary table
table=sailorTable(sailorIndices(ref,mod),roundDigits=2)
print(table.to_string())
table.to_csv(os.path.join(out,"Dragonera_table.csv"),index=False)

# Fig. 5 (left): centred, with the area shown in the right panel
fig,ax=sailorPlot(ref,mod,sfactor=1,docenter=True,
                  xlim=(-10,10),ylim=(-10,10),xlab="Ux (m/s)",ylab="Uy (m/s)",
                  plotmain="Dragonera U10/V10",
                  rmseUnits=" m/s",rmseRounding=1,referenceName="Reference")
ax.plot([-1,1,1,-1,-1],[-1,-1,1,1,-1],color="brown",linewidth=4)
fig.savefig(os.path.join(out,"Dragonera_1.png"),dpi=100)

# Fig. 5 (right): uncentred, ellipses scaled down to show the biases
fig,ax=sailorPlot(ref,mod,sfactor=0.025,docenter=False,
                  xlim=(-1,1),ylim=(-1,1),xlab="Ux (m/s)",ylab="Uy (m/s)",
                  plotmain="Dragonera U10/V10 (scaled)",
                  rmseUnits=" m/s",rmseRounding=1,referenceName="Reference")
fig.savefig(os.path.join(out,"Dragonera_2.png"),dpi=100)
