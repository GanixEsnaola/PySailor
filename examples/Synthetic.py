#!/usr/bin/env python
# Synthetic.py
#
# Reference and four synthetic models showing the response of the Sailor
# diagram to different kinds of error (Saenz et al., 2020, Fig. 4, Table 1)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
# Run from any directory: python examples/Synthetic.py [outputDir]
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

ref=pd.read_csv(os.path.join(data,"Synthetic_ref.csv"))
mod=pd.read_csv(os.path.join(data,"Synthetic_mod.csv"))

# Indices and summary table (Table 1)
table=sailorTable(sailorIndices(ref,mod),roundDigits=4)
print(table.to_string())
table.to_csv(os.path.join(out,"Synthetic_table.csv"),index=False)

# Fig. 4: uncentred (left) and centred at the reference mean (right)
for i,docenter in ((1,False),(2,True)):
    fig,ax=sailorPlot(ref,mod,sfactor=1,docenter=docenter,
                      xlim=(0,15),ylim=(-25,10),xlab="Ux (m/s)",ylab="Uy (m/s)",
                      plotmain="Reference and synthetic models",
                      rmseUnits=" m/s",rmseRounding=2)
    fig.savefig(os.path.join(out,"Synthetic_"+str(i)+".png"),dpi=100)
