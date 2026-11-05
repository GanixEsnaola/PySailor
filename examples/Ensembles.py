#!/usr/bin/env python
# Ensembles.py
#
# Climatological 10 m wind from three CMIP5 ensembles (MIROC, IPSL, HadGEM)
# against ERA5 on a global grid (Saenz et al., 2020, Figs. 8 and 9)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
# Run from any directory: python examples/Ensembles.py [outputDir]
#
import os
import re
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

# All datasets share the grid: lon, lat, <name>_U, <name>_V
ens=pd.read_csv(os.path.join(data,"Ensembles.csv"))
# Area weighting of the grid points
weights=np.sqrt(np.cos(np.radians(ens["lat"])))
def dataset(name,label):
    return pd.DataFrame({"mod":label,"U":ens[name+"_U"]*weights,"V":ens[name+"_V"]*weights})
ref=dataset("ERA5","ERA5")
members=[c[:-2] for c in ens.columns if c.endswith("_U") and c!="ERA5_U"]

# Each ensemble as a single model: members stacked under the model name
# (MIROCr1 -> MIROC) and compared with ERA5 repeated once per member
modPooled=pd.concat([dataset(m,re.sub("r[0-9]+$","",m)) for m in members])
table=sailorTable(sailorIndices(ref,modPooled,pooled=True),roundDigits=3)
print(table.to_string())
table.to_csv(os.path.join(out,"Ensembles_pooled_table.csv"),index=False)
# Fig. 8 (right)
fig,ax=sailorPlot(ref,modPooled,pooled=True,sfactor=0.2,docenter=True,
                  xlim=(-1,1),ylim=(-0.5,1.5),xlab="Ux (m/s)",ylab="Uy (m/s)",
                  plotmain="CMIP5 samples vs ERA5",
                  rmseUnits=" m/s",rmseRounding=2,referenceName="ERA5")
fig.savefig(os.path.join(out,"Ensembles_1.png"),dpi=100)

# Each member as a separate model
modMembers=pd.concat([dataset(m,m) for m in members])
table=sailorTable(sailorIndices(ref,modMembers),roundDigits=1)
print(table.to_string())
table.to_csv(os.path.join(out,"Ensembles_members_table.csv"),index=False)
# Fig. 9 (right), one colour per member
colours=["red","blue","#A2CD5A","orange","#43CD80","gold","purple","pink","magenta",
         "lightblue","#8B7500","darkslateblue","coral","#008B00","#8B0A50"]
fig,ax=sailorPlot(ref,modMembers,colourList=colours,sfactor=0.01,docenter=False,
                  xlim=(-0.15,0.25),ylim=(0.19,0.61),xlab="Ux (m/s)",ylab="Uy (m/s)",
                  plotmain="CMIP5 samples vs ERA5",
                  rmseUnits=" m/s",rmseRounding=2,referenceName="ERA5")
fig.savefig(os.path.join(out,"Ensembles_2.png"),dpi=100)
