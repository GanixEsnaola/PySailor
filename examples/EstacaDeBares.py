#!/usr/bin/env python
# EstacaDeBares.py
#
# Wave energy flux from several models against the Estaca de Bares buoy
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
# Run from any directory: python examples/EstacaDeBares.py [outputDir]
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

ref=pd.read_csv(os.path.join(data,"EstacaDeBares_ref.csv"))
mod=pd.read_csv(os.path.join(data,"EstacaDeBares_mod.csv"))

# Indices and summary table
table=sailorTable(sailorIndices(ref,mod),roundDigits=2)
print(table.to_string())
table.to_csv(os.path.join(out,"EstacaDeBares_table.csv"),index=False)

fig,ax=sailorPlot(ref,mod,sfactor=0.25,docenter=True,
                  xlim=(0,25),ylim=(0,25),
                  plotmain="Wave Energy Flux. Estaca De Bares.",
                  rmseUnits=" W/m",rmseRounding=0,referenceName="Buoy")
fig.savefig(os.path.join(out,"EstacaDeBares_1.png"),dpi=100)
