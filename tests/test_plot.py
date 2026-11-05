# The reference ellipse is drawn once per model, or once with refOnce
import os
import matplotlib
matplotlib.use("Agg")
import pandas as pd
from SailorPlot import sailorPlot

data=os.path.join(os.path.dirname(os.path.abspath(__file__)),"..","data")

def test_refOnce():
    ref=pd.read_csv(os.path.join(data,"Synthetic_ref.csv"))
    mod=pd.read_csv(os.path.join(data,"Synthetic_mod.csv"))
    nModels=len(pd.unique(mod["mod"]))
    for refOnce,nReference in ((False,nModels),(True,1)):
        fig,ax=sailorPlot(ref,mod,refOnce=refOnce,plotLegend=False,plotRMSELegend=False)
        # Ellipse and major axis of every model, then of the reference
        dashed=[l for l in ax.get_lines() if l.get_linestyle()=="--"]
        assert len(dashed)==2*nReference
        assert len(ax.get_lines())-len(dashed)==2*nModels+nModels+1
        matplotlib.pyplot.close(fig)
