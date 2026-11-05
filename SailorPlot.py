#!/usr/bin/env python
# SailorPlot.py
#
# Sailor diagram (Saenz et al., 2020, Geosci. Model Dev. 13, 3221-3240,
# and corrigendum)
#
# Ported from the R package SailoR 1.2 (CRAN, 2020), modified in 2026.
# Copyright (C) 2020 J. Saenz, S. Carreno-Madinabeitia, S. J. Gonzalez-Roji,
#   G. Esnaola, G. Ibarra-Berastegi, A. Ulazia
# Licensed under the GNU General Public License version 3 (see LICENSE).
#
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from SailorIndices import sailorIndices

# Default colours: reference first, then one per model
refColour="#666666"
defaultColours=["red","blue","#A2CD5A","orange","#43CD80","gold","purple","pink"]
# Line types with the names used in the R package, others passed to matplotlib
rLineStyles={"dotdash":"dashdot","longdash":(0,(8,3)),"twodash":(0,(6,2,2,2))}

# Points of the ellipse with semi-axes sigma along the EOFs (columns of E),
# scaled by sfactor and centred at centre. If axisOnly, only the major axis.
def ellipsePoints(sigma,E,centre,sfactor,axisOnly=False):
    if axisOnly:
        Pe=np.array([[sigma[0],0.],[-sigma[0],0.]])
    else:
        theta=np.linspace(0.,2.*np.pi,101)
        Pe=np.column_stack([sigma[0]*np.cos(theta),sigma[1]*np.sin(theta)])
    # From the space of the PCs back to physical space
    return np.matmul(Pe*sfactor,np.transpose(E))+centre[np.newaxis,:]

# Sailor diagram of the models in mod against the reference ref
# (DataFrames with columns mod, U, V; see sailorIndices).
# For each model: its ellipse and major axis (colour), the ellipse and major
# axis of the reference (dashed grey) and its mean (bias). Ellipses are
# centred at the model mean, or at the reference mean if docenter.
# refOnce: draw the reference ellipse and axis once, at the reference mean,
# instead of once per model (useful when the model means nearly coincide)
# colourList: list of colours or a matplotlib colormap name
# lineStyles, biasMarkers: one per model (recycled if shorter)
# style: a matplotlib style name, e.g. "seaborn-whitegrid"
# ax: existing axes to draw on; otherwise a new figure is created
# Returns (fig, ax).
def sailorPlot(ref,mod,pooled=False,colourList=None,sfactor=1.,docenter=False,
               xlim=None,ylim=None,xlab=None,ylab=None,plotmain=None,
               plotLegend=True,plotRMSELegend=True,plotScaleLegend=True,
               rmseUnits="",rmseRounding=2,referenceName=None,
               lineStyles=None,biasMarkers=None,style=None,ax=None,figsize=(8,8),
               refOnce=False):
    UV=sailorIndices(ref,mod,pooled)
    names=list(UV)
    nModels=len(names)
    # Colours, line styles and bias markers per model
    if colourList is None:
        colours=list(defaultColours)
    elif isinstance(colourList,str):
        colours=list(matplotlib.colormaps[colourList](np.linspace(0.,1.,nModels)))
    else:
        colours=list(colourList)
    if len(colours)<nModels:
        # Extend with distinct colours instead of repeating
        extra=matplotlib.colormaps["tab20"](np.linspace(0.,1.,20))
        colours=colours+list(extra[:nModels-len(colours)])
    if lineStyles is None:
        lineStyles=["solid"]
    # Dashed is reserved for the reference
    lineStyles=["solid" if ls in ("dashed","--") else rLineStyles.get(ls,ls) for ls in lineStyles]
    lineStyles=[lineStyles[i%len(lineStyles)] for i in range(nModels)]
    if biasMarkers is None:
        biasMarkers=["o"]
    biasMarkers=[biasMarkers[i%len(biasMarkers)] for i in range(nModels)]
    # Bias markers are open when centred (they sit on the ellipses), filled otherwise
    markerSize=14 if docenter else 10
    # Legends get one column per 8 entries
    legendColumns=int(np.ceil((nModels+1)/8.))

    with plt.style.context(style if style is not None else "default"):
        if ax is None:
            fig,ax=plt.subplots(figsize=figsize)
        else:
            fig=ax.figure
        # Equal scales (as R's asp=1). If both limits are given, the shorter
        # range is widened to a square around its centre
        if xlim is not None and ylim is not None:
            half=max(xlim[1]-xlim[0],ylim[1]-ylim[0])/2.
            xc=(xlim[0]+xlim[1])/2.
            yc=(ylim[0]+ylim[1])/2.
            xlim=(xc-half,xc+half)
            ylim=(yc-half,yc+half)
            ax.set_aspect("equal",adjustable="box")
        else:
            ax.set_aspect("equal",adjustable="datalim")
        ax.grid(True,linestyle=":",alpha=0.6)
        meanU=next(iter(UV.values())).getMeanU()
        # Ellipses and major axes
        for i,name in enumerate(names):
            uv=UV[name]
            centre=meanU if docenter else uv.getMeanV()
            for axisOnly in (False,True):
                P=ellipsePoints(uv.getSigmaV(),uv.getEv(),centre,sfactor,axisOnly)
                ax.plot(P[:,0],P[:,1],color=colours[i],linestyle=lineStyles[i],linewidth=2)
                if not refOnce:
                    P=ellipsePoints(uv.getSigmaU(),uv.getEu(),centre,sfactor,axisOnly)
                    ax.plot(P[:,0],P[:,1],color=refColour,linestyle="dashed",linewidth=2)
        if refOnce:
            uv=UV[names[0]]
            for axisOnly in (False,True):
                P=ellipsePoints(uv.getSigmaU(),uv.getEu(),meanU,sfactor,axisOnly)
                ax.plot(P[:,0],P[:,1],color=refColour,linestyle="dashed",linewidth=2)
        # Model means (bias) on top of the ellipses, then the reference mean
        for i,name in enumerate(names):
            meanV=UV[name].getMeanV()
            ax.plot(meanV[0],meanV[1],linestyle="none",marker=biasMarkers[i],markersize=markerSize,
                    markeredgewidth=2.5,color=colours[i],
                    markerfacecolor="none" if docenter else colours[i])
        ax.plot(meanU[0],meanU[1],linestyle="none",marker="s",markersize=10,color=refColour)
        if xlim is not None:
            ax.set_xlim(xlim)
        if ylim is not None:
            ax.set_ylim(ylim)
        ax.set_xlabel(xlab if xlab is not None else "")
        ax.set_ylabel(ylab if ylab is not None else "")
        if plotmain is not None:
            ax.set_title(plotmain,fontweight="bold")
        # Legend of the datasets
        if plotLegend:
            refName=referenceName if referenceName is not None else str(ref["mod"].iloc[0])
            handles=[plt.Line2D([],[],color=refColour,linestyle="dashed",linewidth=2)]
            handles+=[plt.Line2D([],[],color=colours[i],linestyle=lineStyles[i],linewidth=2)
                      for i in range(nModels)]
            legend=ax.legend(handles,[refName]+[str(n) for n in names],loc="upper right",
                             ncol=legendColumns)
            ax.add_artist(legend)
        # Legend of the RMSE, with the same markers as in the plot
        if plotRMSELegend:
            handles=[plt.Line2D([],[],color=colours[i],linestyle=lineStyles[i],linewidth=2,
                                marker=biasMarkers[i],markeredgewidth=2,
                                markerfacecolor="none" if docenter else colours[i])
                     for i in range(nModels)]
            labels=[str(n)+": "+"{:.{}f}".format(UV[n].getRMSE(),rmseRounding)+rmseUnits for n in names]
            legend=ax.legend(handles,labels,loc="lower right",title="RMSE",ncol=legendColumns)
            legend.get_title().set_fontweight("bold")
            ax.add_artist(legend)
        # Scale factor
        if plotScaleLegend:
            ax.text(0.02,0.98,"Scaled with "+"{:g}".format(sfactor),transform=ax.transAxes,
                    verticalalignment="top",fontweight="bold")
    return fig,ax
