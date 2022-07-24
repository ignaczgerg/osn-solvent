from rdkit import Chem

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import seaborn as sns
from matplotlib.lines import Line2D
from sklearn.metrics import r2_score
import matplotlib.colors as colors
import matplotlib.cm as cmx

from matplotlib.patches import Circle
from matplotlib.offsetbox import (TextArea, DrawingArea, OffsetImage,
                                  AnnotationBbox)
from matplotlib.cbook import get_sample_data
import matplotlib.patches as mpatches

def smiles_concat(first_smiles: list, second_smiles: list):
    """
    Given two list of smiles, returns the their list with the concatenated version. 
    Examples: A='CO' and B='CN' returns CO.CN
    """
    if bool(first_smiles) is False:
        raise ValueError("You should know not to pass an empty list, dumbass.")
    if bool(second_smiles) is False:
        raise ValueError("You should know not to pass an empty list, dumbass.")

    full_smiles = []
    for (smile_1, smile_2) in zip(first_smiles, second_smiles):
        instance = type(Chem.rdchem.Mol())
        if (isinstance(Chem.MolFromSmiles(smile_1), instance), isinstance(Chem.MolFromSmiles(smile_2), instance)) == (True, True):
            full_smiles.append(smile_1 + "." + smile_2)
        else:
           raise ValueError(f"Invalid smiles on on {smiles_1} or {smiles_2}.")
        
    return full_smiles



def rejection_diagram(x: str, y: str, data: pd.DataFrame, x_axis: str, y_axis: str, group="solvent_name", save=None):
    """
    :x: measured data
    :y: predicted data
    :data: pandas dataframe
    :x_axis: Label name on the X-axis
    :y_axis: Label name on the Y-axis
    """
    plt.figure(figsize=(6,6), tight_layout=True)
    ax = plt.axes()
    ax.set(facecolor = "white")
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    _data = data
    x_values=_data[x].astype('float64')
    y_values=_data[y].astype('float64')
    
    # Get unique names of species
    uniq = list(set(_data[group]))

    # Set the color map to match the number of species
    z = range(1,len(uniq))
    hot = plt.get_cmap('tab20')
    cNorm  = colors.Normalize(vmin=0, vmax=len(uniq))
    scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=hot)

    markers = ["o", "v", "^", "<", ">", "1", "2",
                "3", "4", "8", "s", "p", "P", "*",
                "h", "H", "+", "x", "X", "D", "d"]

    for enum, (k, i) in enumerate(zip(uniq, markers)):
        indx = _data[group] == k
        plt.scatter(x_values[indx], y_values[indx], s=50, color=scalarMap.to_rgba(enum), label=k, marker=i)

    plt.legend(loc='lower right', bbox_to_anchor=(1, 0))
    # plt.ylim(-1,3)

    z = np.polyfit(_data[x], _data[y], 1)
    p = np.poly1d(z)
    y_hat = np.poly1d(z)(_data[x])

    r_square = r2_score(_data[y], y_hat)
    plt.plot(_data[x],p(_data[x]),"-")
    
    
    text = f"$y={z[0]:0.3f}\;x{z[1]:+0.3f}$\n$R^2 = {r2_score(_data[y],y_hat):0.3f}$"
    plt.gca().text(0.05, 0.95, text,transform=plt.gca().transAxes,
        fontsize=10, verticalalignment='center', horizontalalignment='left')

    plt.xlabel(x_axis)
    plt.ylabel(y_axis)

    if save is not None:
        plt.savefig(save)