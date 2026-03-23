
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os
import pandas as pd

def save_show(save_dir, save_name, show):
    print("save_show called")
    # Save if save_name is provided
    if save_name is not None:
        # Use current directory if save_path not provided
        if save_dir is None:
            save_dir = "."
        # Create directory if it doesn't exist
        os.makedirs(save_dir, exist_ok=True)
        # Save as PDF
        plt.savefig(os.path.join(save_dir, f"{save_name}.pdf"), bbox_inches='tight')

    if show:
        plt.show()
    plt.close()  
    return

def plot_densities(vect_list, 
                   group_label_list, 
                   title, 
                   xaxis, 
                   save_name=None, 
                   save_dir=None,
                   show=True,
                   xlim=[0, 1],
                   smooth=1,
                   alpha=0.2,
                   color_dict=None,
                   hist=False,       # NEW ARG
                   bins=30):         # optional bins for histogram
    """
    Plot density (KDE) and optionally histogram for multiple groups.
    """
    plt.figure(figsize=(8, 4))
    
    for vect, group_label in zip(vect_list, group_label_list):
        color = color_dict[group_label] if color_dict and group_label in color_dict else None
        
        if hist:
            # Plot histogram
            plt.hist(vect, bins=bins, density=True, alpha=0.5, 
                     label=group_label.replace("_", " "), color=color, edgecolor="black")
        
        # Always plot KDE (overlay)
        sns.kdeplot(vect, fill=True, bw_adjust=smooth, alpha=alpha, 
                    label=None if hist else group_label.replace("_", " "), 
                    color=color)
    
    plt.xlabel(xaxis)
    plt.ylabel("Density")
    plt.title(title)
    if xlim is not None:
        margin = (xlim[1] - xlim[0]) * 0.05  # 5% padding
        plt.xlim(xlim[0] - margin, xlim[1] + margin)
    plt.legend(loc="upper left")    
    save_show(save_name=save_name, 
              save_dir=save_dir, 
              show=show)
    return
