from  matplotlib.patches import ArrowStyle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib import rcParams
from scipy.stats import circstd

def plot_CJ(fA_var, 
fA_pointSize = 20, 
fA_arrowSize = "S", 
fA_alphaC = 0.5, 
fA_alphaP = 1.0, 
fA_levels = 2, 
fA_thresh = 0.8):
    cartData_fromR_list = fA_var
    
    cartData_pd_20P = pd.DataFrame(data = cartData_fromR_list[0])
    cartData_pd_20P = cartData_pd_20P.transform(lambda x: x-0.5)
    cartData_pd_20P["Stage"] = "20P"
    
    cartData_pd_80P = pd.DataFrame(data = cartData_fromR_list[1])
    cartData_pd_80P = cartData_pd_80P.transform(lambda x: x-0.5)
    cartData_pd_80P["Stage"] = "80P"
    
    cartData_pd_concat2080P = pd.concat( [cartData_pd_20P,cartData_pd_80P]  )
    
    print(cartData_pd_concat2080P.head())
    
    def cart2pol(x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        # circVar = np.sqrt(np.sum())
        return(rho, phi)
    
    conv_rho, conv_phi = cart2pol(cartData_pd_concat2080P["X"],cartData_pd_concat2080P["Y"])
    
    conv_comp2080P_df = pd.DataFrame(data = {'rho' : conv_rho, 'phi': conv_phi, 'Stage': cartData_pd_concat2080P["Stage"]  })
    
    fmv_2080P_polar = conv_comp2080P_df.groupby("Stage").mean()
    
    conv_comp2080P_df_naD = conv_comp2080P_df.dropna(how = "any")

    # 20P data
    fmv_20P_phi = conv_comp2080P_df_naD[conv_comp2080P_df_naD['Stage'] == "20P"]['phi']
    fmv_20P_R = np.sqrt(np.sum(np.cos(fmv_20P_phi))**2 + np.sum(np.sin(fmv_20P_phi))**2)
    fmv_20P_Rav = fmv_20P_R/len(fmv_20P_phi)
    fmv_20P_CircVar = 1- fmv_20P_Rav
    
    # scipy to find circ std
    fmv_20P_rho = conv_comp2080P_df_naD[conv_comp2080P_df_naD['Stage'] == "20P"]['rho']
    fmv_20P_circular_std_dev = circstd(fmv_20P_phi)
    # print(fmv_20P_circular_std_dev*180/np.pi)
    

    # 80P data
    fmv_80P_phi = conv_comp2080P_df_naD[conv_comp2080P_df_naD['Stage'] == "80P"]['phi']
    fmv_80P_R = np.sqrt(np.sum(np.cos(fmv_80P_phi))**2 + np.sum(np.sin(fmv_80P_phi))**2)
    fmv_80P_Rav = fmv_80P_R/len(fmv_80P_phi)
    fmv_80P_CircVar = 1- fmv_80P_Rav
    
    # scipy to find circ std
    fmv_80P_rho = conv_comp2080P_df_naD[conv_comp2080P_df_naD['Stage'] == "80P"]['rho']
    fmv_80P_circular_std_dev = circstd(fmv_80P_phi)
    # print(fmv_80P_circular_std_dev*180/np.pi)
    
    # Assigning Arrow Props
    if fA_arrowSize == "S":
      fV_arrowProps_1 = dict(facecolor='C0',  
        width = 2.1,  
        headwidth = 7,  
        alpha = 0.7,  
        ec = "black",  
        )
      fV_arrowProps_2 = dict(facecolor='C1',  
        width = 2.1,  
        headwidth = 7,  
        alpha = 0.7,  
        ec = "black",  
        )
    elif fA_arrowSize == "M":
      fV_arrowProps_1 = dict(facecolor='C0',  
        width = 5.1,  
        headwidth = 10,  
        alpha = 0.7,  
        ec = "black",  
        )
      fV_arrowProps_2 = dict(facecolor='C1',  
        width = 5.1,  
        headwidth = 10,  
        alpha = 0.7,  
        ec = "black",  
        )
    else:
      fV_arrowProps_1 = dict(facecolor='C0',  
        width = 7.1,  
        headwidth = 14,  
        alpha = 0.7,  
        ec = "black",  
        )
      fV_arrowProps_2 = dict(facecolor='C1',  
        width = 7.1,  
        headwidth = 14,  
        alpha = 0.7,  
        ec = "black",  
        )
      
    
    # cleaning the axes and plots
    plt.clf()
    plt.cla()
    
    # resetting the dataframe index
    conv_comp2080P_df  = conv_comp2080P_df.reset_index()
    
    # extracting figure and axes for the new plot
    fig_p, ax_p = plt.subplots(subplot_kw={'projection': 'polar'})
   
    
    # Doing a kde plot of data
    sns.kdeplot(data = conv_comp2080P_df, 
    y = "rho", 
    x = "phi", 
    hue = "Stage", 
    fill = True, 
    levels = fA_levels,
    thresh = fA_thresh, 
    alpha = fA_alphaC )
    
    # Scatter plot
    sns.scatterplot( data = conv_comp2080P_df, y= "rho", x = "phi", s = fA_pointSize, hue = "Stage", alpha = fA_alphaP )
    
    # setting axes limits, ticks and tick labels
    ax_p.set_rmax(0.6)
    ax_p.set_rticks([0.25, 0.5], labels = ["0.5", "1.0"])  # Less radial ticks
    ax_p.set_thetagrids([0,90,180,270])
    ax_p.xaxis.set_ticks([0,np.pi/2,np.pi, 3*np.pi/2], labels = ["","","",""])
    ax_p.set_rlim(bottom=0)
    
    # turning grid on
    ax_p.grid(True)
    
    # Changing label position
    ax_p.set_rlabel_position(22.5)  # Move radial labels away from plotted line
    
    print("All good - now annotating")
    
    print(fmv_20P_Rav/2)
    # Defining arrows to the mean of the values
    ax_p.annotate('',  
        xy=(fmv_2080P_polar["phi"]["20P"],  
        fmv_20P_Rav/2 ),
        # fmv_2080P_polar["rho"]["20P"] ),
        xytext=(0,  0),
        xycoords='data',  
        arrowprops= fV_arrowProps_1)

    print("All good - ann 1 done")
    
    ax_p.annotate('',
    xy=(fmv_2080P_polar["phi"]["80P"],
    # fmv_2080P_polar["rho"]["80P"] ),  # theta, radius
    fmv_80P_Rav/2 ),  # theta, radius
                # xytext=(0.01, 0.01),    # fraction, fraction
                xytext=(0,  0),    # fraction, fraction
                # textcoords='figure fraction',
                xycoords='data',
                arrowprops=fV_arrowProps_2
                )
    print("All good - ann2 done")
    # removing axis titles (label in python)
    plt.xlabel("")
    plt.ylabel("")
    print("All good - label done")
    # Change the position of the legend to right, center
    ax_p.get_legend
    ax_p.get_legend().remove()
    
    print("All good - legend done")
    
    # Setting DPI
    plt.rcParams['figure.dpi'] = 300
    # Showing the plot
    return(fig_p)
  
def printPlot(fA_figure):
  plt.show(fA_figure)


def py_circstd(fA_angles):
  return(circstd(fA_angles))

### Copied from PythonMatplolib RMD in prudencio account
# import  matplotlib.patches as mpatches
# import matplotlib.pyplot as plt
# import seaborn as sns
# import pandas as pd
# import numpy as np
# import sys
# 
# fA_var = r.fA_var
# 
# cartData_fromR_list = fA_var
# 
# cartData_pd_20P = pd.DataFrame(data = cartData_fromR_list[0])
# cartData_pd_20P = cartData_pd_20P.transform(lambda x: x-0.5)
# cartData_pd_20P["Stage"] = "20P"
# 
# cartData_pd_80P = pd.DataFrame(data = cartData_fromR_list[1])
# cartData_pd_80P = cartData_pd_80P.transform(lambda x: x-0.5)
# cartData_pd_80P["Stage"] = "80P"
# 
# cartData_pd_concat2080P = pd.concat( 
#   [cartData_pd_20P,cartData_pd_80P]  )
# 
# 
# def cart2pol(x, y):
#   rho = np.sqrt(x**2 + y**2)
#   phi = np.arctan2(y, x)
#   # circVar = np.sqrt(np.sum(np.cos(phi))**2 + np.sum(np.sin(phi))**2)
#   return(rho, phi)
# 
# conv_rho, conv_phi  = cart2pol(cartData_pd_concat2080P["X"],cartData_pd_concat2080P["Y"])
# 
# conv_comp2080P_df = pd.DataFrame(data = {'rho' : conv_rho, 'phi': conv_phi, 'Stage': cartData_pd_concat2080P["Stage"]  })
# 
# conv_comp2080P_df_naD = conv_comp2080P_df.dropna(how = "any")
# 
# fmv_2080P_polar = conv_comp2080P_df.groupby("Stage", dropna = True).mean()
# fmv_2080P_polar["phi"]["20P"]
# 
# fmv_20P_phi = conv_comp2080P_df_naD[conv_comp2080P_df_naD['Stage'] == "20P"]['phi']
# fmv_20P_R = np.sqrt(np.sum(np.cos(fmv_20P_phi)**2) + np.sum(np.sin(fmv_20P_phi)**2))
# fmv_20P_Rav = fmv_20P_R/len(fmv_20P_phi)
# fmv_20P_CircVar = 1- fmv_20P_Rav
# 
# 
# fmv_80P_phi = conv_comp2080P_df_naD[conv_comp2080P_df_naD['Stage'] == "80P"]['phi']
# fmv_80P_R = np.sqrt(np.sum(np.cos(fmv_80P_phi)**2) + np.sum(np.sin(fmv_80P_phi)**2))
# fmv_80P_Rav = fmv_80P_R/len(fmv_80P_phi)
# fmv_80P_CircVar = 1- fmv_80P_Rav
# 
# 
# # fmv_2080P_polar["rho"]
# 
# # fmv_2080P_VarCirc = 1- conv_circVar/len(conv_rho)
# # conv_comp2080P_df.groupby('Stage')['phi'].transform(lambda x: np.sqrt(np.cos(x[~np.isnan(x)])))
# 
# # conv_comp2080P_df['R'] = conv_comp2080P_df.groupby('Stage', dropna = True)['phi'].transform(lambda x: np.sqrt((np.sum(np.cos(x)))**2 +  (np.sum(np.sin(x)))**2))
# # 
# # conv_comp2080P_df.groupby('Stage').count()
# # conv_comp2080P_df.groupby('Stage').size()
# # fmv_2080P_polar
# # 
# # conv_comp2080P_df_naD.filter('Stage' == "20P")
