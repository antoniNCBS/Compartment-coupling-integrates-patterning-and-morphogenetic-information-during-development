from  matplotlib.patches import ArrowStyle
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from matplotlib.backends.backend_svg import FigureCanvasSVG as FCsvg



def savePlotSVG(fA_figure, fA_filename):
  FCsvg(fA_figure).print_svg(fA_filename)
  
