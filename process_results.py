from pylab import *
import dnest4
import corner
from matplotlib import rcParams

import dnest4.classic as dn4
dn4.postprocess(single_precision=True)
import display
import corner
import os
import pyfits


rcParams["font.size"] = 16
rcParams["font.family"] = "sans-serif"
rcParams["font.sans-serif"] = ["Computer Modern Sans"]
rcParams["text.usetex"] = True
rcParams["text.latex.preamble"] = r"\usepackage{cmbright}"




# Piecewise linear stretch
def stretch(x):
	y = x.copy()
	y = (y - y.min())/(y.max() - y.min())
	y[y > 0.1] = 0.1 + 0.05*(y[y > 0.1] - 0.1)
	return y

saveFrames = False # For making movies
if saveFrames:
	os.system('rm Frames/*.png')
