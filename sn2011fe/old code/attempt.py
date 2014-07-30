import matplotlib.pyplot as plt
import numpy as np
import glob
import sqlite3 as sq3
from scipy import interpolate as intp
import math
from astropy.table import Table

class supernova(object):
    """Attributes can be added"""

SN_Array = []
compare_spectrum = []

spec1=Table.read('Data/data/sn2011fe-visit6-fuv.flm',format='ascii')
spec2=Table.read('Data/data/sn2011fe-visit6-uv.flm',format='ascii')

wave1=spec1["col1"]
flux1=spec1["col2"]


wave2=spec2["col1"]
flux2=spec2["col2"]

print max(wave1), min(wave2)