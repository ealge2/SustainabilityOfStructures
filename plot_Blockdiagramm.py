import struct_analysis  # file with code for structural analysis
import struct_optimization  # file with code for structural optimization
import sqlite3  # import modul for SQLite
import random
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
from shapely.geometry import Polygon

# plot Blockdiagrams of the composition of emission of cross-sections
# ----------------------------------------------------------------------------------------------------------------------

def plot_Blockdiagramm(lengths, database_name, criteria, optima, floorstruc, requirements, crsec_type, mat_names,
                 g2k=0.75, qk=2.0, max_iter=100, idx_vrfctn=-1):
