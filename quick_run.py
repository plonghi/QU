#!/usr/bin/env python
"""
Just start computing data for the T3 graph
"""

from mcsn import MCSN 
from config import MCSNConfig 
from soliton_data import NetworkSolitonContent
cf = MCSNConfig(file_path='graph_data/T_3_revisited.ini') 
w=MCSN(branch_points=cf['branch_points'], streets=cf['streets'], 
    joints=cf['joints'], homology_classes=cf['homology_classes'])
w_solitons = NetworkSolitonContent(network=w, iterations=3)  
