import sys

from mcsn import MCSN
from soliton_data import NetworkSolitonContent
from config import MCSNConfig



config_file_name = 'argyres_douglas_3'

sys.stdout = open('results/soliton_data_'+config_file_name+'.txt', 'w')

cf = MCSNConfig(file_path='config/'+config_file_name+'.ini')

w = MCSN(
  branch_points=cf['branch_points'], 
  streets=cf['streets'], 
  joints=cf['joints'], 
  homology_classes=cf['homology_classes']
)

w_solitons = NetworkSolitonContent(network=w, iterations=cf['iterations'])

# #------ Finished creating network and computing soliton data-------

print '\n\n===============\nNETWORK DATA\n==============='
w.print_info()

print '\n\n===============\nSOLITON DATA\n==============='
w_solitons.print_info()
