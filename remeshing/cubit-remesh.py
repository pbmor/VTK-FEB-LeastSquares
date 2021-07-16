import os
import sys
sys.path.append('/Applications/Coreform-Cubit-2020.2.app/Contents/MacOS/')
import cubit
cubit.init(['cubit','-nojournal'])
filename = os.path.abspath(os.getcwd()) + '/test2.stl'
cubit.cmd('import stl "'+filename+ '" feature_angle 135.00 merge')
cubit.cmd('surface 1 size auto factor 1')
cubit.cmd('mesh surface 1')
cubit.cmd('set exodus netcdf4 off')
cubit.cmd('set large exodus file on')
cubit.cmd('block 1 add surface 1')
filename = os.path.abspath(os.getcwd()) + '/test-remesh2.e'
cubit.cmd('export mesh "'+filename+'"  dimension 3  overwrite ')
