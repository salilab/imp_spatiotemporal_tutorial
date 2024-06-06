#!/bin/bash                         
# Script to write GMMs for each protein density. 
# Uses 40 Gaussians for each time-dependent ET Map
# See imp/modules/isd/pyext/src/create_gmm.py for all options.

# 0min --------------------------------------------------------------------------------------------------
python /.../imp/modules/isd/pyext/src/create_gmm.py 0min_fitted.mrc 40 0min_gmm.mrc.txt -s 0.001 -m 0min_gmm.mrc -i 200 -a 5.0

# 1min --------------------------------------------------------------------------------------------------
python /.../imp/modules/isd/pyext/src/create_gmm.py 1min_fitted.mrc 40 1min_gmm.mrc.txt -s 0.001 -m 1min_gmm.mrc -i 200 -a 5.0

# 2min --------------------------------------------------------------------------------------------------
python /.../imp/modules/isd/pyext/src/create_gmm.py 2min_fitted.mrc 40 2min_gmm.mrc.txt -s 0.001 -m 2min_gmm.mrc -i 200 -a 5.0

