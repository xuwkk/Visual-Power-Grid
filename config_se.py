"""
Contains the descriptions on the configuration of state estimation and optimal power flow
"""

from pypower.api import ppoption

# SE settings: output
se_config = {
    'tol' : 1e-5,       # the tolerance on the minimum jacobian matrix norm changes before considered as converged
    'max_it' : 100,     # maximum iteration
    'verbose' : 0       # description settings on the output, 0 for nothing to show; 1 for showing the norm loss in the Newton-Raphson method
}

# OPF settings: no output
opt = ppoption()              # OPF options
opt['VERBOSE'] = 0
opt['OUT_ALL'] = 0
opt['OPF_FLOW_LIM'] = 1       # Constraint on the active power flow