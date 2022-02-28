"""
Contains the descriptions on the configuration of state estimation
"""

se_config = {
    'tol' : 1e-5,       # the tolerance on the minimum jacobian matrix norm changes before considered as converged
    'max_it' : 100,     # maximum iteration
    'verbose' : 0       # description settings on the output, 0 for nothing to show; 1 for showing the norm loss in the Newton-Raphson method
}