"""
Copied from https://github.com/xuwkk/steady-state-power-system.
The code here may not be the same as the original repo.

A python library to generate FDI attack
Inherited from the state estimation tool run_AC_SE
Author: W XU 
"""

from class_se import SE
from pypower.api import *
import numpy as np
import warnings

class FDI(SE):
    # Inherit from the Class SE
    def __init__(self, case, noise_sigma, idx, fpr):
        super().__init__(case, noise_sigma, idx, fpr)
    
    def gen_ran_att(self, z_noise, att_ratio_max):
        """
        Generate a random attack without using the knowledge of model
        att_ratio_max: the maximum change ratio of each measurement
        """
        att_ratio = -att_ratio_max + att_ratio_max*2*np.random.rand(z_noise.shape[0])
        att_ratio = np.expand_dims(att_ratio, axis = 1)
        z_att_noise = z_noise * (1+att_ratio)

        return z_att_noise

    def gen_fdi_att(self, v_est, att_spec):
        """
        Generate a single FDI attack based on the Given state
        Only the non-reference bus can be attacked
        v_est: estimated state
        att_spec: a dictionary contains the attack information
        'ang_posi': 
        'mag_posi':
        'ang_str':
        'mag_str':
        """
        
        ang_posi = att_spec['ang_posi']
        ang_str = att_spec['ang_str']
        mag_posi = att_spec['mag_posi']
        mag_str = att_spec['mag_str']
        
        # Raise an error if the reference bus is attacked
        if self.ref_index[0] in ang_posi or self.ref_index[0] in mag_posi:
            warnings.warn('The reference bus is attacked. Consider reconfiguring the attack positions.')
        
        vang_est = np.angle(v_est)
        vmag_est = np.abs(v_est)
        
        vang_att = vang_est.copy()
        vmag_att = vmag_est.copy()
        
        vang_att[ang_posi] = vang_est[ang_posi] * (1 + ang_str)
        vmag_att[mag_posi] = vmag_est[mag_posi] * (1 + mag_str)
        
        v_att = vmag_att * np.exp(1j*vang_att)
        
        return v_att, ang_posi
    