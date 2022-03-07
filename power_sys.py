"""
The main file to run the power system steady state simulation
Inherit from the FDI class
"""

import numpy as np
from pypower.api import case14, ppoption, runopf
from pypower.idx_bus import PD, QD
from pypower.idx_brch import RATE_A, BR_X
from fdi_att import FDI
from config_mea_idx import define_mea_idx_noise
from config_se import se_config, opt
import copy
from gen_data import gen_case, gen_load

class power_env(FDI):
    def __init__(self, case, case_name, noise_sigma, idx, fpr):
        # Inherit from FDI which inherits from SE 
        super().__init__(case, noise_sigma, idx, fpr)
        
        """
        Modify the grid details
        """
        if case_name == 'case14':
            case['branch'][:,RATE_A] = case['branch'][:,RATE_A]/2.5  # Further reduce the active line flow limit
            small_load_line = [15,18,19]                             # Reduce the active line flow limit in lines with low power loading rate
            case['branch'][small_load_line,RATE_A] = 20
            case['gencost'][1:,4] = 30          # Set the non-ref generator bus linear cost
            
        # Add new attribute
        # self.case = case
        self.f_bus = self.case_int['branch'][:, 0].astype('int')        # list of "from" buses
        self.t_bus = self.case_int['branch'][:, 1].astype('int')        # list of "to" buses
        
        # Test the case
        result = self.run_opf()
        if result['success']:
            print('Initial OPF tests ok.')
        else:
            print('The OPF under this load condition is not converged! Decrease the load.')

        print('*'*60)
        
        # Load
        load_active_dir = f"src/{case_name}/load/load_active.npy"
        load_reactive_dir = f"src/{case_name}/load/load_reactive.npy"

        self.load_active = np.load(load_active_dir)
        self.load_reactive = np.load(load_reactive_dir)        

        # The default reactance in case file
        self.reactance_ori = copy.deepcopy(self.case['branch'][:,BR_X])   
    
    def run_opf(self, **kwargs):
        """
        Rewrite the run_opf function to directly read OPF index
        """
        
        case_opf = copy.deepcopy(self.case)
        if 'opf_idx' in kwargs.keys():
            print(f'Load index: {10*kwargs["opf_idx"]}.')
            # Specify the index
            load_active_ = self.load_active[int(10*kwargs['opf_idx'])]
            load_reactive_ = self.load_reactive[int(10*kwargs['opf_idx'])]
        
            case_opf['bus'][:,PD] = load_active_
            case_opf['bus'][:,QD] = load_reactive_
        else:
            # run the default
            print('Run on the default load condition.')
            pass
        
        result = runopf(case_opf, opt)
        
        return result

"""
Test
"""
if __name__ == "__main__":
    # Instance power env
    case_name = 'case14'
    case = case14()
    case = gen_case(case, 'case14')  # Modify the case
    
    # Define measurement index
    mea_idx, no_mea, noise_sigma = define_mea_idx_noise(case, 'RTU')
    
    # Instance the class
    case_env = power_env(case = case, case_name = case_name, noise_sigma = noise_sigma, idx = mea_idx, fpr = 0.05)
    
    # Generate load if it does not exist
    _, _ = gen_load(case, 'case14')
    
    # Run opf 
    result = case_env.run_opf()
    
    # Construct the measurement
    z, z_noise, vang_ref, vmag_ref = case_env.construct_mea(result) # Get the measurement
    
    # Run AC-SE
    se_config['verbose'] = 1        # Default is 0
    v_est, _ = case_env.ac_se_pypower(z_noise, vang_ref, vmag_ref, config = se_config)

    # BDD
    residual = case_env.bdd_residual(z_noise, v_est)    
    print(f'BDD threshold: {case_env.bdd_threshold}')
    print(f'residual: {residual}')

    print('*'*60)
    
    # Run OPF on a given load
    result = case_env.run_opf(opf_idx = 1)
    # Construct the measurement
    z, z_noise, vang_ref, vmag_ref = case_env.construct_mea(result) # Get the measurement
    v_est, _ = case_env.ac_se_pypower(z_noise, vang_ref, vmag_ref, config = se_config)
    residual = case_env.bdd_residual(z_noise, v_est)    
    print(f'BDD threshold: {case_env.bdd_threshold}')
    print(f'residual: {residual}')
    
