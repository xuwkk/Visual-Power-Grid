"""
The main file to run the power system steady state simulation
"""

import numpy as np
from os.path import exists
import pandas as pd
from pypower.api import case14, ppoption, runopf
from pypower.idx_bus import PD, QD
from pypower.idx_brch import RATE_A
from run_AC_SE import SE
from fdi_att import FDI
from mea_idx import define_mea_idx_noise
from se_config import se_config
import copy


"""
Construct the power system environment based on SE and FDI.
"""

class power_env(FDI):
    def __init__(self, case, case_name, noise_sigma, idx, fpr):
        # Inherit from FDI which inherits from SE 
        super().__init__(case, noise_sigma, idx, fpr)
        
        """
        Modify the grid details
        """
        
        if case_name == 'case14':
            # Modify the load
            for i in range(len(case['bus'])):
                if case['bus'][i,PD] == 0 and i != self.ref_index:
                    # Add load on the zero load non ref bus
                    case['bus'][i,PD] = 30
                    case['bus'][i,QD] = 10
                if case['bus'][i,PD] <= 15:
                    # increase the load on small loaded bus
                    case['bus'][i,PD] = case['bus'][i,PD] *2.
            
            # Add constraints on the branch power flow RATE_A
            case['branch'][:,RATE_A] = case['branch'][:,RATE_A]/250  # The default is 9900, we change it to 99
            small_load_line = [10,16,17,18,19]                       # lines with low power loading rate
            case['branch'][small_load_line,RATE_A] = 20            
            
            # Change the cost of the generators
            case['gencost'][2:,5] = 30
            
            # Add as attribute
            self.case = case
            self.f_bus = self.case_int['branch'][:, 0].astype('int')        # list of "from" buses
            self.t_bus = self.case_int['branch'][:, 1].astype('int')        # list of "to" buses
        
        # Test the case
        result = self.run_opf()
        if result['success']:
            print('Initial OPF tests ok.')
        else:
            print('The OPF under this load condition is not converged! Decrease the load.')
    
        """
        Generate the Load Dataset
        """
        # Load address
        load_active_dir = f"src/load/{case_name}/load_active.npy"
        load_reactive_dir = f"src/load/{case_name}/load_reactive.npy"
        
        if exists(load_active_dir):
            # Load data
            self.load_active = np.load(load_active_dir)
            self.load_reactive = np.load(load_reactive_dir)      
            
        else:
            # Generate load
            load_raw_dir = "src\load\load_normalize_clean.csv"
            self._generate_load(load_raw_dir, load_active_dir, load_reactive_dir)
    
    def run_opf(self, **kwargs):
        # Rewrite the run_opf function
        
        # OPF settings
        opt = ppoption()              # OPF options
        opt['VERBOSE'] = 0
        opt['OUT_ALL'] = 0
        opt['OPF_FLOW_LIM'] = 1       # Constraint on the active power flow
        
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
    
    def _generate_load(self,load_raw_dir, load_active_dir, load_reactive_dir):
        """
        The function that is used to generate the load data.
        Define the active load wisely in the self.case['bus'][:,PD] so that the generated load is not raising non-convergency in OPF
        """
        load_raw = pd.read_csv(load_raw_dir)
        load_active = []
        load_reactive = []
        for i in range(len(self.no_bus)):
            # Active power
            load_active_ = load_raw.iloc[:,i+1].values
            load_active_max_ = np.max(load_active_)
            load_active_max = self.case['bus'][i,PD]     # rated power in the case file
            print(load_active_max)
            # Rescale
            load_active_ = load_active_ * load_active_max/load_active_max_
            load_active.append(load_active_)
            # Reactive power
            pf = 0.9+2*(1-0.95)*np.random.rand(len(load_raw)) # PF: 0.9-1.0
            Q_P = np.tan(np.arccos(pf))                       # Q to P ratio
            load_reactive_= load_active_*Q_P
            
            load_reactive.append(load_reactive_)
        
        load_active = np.array(load_active)
        load_reactive = np.array(load_reactive)
        
        np.save(load_active_dir, load_active.T, allow_pickle=True)
        np.save(load_reactive_dir, load_reactive.T, allow_pickle=True)     # (sample,bus)
        
        self.load_active = load_active
        self.load_reactive = load_reactive
        
# class power_env:
#     def __init__(self, case_name):
        
#         """
#         OPF OPTION
#         """
#         opt = ppoption()
#         opt['VERBOSE'] = 0
#         opt['OUT_ALL'] = 0
#         opt['OPF_FLOW_LIM'] = 1  # Constraint on the active power flow
#         small_load_line = [10,16,17,18,19]  # lines with low power loading rate
        
#         """
#         GRID
#         """
#         # Load and Modify the case file
#         if case_name == "case14":
#             case = case14()
#             # Determine the bus type
#             case_int = ext2int(case)  # Convert to internal representation
#             slack_index, PV_index, PQ_index = bustypes(case_int['bus'], case_int['gen'])

#             print(f'slack_index: {slack_index}')
#             print(f'PV_index: {PV_index}')
#             print(f'PQ_index: {PQ_index}')

#             # Modify the load
#             for i in range(len(case['bus'])):
#                 if case['bus'][i,PD] == 0 and i != slack_index:
#                     case['bus'][i,PD] = 30
#                     case['bus'][i,QD] = 10
#                 if case['bus'][i,PD] <= 15:
#                     case['bus'][i,PD] = case['bus'][i,PD] *3
            
#             # Add constraints on the branch power flow RATE_A
            
#             case['branch'][:,RATE_A] = case['branch'][:,RATE_A]/250  # The default is 9900, we change it to 99
#             case['branch'][small_load_line,RATE_A] = 20            
            
#             # Change the cost of the generators
#             case['gencost'][2:,5] = 30

#         # Test the case
#         result = runopf(case, opt)
#         if result['success']:
#             print('Initial OPF tests ok.')
#         else:
#             print('The OPF under this load condition is not converged! Decrease the load.')
        
#         """
#         LOAD DATASET
#         """
        
#         # Load dataset
#         load_active_dir = f"src/load/{case_name}/load_active.npy"
#         load_reactive_dir = f"src/load/{case_name}/load_reactive.npy"
        
#         if exists(load_active_dir):
#             # Load data
#             load_active = np.load(load_active_dir)
#             load_reactive = np.load(load_reactive_dir)      
            
#         else:
#             # Generate load
#             load_raw_dir = "src\load\load_normalize_clean.csv"
#             load_raw = pd.read_csv(load_raw_dir)
#             load_active = []
#             load_reactive = []
#             for i in range(len(case['bus'])):
#                 # Active power
#                 load_active_ = load_raw.iloc[:,i+1].values
#                 load_active_max_ = np.max(load_active_)
#                 load_active_max = case['bus'][i,PD]     # rated power in the case file
#                 print(load_active_max)
#                 # Rescale
#                 load_active_ = load_active_ * load_active_max/load_active_max_
#                 load_active.append(load_active_)
                
#                 # Reactive power
#                 pf = 0.9+2*(1-0.95)*np.random.rand(len(load_raw)) # PF: 0.9-1.0
#                 Q_P = np.tan(np.arccos(pf))                       # Q to P ratio
#                 load_reactive_= load_active_*Q_P
               
#                 load_reactive.append(load_reactive_)
            
#             load_active = np.array(load_active)
#             load_reactive = np.array(load_reactive)
            
#             np.save(load_active_dir, load_active.T, allow_pickle=True)
#             np.save(load_reactive_dir, load_reactive.T, allow_pickle=True)     # (sample,bus)
        
#         """
#         INCIDENCE MATRIX
#         """
        
#         # Incidence matrix: bus->branch
#         f_bus = np.array(case_int['branch'][:,F_BUS], dtype='int')
#         t_bus = np.array(case_int['branch'][:,T_BUS], dtype=int)
#         #print(f_bus)
#         incidence_matrix = np.zeros((len(case['branch']),len(case['bus'])))
#         for i in range(len(case['branch'])):
#             incidence_matrix[i,f_bus[i]] = 1
#             incidence_matrix[i,t_bus[i]] = -1
        
#         # Incidence matrix: gen->bus
#         gen_idx = (case['gen'][:,0]).astype('int')
#         print(gen_idx)
#         gen_incidence_matrix = np.zeros((len(case['bus']), len(case['gen'])))
#         for i in range(len(case['gen'])):
#             gen_incidence_matrix[gen_idx[i],i] = 1
        
#         """
#         SELF
#         """
        
#         self.case_name = case_name
#         self.case = case
#         self.load_active = load_active
#         self.load_reactive = load_reactive
#         self.opt = opt
#         self.incidence_mx = incidence_matrix
#         self.f_bus_index = f_bus
#         self.t_bus_index = t_bus
#         self.gen_incidence_matrix = gen_incidence_matrix

#     def run_opf(self, idx):
#         load_active_ = self.load_active[int(10*idx)]
#         load_reactive_ = self.load_reactive[int(10*idx)]
        
#         self.case['bus'][:,PD] = load_active_
#         self.case['bus'][:,QD] = load_reactive_
        
#         results = runopf(self.case, self.opt)
#         flag = results['success']
        
#         res_brh = results['branch']
#         res_gen = results['gen']
#         res_bus = results['bus']

#         pf = res_brh[:,PF]
#         pt = res_brh[:,PT]
#         qf = res_brh[:,QF]
#         qt = res_brh[:,QT]
        
#         pl = np.sum(np.abs(pf+pt))
        
#         # Construct the measurement for state estimation
#         pi = self.gen_incidence_matrix@res_gen[:,1] - res_bus[:,2]  # active power injection
#         qi = self.gen_incidence_matrix@res_gen[:,2] - res_bus[:,3]  # reactive power injection
        
#         z = np.concatenate([pi, pf, qi, qf], axis = -1)/self.case['baseMVA']
        
#         # Add noise on the measurement
#         sigma = 0.005
#         R = sigma**2 * np.diag(np.ones(len(z)))
#         z_noise = z + np.random.multivariate_normal(mean = np.zeros((len(z),)), cov = R)
#         print(z_noise)
#         print(R)
#         self.z_noise = z_noise
        
#         return load_active_, load_reactive_, pf, pl, results['success']
    
#     def run_ac_se(self):
#         """
#         AC State Estimation
#         """
#         # Construct the measurement
        


# Test
if __name__ == "__main__":
    case_name = 'case14'
    case = case14()
    
    # Define measurement index
    mea_idx, no_mea, noise_sigma = define_mea_idx_noise(case, 'HALF_RTU')
    # Instance the class
    case_env = power_env(case = case, case_name = case_name, noise_sigma = noise_sigma, idx = mea_idx, fpr = 0.02)
    # Run opf
    print(case_env.run_opf)   
    result = case_env.run_opf()
    # Construct the measurement
    z, z_noise, vang_ref, vmag_ref = case_env.construct_mea(result) # Get the measurement
    # Run AC-SE
    se_config['verbose'] = 1        # Default is 0
    v_est = case_env.ac_se_pypower(z_noise, vang_ref, vmag_ref, config = se_config)
    residual = case_env.bdd_residual(z_noise, vang_ref, vmag_ref)    
    print(f'BDD threshold: {case_env.bdd_threshold}')
    print(f'residual: {residual}')
    
    # 
    result = case_env.run_opf(opf_idx = 1)
    # Construct the measurement
    z, z_noise, vang_ref, vmag_ref = case_env.construct_mea(result) # Get the measurement
    v_est = case_env.ac_se_pypower(z_noise, vang_ref, vmag_ref, config = se_config)
    residual = case_env.bdd_residual(z_noise, vang_ref, vmag_ref)    
    print(f'BDD threshold: {case_env.bdd_threshold}')
    print(f'residual: {residual}')
    
