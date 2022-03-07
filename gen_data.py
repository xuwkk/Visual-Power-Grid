"""
This module is used to 
1. Modify the case
2. generate the load conditions for a given system configuration
"""
import pandas as pd
import numpy as np
import copy
from pypower.idx_bus import PD, QD
from pypower.idx_brch import RATE_A
from pypower.api import ext2int, bustypes
from os.path import exists

np.random.seed(10)

def gen_case(case, case_name):
    """
    The function that is used to generate the load data.
    You have to weite in this section by yourself if the case is not listed.
    In general, define the active load wisely in the self.case['bus'][:,PD] so that the generated load is not raising non-convergency in OPF

    case: an initiated case to be modified by some previous defined settings.
    Return: a new modified case.
    """
    case = copy.copy(case)
    
    # Determine the grid details
    # TODO: Merge this part into the SE
    
    # Numbers
    no_bus = len(case['bus'])
    # Bus types
    case_int = ext2int(case)
    # Determine the bus type
    ref_index, pv_index, pq_index = bustypes(case_int['bus'], case_int['gen'])  # reference bus (slack bus), pv bus, and pq (load bus)
    non_ref_index = list(pq_index) + list(pv_index)                        # non reference bus
    non_ref_index.sort()

    if case_name == 'case14':
        """
        case14:
            1. Add load on the zero load non ref bus;
            2. Increase the load on small loaded bus;
            3. Change into linear cost
        """

        # Modify the load
        for i in range(len(case['bus'])):
            if case['bus'][i,PD] == 0 and i != ref_index:
                # Add load on the zero load non ref bus
                case['bus'][i,PD] = 30
                case['bus'][i,QD] = 10     # An arbitraty value
            if case['bus'][i,PD] <= 15:
                # Increase the load on small loaded bus
                case['bus'][i,PD] = case['bus'][i,PD] * 2.
        
        # Add constraints on the branch power flow RATE_A
        case['branch'][:,RATE_A] = case['branch'][:,RATE_A]/100  # The default is 9900           
        
        # Change the cost of the generators: a linear cost is considered
        case['gencost'][:,3] = 2   # Linear cost
        case['gencost'][:,4] = case['gencost'][:,5]  # Linear cost
        case['gencost'][:,5] = case['gencost'][:,6]  # Constant cost
    
    return case

def gen_load(case, case_name):
    """
    case: is a modified case file from, e.g. gen_case
    """

    # Load address
    load_active_dir = f"src/{case_name}/load/load_active.npy"
    load_reactive_dir = f"src/{case_name}/load/load_reactive.npy"

    # Test if the data already exists.
    if exists(load_active_dir):
        # Load data
        print(f'Load file found, loading the data...')
        load_active = np.load(load_active_dir)
        load_reactive = np.load(load_reactive_dir)      
        
    else:
        # Generate load
        print(f'No load data found, generating load...')
        load_raw_dir = "src/load_normalize_clean.csv"
        load_raw = pd.read_csv(load_raw_dir)
        load_active = []
        load_reactive = []
        for i in range(len(case['bus'])):
            # Active power
            load_active_ = load_raw.iloc[:,i+1].values
            load_active_max_ = np.max(load_active_)
            load_active_max = case['bus'][i,PD]     # rated power in the case file
            print(f'Defualt load in case file:" {load_active_max}')
            # Rescale
            load_active_ = load_active_ * load_active_max/load_active_max_
            load_active.append(load_active_)
            # Reactive power
            pf = 0.9+2*(1-0.95)*np.random.rand(len(load_raw)) # PF: 0.9-1.0
            Q_P = np.tan(np.arccos(pf))                       # Q to P ratio
            load_reactive_= load_active_*Q_P
            
            load_reactive.append(load_reactive_)
        
        load_active = np.array(load_active).T
        load_reactive = np.array(load_reactive).T
        
        np.save(load_active_dir, load_active, allow_pickle=True)
        np.save(load_reactive_dir, load_reactive, allow_pickle=True)     # (sample,bus)
    
    return load_active, load_reactive