"""
Examples of measurement index and measurement noise
The measurement should always be arranged as z = [pf, pt, pi, vang, qf, qt, qi, vmag]
"""

import numpy as np
import random

def define_mea_idx_noise(case, choice):
    
        
    """
    noise setting
    """
    mea_noise_std = {}
    rtu_noise_std = 0.005
    pmu_noise_std = 0.001
    mea_noise_std['pf'] = rtu_noise_std
    mea_noise_std['pt'] = rtu_noise_std
    mea_noise_std['pi'] = rtu_noise_std
    mea_noise_std['vang'] = pmu_noise_std
    mea_noise_std['qf'] = rtu_noise_std
    mea_noise_std['qt'] = rtu_noise_std
    mea_noise_std['qi'] = rtu_noise_std
    mea_noise_std['vmag'] = pmu_noise_std
    
    """
    idx
    """
    idx = {}
    
    if choice == 'FULL':
        """
        Full Measurement
        """
        idx['pf'] = np.arange(len(case['branch']))
        idx['pt'] = np.arange(len(case['branch']))
        idx['pi'] = np.arange(len(case['bus']))
        idx['vang'] = np.arange(len(case['bus']))
        idx['qf'] = np.arange(len(case['branch']))
        idx['qt'] = np.arange(len(case['branch']))
        idx['qi'] = np.arange(len(case['bus']))
        idx['vmag'] = np.arange(len(case['bus']))

    if choice == 'RTU':
        """
        Only the RTU measurement
        """
        idx['pf'] = np.arange(len(case['branch']))
        idx['pt'] = np.arange(len(case['branch']))
        idx['pi'] = np.arange(len(case['bus']))
        idx['vang'] = []                             # No voltage phase angle measurement
        idx['qf'] = np.arange(len(case['branch']))
        idx['qt'] = np.arange(len(case['branch']))
        idx['qi'] = np.arange(len(case['bus']))
        idx['vmag'] = np.arange(len(case['bus']))
    
    if choice == 'RTU_POWER':
        """
        All power injections and flows
        """
        idx['pf'] = np.arange(len(case['branch']))
        idx['pt'] = np.arange(len(case['branch']))
        idx['pi'] = np.arange(len(case['bus']))
        idx['vang'] = []                             # No voltage phase angle measurement
        idx['qf'] = np.arange(len(case['branch']))
        idx['qt'] = np.arange(len(case['branch']))
        idx['qi'] = np.arange(len(case['bus']))
        idx['vmag'] = []

    if choice == 'HALF_RTU':
        """
        Only pf, pi, qf, qi
        """

        idx['pf'] = np.arange(len(case['branch']))
        idx['pt'] = []
        idx['pi'] = np.arange(len(case['bus']))
        idx['vang'] = []                             # No voltage phase angle measurement
        idx['qf'] = np.arange(len(case['branch']))
        idx['qt'] = []
        idx['qi'] = np.arange(len(case['bus']))
        idx['vmag'] = []
    
    if choice == 'RANDOM':
        """
        Randomly pick the measurement, may not be observable 
        """
        idx['pf'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])))
        idx['pt'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])))
        idx['pi'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])))
        idx['vang'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])))
        idx['qf'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])))
        idx['qt'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])))
        idx['qi'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])))
        idx['vmag'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])))        
    
    
    if choice == 'UNOBSERVABLE':
        """
        A case that likely gives out unobservable state estimation
        """
        idx['pf'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])-5))
        idx['pt'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])-5))
        idx['pi'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])-5))
        idx['vang'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])-5))
        idx['qf'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])-5))
        idx['qt'] = random.sample(range(len(case['branch'])), np.random.randint(len(case['branch'])-5))
        idx['qi'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])-5))
        idx['vmag'] = random.sample(range(len(case['bus'])), np.random.randint(len(case['bus'])-5)) 
        
    """
    no_mea
    """
    no_mea = 0
    for key in idx.keys():
        no_mea = no_mea + len(idx[key])

    """
    Noise vector
    """
    noise_sigma = []
    for key in idx.keys():
        noise_sigma = noise_sigma + len(idx[key])*[mea_noise_std[key]]
    
    return idx, no_mea, np.array(noise_sigma)