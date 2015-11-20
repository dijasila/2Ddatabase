import numpy as np
from ase.units import Hartree

def get_vacuum_level(calc):
    calc.restore_state()
    vHt_g = calc.hamiltonian.pd3.ifft(calc.hamiltonian.vHt_q) * Hartree
    vHt_z = np.mean(np.mean(vHt_g, axis=0), axis=0)
    return vHt_z[0]
