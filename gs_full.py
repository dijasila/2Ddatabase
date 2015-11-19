import os.path
from math import pi
import numpy as np
from optparse import OptionParser
from ase.units import Hartree, Bohr
from gpaw import GPAW

nbecut = None

parser = OptionParser()
parser.add_option('-b', '--nbecut', type='float', dest='nbecut')
parser.add_option('-k', '--nk', type='string', dest='nk')
parser.add_option('-s', '--nonsym', action='store_true',
                  dest='non_symmorphic', default=False)
parser.add_option('-p', '--scalapack', action='store_true', dest='scalapack',
                  default=False)

opts, args = parser.parse_args()
gsfile = args[0]
if opts.nbecut:
    nbecut = opts.nbecut

filename = '.'.join(os.path.basename(gsfile).split('.')[:-1])

calc = GPAW(gsfile, h=0.15, fixdensity=True,
            txt=filename + '_full.txt',
            parallel={'band': 1},
            symmetry={'point_group': True, 'time_reversal': True,
                      'symmorphic': not opts.non_symmorphic})

if opts.nk is not None:
    nks = [int(x) for x in opts.nk.split(',')]
    if len(nks) == 1:
        nkx = nks[0]
        nky = nks[0]
    else:
        nkx = nks[0]
        nky = nks[1]
    calc.set(kpts={'size': (nkx, nky, 1), 'gamma': True})

calc.get_potential_energy()

nbands = None
if nbecut is not None:
    atoms = calc.get_atoms()
    vol = atoms.get_volume() / Bohr**3
    #vol = abs(np.linalg.det(calc.wfs.gd.cell_cv))
    nbands = int(vol * (nbecut / Hartree)**1.5 * 2**0.5 / 3 / pi**2)
    nbands += nbands % 2 # Make nbands an even number - easier to paralellize over

scalapack = None
if opts.scalapack:
    scalapack = opts.scalapack

calc.diagonalize_full_hamiltonian(nbands=nbands, scalapack=scalapack)
calc.write(filename + '_full.gpw', mode='all')
