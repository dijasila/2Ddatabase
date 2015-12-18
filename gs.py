import sys
from optparse import OptionParser
from math import pi, sqrt, acos
import pickle
import numpy as np
import ase.io
from ase.atoms import Atoms
from ase.parallel import paropen, parprint
from ase.units import Hartree, Bohr
from ase.dft.bandgap import get_band_gap
from gpaw import GPAW, PW, FermiDirac
from gpaw.symmetry import Symmetry
from gpaw.occupations import *#get_homo_lumo
from utils import get_vacuum_level

nkx = 12
nky = 12
xc = 'PBE'
ecut = 600.
vacuum = 5.
minvacuum = 10.
converge_vacuum = True
spinpol = False
tol = 0.01

"""Sweet command line interface.
Usage:
python gs.py [OPTIONS] NAME ATOMSFILE

NAME = name (determines file name),
ATOMSFILE = File containing relaxed structure.
OPTIONS:
-k N or -k NX,NY: Number of k-points either (N, N) or (NX. NY)
-x XC: XC functional
-c ECUT: PW cut-off energy
-v VAC: Force the calculation to be done with VAC vacumm = interlayer distance.
        If not set, the calculation (absolute fermi-level and band gap) will be
        converged with respect to vacuum, ensuring well-converged properties
        with minimum amoount of vacuum.
-s: Force calculation to be spin-polarized.
-n: Use non-symmorphic symmetries - for testing only.
"""
root = '/home/niflheim2/kiran/2ddatabase/'

parser = OptionParser()
parser.add_option('-k', '--nk', type='string', dest='nk')
parser.add_option('-x', '--xc', type='string', dest='xc')
parser.add_option('-c', '--ecut', type='float', dest='ecut')
parser.add_option('-v', '--vac', type='float', dest='vac')
parser.add_option('-s', '--spinpol', action='store_true', dest='spinpol',
                  default=False)
parser.add_option('-n', '--nonsymm', action='store_true',
                  dest='non_symmorphic', default=False)

opts, args = parser.parse_args()

name = args[0]
atomsfile = args[1]
filename = '%s_gs' % (xc) # name gs ####
if opts.nk:
    nks = [int(nk) for nk in opts.nk.split(',')]
    if len(nks) > 1:
        nkx = nks[0]
        nky = nks[1]
    else:
        nkx = nks[0]
        nky = nks[0]
if opts.xc:
    xc = opts.xc
if opts.ecut:
    ecut = opts.ecut
if opts.vac:
    vacuum = opts.vac
    converge_vacuum = False
if opts.spinpol:
    spinpol = opts.spinpol

atoms = ase.io.read(atomsfile)
if converge_vacuum:
    pos_av = atoms.get_positions()
    h = np.amax(pos_av[:, 2]) - np.amin(pos_av[:, 2])
    vacuum = max(2 * h + 4, minvacuum - 1.)

del atoms.constraints
atoms.center(vacuum=0.5 * vacuum, axis=2)
spos_ac = atoms.get_scaled_positions()
spos_av = atoms.get_positions()

# Relaxation may have introduced small numerical errors that leace the crystal
# slightly asymmetric - use tolerances to symmetrizes the unit cell.
cell_cv = atoms.get_cell()
a1 = cell_cv[0]
a2 = cell_cv[1]
a3 = cell_cv[2]
l1 = np.linalg.norm(a1)
l2 = np.linalg.norm(a2)
ang = 180 / pi * acos(np.dot(a1, a2) / (l1 * l2))
if abs(ang - 120) < 0.5 and abs(l1 - l2) < tol:
    # Hexagonal unit cell
    a1 = np.array([l1, 0, 0])
    a2 = np.array([- l1 / 2, l1 * sqrt(3) / 2, 0])
elif abs(ang - 60) < 0.5 and abs(l1 - l2) < tol:
    # Hexagonal unit cell
    a1 = np.array([l1, 0, 0])
    a2 = np.array([l1 / 2, l1 * sqrt(3) / 2, 0])
elif abs(ang - 90) < 0.5 and abs(l1 - l2) < tol:
    a1 = np.array([l1, 0, 0])
    a2 = np.array([0, l1, 0])
    xmin = np.amin(spos_ac[:, 0])
    ymin = np.amin(spos_ac[:, 1])
    spos_ac -= np.array([xmin, ymin, 0])
elif abs(ang - 90) < 0.5:
    a1 = np.array([l1, 0, 0])
    a2 = np.array([0, l2, 0])
    xmin = np.amin(spos_ac[:, 0])
    ymin = np.amin(spos_ac[:, 1])
    spos_ac -= np.array([xmin, ymin, 0])

cell_cv = np.array([a1, a2, a3])
atoms.set_cell(cell_cv)

if not spinpol:
    atoms.set_initial_magnetic_moments([0 for a in range(len(atoms))])

# We now symmetrize the atomic positions within some tolerance.
symm = Symmetry(atoms.get_atomic_numbers(), atoms.get_cell() / Bohr,
                tolerance=tol, symmorphic=not opts.non_symmorphic)
symm.analyze(spos_ac)
#symm.print_symmetries(sys.stdout)
spos_ac = symm.symmetrize_positions(spos_ac)
atoms.set_scaled_positions(spos_ac)


calc = GPAW(mode=PW(ecut),
            h=0.15, # This is relevant for FFT's. Maybe we should use a higher number for speed-up?
            xc=xc,
            kpts={'size': (nkx, nky, 1), 'gamma': True},
            occupations=FermiDirac(0.01), # Tight smearing for accuracy
            nbands=-20,
            convergence={'bands': -15}, # Converge some unoccupied band so the calculation can be used for band gaps
            spinpol=spinpol,
            symmetry={'point_group': True, 'time_reversal': True,
                      'symmorphic': not opts.non_symmorphic},
            txt=filename + '.txt')
atoms.set_calculator(calc)

atoms.center(vacuum=0.5 * vacuum, axis=2)

if converge_vacuum:
    E0 = atoms.get_potential_energy()
    Vvac0 = get_vacuum_level(calc) # Get the Hartree potential in vacuum
# ef0 = calc.get_fermi_level() - Vvac0 # Absolute Fermi level
    valbandef0 = calc.get_homo_lumo()[0] - Vvac0 # Absolute valence band energy level
    bandgap0, k1, k2 = get_band_gap(calc, output=None)

# If vacuum is not set, increase vacuum until the Fermi level and band gap are
# converged. This should ensure smallest possible unit cell.

while converge_vacuum:
    vacuum += 1.0
    atoms.center(vacuum=0.5 * vacuum, axis=2)
    E = atoms.get_potential_energy()
    Vvac = get_vacuum_level(calc)
    #ef = calc.get_fermi_level() - Vvac
    valband = calc.get_homo_lumo()[0]
    valbandef = valband - Vvac
    bandgap, k1, k2 = get_band_gap(calc, output=None)
    diff = bandgap - bandgap0
    parprint('vacuum=%s, Valence band level=%s, Epot change=%s, Ef change=%s, band gap change=%s' %
             (vacuum, valband, abs(E - E0), abs(valbandef - valbandef0), diff))
    if abs(E - E0) < 0.01 and abs(valbandef - valbandef0) < 0.01 and abs(diff) < 0.01:
            converge_vacuum = False

    E0 = E
    #ef0 = ef
    valbandef0 = valbandef
    bandgap0 = bandgap
    #calc.write(filename + '_vacuum_%s.gpw'%vacuum)
    #ase.io.write(filename + '.traj', atoms)

parprint('Using vacuum=%s, calculating forces' % vacuum)

# Get forces and stress to store in GPW to put in database
atoms.get_forces()
if xc != 'GLLBSC':
    # Stress doesn't work for GLLB
    parprint('Calculating stress')
    atoms.get_stress()

if xc == 'GLLBSC':
    # Calculate delta xc
    response = calc.hamiltonian.xc.xcs['RESPONSE']
    response.calculate_delta_xc()
    ksgap, dxc = response.calculate_delta_xc_perturbation()
    data = {'ksgap': ksgap, 'dxc': dxc}
    pickle.dump(data, paropen(filename + '_dxc.pckl', 'wb'))

calc.write(filename + '.gpw')
ase.io.write(filename + '.traj', atoms)
