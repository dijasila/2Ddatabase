from math import pi, sqrt
import sys
import os
import optparse
import pickle
import numpy as np
import ase.io
from ase.io.trajectory import Trajectory
from ase.optimize import BFGSLineSearch, BFGS, MDMin
from ase.constraints import UnitCellFilter, StrainFilter, FixedLine
from ase.parallel import parprint, paropen
from ase.data import ground_state_magnetic_moments
from gpaw import GPAW, PW, FermiDirac
from gpaw.mpi import world
from scaledfilters import ScaledStrainFilter, ScaledUnitCellFilter

nkx = 6
nky = 6
xc = 'PBE'
ecut = 600.0
vacuum = 10.
fixxy = False
fmax = 0.01
smax = 0.01
suffix = ''
nlat = 1
strain = 0.0

"""Nice command line interface for relaxation of 2D structures.

This script works like:

python relax.py [OPTIONS] NAME ATOMSFILE

where OPTIONS are any of the options defined below, NAME is the of the material
like Phosphorene or H-MoS2 and ATOMSFILE is some file with the initial atomic
structure readable by ase.io.read(), i.e. .cif, .xyz, .traj, .gpw etc.
"""
parser = optparse.OptionParser()
parser.add_option('-k', '--nk', type='string', dest='nk')
parser.add_option('-x', '--xc', type='string', dest='xc')
parser.add_option('-c', '--ecut', type='float', dest='ecut')
parser.add_option('-v', '--vacuum', type='float', dest='vacuum')
#parser.add_option('-n', '--nlat', type='int', dest='nlat')
parser.add_option('-f', '--fmax', type='float', dest='fmax')
parser.add_option('-s', '--smax', type='float', dest='smax')
parser.add_option('-z', '--fixxy', action='store_true', default=False,
                  dest='fixxy')
parser.add_option('-u', '--useucf', action='store_true', dest='use_ucf',
                  default=False)
parser.add_option('-e', '--suffix', type='string', dest='suffix')
parser.add_option('-b', '--bandpar', action='store_true', dest='bandpar',
                  default=False)

opts, args = parser.parse_args()
name = args[0]
atomsfile = args[1]

if opts.nk:
    nks = opts.nk.split(',')
    if len(nks) == 1:
        nkx = int(nks[0])
        nky = int(nks[0])
    else:
        nkx = int(nks[0])
        nky = int(nks[1])

if opts.xc:
    xc = opts.xc
if opts.ecut:
    ecut = opts.ecut
if opts.vacuum:
    vacuum = opts.vacuum
if opts.fmax:
    fmax = opts.fmax
if opts.smax:
    smax = opts.smax
if opts.fixxy:
    fixxy = opts.fixxy
if opts.suffix:
    suffix = opts.suffix
bandpar = opts.bandpar
use_ucf = opts.use_ucf

filename = '%s_%s_relax' % (name, xc)

if suffix and suffix != '':
    filename += suffix

parprint('Relaxation calculation of %s' % name)
parprint('  atom file: %s' % atomsfile)
parprint('Parameters:')
parprint('  kpts: (%d, %d, 1)' % (nkx, nky))
parprint('  xc: %s' % xc)
parprint('  ecut: %s' % ecut)
parprint('  vacuum: %s' % vacuum)
parprint('  fmax: %s' % fmax)
parprint('  smax: %s' % smax)
parprint('  fixxy: %s' % fixxy)
parprint('  suffix: %s' % suffix)
parprint('')

atoms = ase.io.read(atomsfile)

# Some semi-arbitrary way of determining initial magnetic moments from the
# isolated compound magnetic moments. Maybe not the best way.
if np.allclose(atoms.get_initial_magnetic_moments(), 0):
    magmoms = []
    Ns = atoms.get_atomic_numbers()
    for a in range(len(atoms)):
        N = Ns[a]
        magmoms.append(0.8 * ground_state_magnetic_moments[N])
    atoms.set_initial_magnetic_moments(magmoms)
    parprint('using initial magmoms: %s' % magmoms)

atoms.set_pbc(True)
atoms.center(vacuum=vacuum, axis=2)

cell0 = atoms.get_cell()
pos0 = atoms.get_scaled_positions()
c = cell0[2, 2] # Cell height

# Get initial 2D lattice parameters and angle
vec1 = cell0[0, :]
vec2 = cell0[1, :]
vec1l = np.linalg.norm(vec1)
vec2l = np.linalg.norm(vec2)
ang = np.arccos(np.dot(vec1, vec2) / vec1l / vec2l)

# Some remnant of a method of relaxing unit cell by varying lattice constants
# and fitting to obtain the optimal method. It is not implemented anymore, but
# I'll leave the code here in case it gets relevant.
# Linearly space lattice constants between +/- strain:
pctlo = 1. - strain / 100.
pcthi = 1. + strain / 100.
if nlat == 1:
    a_list = [vec1l]
    b_list = [vec2l]
    a_step = 0
    b_step = 0
else:
    a_list = np.linspace(pctlo * vec1l, pcthi * vec1l, nlat)
    b_list = np.linspace(pctlo * vec2l, pcthi * vec2l, nlat)
    a_step = (pcthi - pctlo) * vec1l / (nlat - 1.0)
    b_step = (pcthi - pctlo) * vec2l / (nlat - 1.0)

# Determine if crystal is hexagonal, tetragonal or something else. If hexagonal
# or tetragonal we will fix the angle by setting the mask. This method is
# somewhat tolerant to a small numerical deviations.
relax2d = False
mask = [1, 1, 1, 1, 1, 1]
lat_list = []
if abs(ang - 2 * pi / 3) < 1e-3:
    bx = -1. / 2
    by = np.sqrt(3.) / 2.
    mask = [1, 1, 0, 0, 0, 0]
    if abs(vec1l - vec2l) < 1e-3:
        # Hexagonal unit cell
        lat_list = [(a, a) for a in a_list]
    else:
        lat_list = [(a, b) for a in a_list for b in b_list]
elif abs(ang - pi / 2) < 1e-3:
    bx = 0.0
    by = 1.0
    mask = [1, 1, 0, 0, 0, 0]
    if abs(vec1l - vec2l) < 1e-3:
        # Tetragonal unit cell
        lat_list = [(a, a) for a in a_list]
    else:
        lat_list = [(a, b) for a in a_list for b in b_list]
else:
    relax2d = True
    bx = np.cos(ang)
    by = np.sin(ang)
    mask = [1, 1, 0, 0, 0, 1]
    lat_list = [(a, b) for a in a_list for b in b_list]

# Initializes the cell in the high symmetry setup.
a0 = vec1l
b0 = vec2l
cell = np.array([[a0, 0, 0],
                 [bx * b0, by * b0, 0],
                 [0, 0, c]])
atoms.set_cell(cell)
atoms.set_scaled_positions(pos0)

# Be sure to remove any previous constraints
del atoms.constraints

# If fixxy we only allow the atoms to relax in the z-direction. This should
# make calculations faster (?) for simple structures with high symmetry.
if fixxy:
    zdir = np.array([0, 0, 1])
    constraints = []
    for a, atom in enumerate(atoms):
        constraints.append(FixedLine(a, zdir))
    atoms.set_constraint(constraints)

# For huge unit cells it may be necessary to use scalapack to parallelize over
# bands in the diagonalization since we have a ton of PWs.
if not bandpar:
    parallel = {'kpt': None, 'order': 'kdb'}
else:
    parallel = {'kpt': 1, 'band': world.size, 'sl_auto': True}

calc = GPAW(mode=PW(ecut), xc=xc, spinpol=True,
            kpts={'size': (nkx, nky, 1), 'gamma': True},
            occupations=FermiDirac(0.05), # Use broad smearing for convergence
            symmetry='off', # Do not use symmetries. We could allow this for structures with known symmetries.
            parallel=parallel,
            eigensolver='rmm-diis', # Allows for band parallelization.
            txt=filename + '.txt')
calc.attach(calc.write, 0, filename + '.gpw')
atoms.set_calculator(calc)

traj = Trajectory(filename + '.traj', atoms=atoms, mode='a')


if use_ucf:
    # Unitcell filter is a bad idea since the strain and atomic positions are
    # scaled differently and the force threshold thus doesn't make sense.
    restartfile = filename + '_restart.pckl'
    ucf = ScaledUnitCellFilter(atoms, mask=mask)
    dyn = BFGSLineSearch(ucf)
    
    #if not os.path.isfile(restartfile):
    #    try:
    #        replay_traj = Trajectory(filename + '.traj', 'r')
    #        dyn.replay_trajectory(replay_traj)
    #    except Exception:
    #        pass
    
    dyn.attach(traj)
    dyn.run(fmax=fmax)
else:
    # Relax strain and atomic positions independently. This is not as effecient
    # but ensures a more stable relaxation which keep symmetries and the
    # force threshold makes sense.
    # We have to use my home made ScaledStrainFilter class for this:
    strain_filter = ScaledStrainFilter(atoms, mask=mask)
    #lat_opt = BFGSLineSearch(strain_filter)
    # Optimizer for unit cell:
    lat_opt = MDMin(strain_filter) # This appears faster??
    lat_opt.attach(traj)
    # Optimizer for atomic positions:
    atom_opt = BFGSLineSearch(atoms, trajectory=traj)

    max_force = 1000
    max_stress = 1000
    # Loop until convergence. fmax and smax should really have the same size...
    while max_force >= fmax or max_stress >= smax:
        parprint('Starting optimization iteration')
        atom_opt.run(fmax=fmax)
        parprint('Atoms optimization done')
        lat_opt.run(fmax=smax)
        parprint('Lattice optimization done')
        # The max_stress is actually the maximum "average force" per atom due
        # due to stress. For definition see:
        # http://scitation.aip.org/content/aip/journal/jcp/136/7/10.1063/1.3684549
        max_stress = np.amax(strain_filter.get_forces())
        max_force = np.amax(np.sum(atoms.get_forces()**2, axis=1)**0.5)
        # max_stress and max_force are now of the same units and same order of
        # magnitude.
        parprint('max stress force: %s eV/Å, ' % max_stress +
                 'max force=%s eV/Å' % max_force)

calc.write(filename + '.gpw')
