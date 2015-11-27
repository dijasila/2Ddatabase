import numpy as np
from ase.constraints import Filter

class ScaledStrainFilter(Filter):

    def __init__(self, atoms, mask=None):
        """Create a filter applying a homogeneous strain to a list of atoms.

        The first argument, atoms, is the atoms object.

        The optional second argument, mask, is a list of six booleans,
        indicating which of the six independent components of the
        strain that are allowed to become non-zero.  It defaults to
        [1,1,1,1,1,1].

        """

        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.index = np.arange(len(atoms))
        self.n = self.index.sum()

        self.origcell = atoms.get_cell()
        vol = abs(np.linalg.det(self.origcell))
        natoms = len(self.atoms)
        self.jacobian0 = vol**(1. / 3) * natoms**(1. / 6)
    
    def get_positions(self):
        return self.strain.reshape((2, 3)).copy() * self.jacobian0

    def set_positions(self, new):
        new = new.ravel() * self.mask / self.jacobian0
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new
    """
    def get_positions(self):
        return self.strain.reshape((2, 3)).copy() * self.jacobian0

    def set_positions(self, new):
        new = new.ravel() * self.mask
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]]) / self.jacobian0

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new / self.jacobian0
    """
    def get_forces(self):
        vol = self.atoms.get_volume()
        natoms = len(self.atoms)
        jacobian = vol**(1. / 3) * natoms**(1. / 6)
        
        stress = self.atoms.get_stress()
        
        return -vol / jacobian * (stress * self.mask).reshape((2, 3))

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return 2

        
class ScaledUnitCellFilter(Filter):
    """Modify the supercell and the atom positions. """
    def __init__(self, atoms, mask=None):
        """Create a filter that returns the atomic forces and unit cell
        stresses together, so they can simultaneously be minimized.

        The first argument, atoms, is the atoms object. The optional second
        argument, mask, is a list of booleans, indicating which of the six
        independent components of the strain are relaxed.

        - True = relax to zero
        - False = fixed, ignore this component

        You can still use constraints on the atoms, e.g. FixAtoms, to control
        the relaxation of the atoms.

        >>> # this should be equivalent to the StrainFilter
        >>> atoms = Atoms(...)
        >>> atoms.set_constraint(FixAtoms(mask=[True for atom in atoms]))
        >>> ucf = UnitCellFilter(atoms)

        You should not attach this UnitCellFilter object to a
        trajectory. Instead, create a trajectory for the atoms, and
        attach it to an optimizer like this:

        >>> atoms = Atoms(...)
        >>> ucf = UnitCellFilter(atoms)
        >>> qn = QuasiNewton(ucf)
        >>> traj = Trajectory('TiO2.traj', 'w', atoms)
        >>> qn.attach(traj)
        >>> qn.run(fmax=0.05)

        Helpful conversion table:

        - 0.05 eV/A^3   = 8 GPA
        - 0.003 eV/A^3  = 0.48 GPa
        - 0.0006 eV/A^3 = 0.096 GPa
        - 0.0003 eV/A^3 = 0.048 GPa
        - 0.0001 eV/A^3 = 0.02 GPa
        """

        Filter.__init__(self, atoms, indices=range(len(atoms)))

        self.atoms = atoms
        self.strain = np.zeros(6)

        if mask is None:
            self.mask = np.ones(6)
        else:
            self.mask = np.array(mask)

        self.origcell = atoms.get_cell()
        self.copy = self.atoms.copy
        self.arrays = self.atoms.arrays

        vol = abs(np.linalg.det(self.origcell))
        natoms = len(self.atoms)
        self.jacobian0 = vol**(1. / 3) * natoms**(1. / 6)

    def get_positions(self):
        '''
        this returns an array with shape (natoms + 2,3).

        the first natoms rows are the positions of the atoms, the last
        two rows are the strains associated with the unit cell
        '''

        atom_positions = self.atoms.get_positions()
        strains = self.strain.reshape((2, 3))

        natoms = len(self.atoms)
        all_pos = np.zeros((natoms + 2, 3), np.float)
        all_pos[0:natoms, :] = atom_positions
        all_pos[natoms:, :] = strains * self.jacobian0

        return all_pos

    def set_positions(self, new):
        '''
        new is an array with shape (natoms+2,3).

        the first natoms rows are the positions of the atoms, the last
        two rows are the strains used to change the cell shape.

        The atom positions are set first, then the unit cell is
        changed keeping the atoms in their scaled positions.
        '''

        natoms = len(self.atoms)

        atom_positions = new[0:natoms, :]
        self.atoms.set_positions(atom_positions)

        new = new[natoms:, :]  # this is only the strains
        new = new.ravel() * self.mask / self.jacobian0
        eps = np.array([[1.0 + new[0], 0.5 * new[5], 0.5 * new[4]],
                        [0.5 * new[5], 1.0 + new[1], 0.5 * new[3]],
                        [0.5 * new[4], 0.5 * new[3], 1.0 + new[2]]])

        self.atoms.set_cell(np.dot(self.origcell, eps), scale_atoms=True)
        self.strain[:] = new

    def get_forces(self, apply_constraint=False):
        '''
        returns an array with shape (natoms+2,3) of the atomic forces
        and unit cell stresses.

        the first natoms rows are the forces on the atoms, the last
        two rows are the stresses on the unit cell, which have been
        reshaped to look like "atomic forces". i.e.,

        f[-2] = -vol*[sxx,syy,szz]*mask[0:3]
        f[-1] = -vol*[syz, sxz, sxy]*mask[3:]

        apply_constraint is an argument expected by ase
        '''

        stress = self.atoms.get_stress()
        atom_forces = self.atoms.get_forces()

        vol = self.atoms.get_volume()
        natoms = len(self.atoms)
        all_forces = np.zeros((natoms + 2, 3), np.float)
        all_forces[0:natoms, :] = atom_forces

        jacobian = vol**(1. / 3) * natoms**(1. / 6)
        stress_forces = -vol / jacobian * (stress * self.mask).reshape((2, 3))
        all_forces[natoms:, :] = stress_forces
        return all_forces

    def get_potential_energy(self):
        return self.atoms.get_potential_energy()

    def has(self, x):
        return self.atoms.has(x)

    def __len__(self):
        return (2 + len(self.atoms))

        
