from openmm.app import ForceField, Modeller, PME, HBonds
from openmm import Vec3
from openmm.unit import nanometers

class SystemSetup:
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.forcefield = ForceField('amber14-all.xml', 'amber14/tip3p.xml')
        self.modeller = Modeller(self.pdb_file.topology, self.pdb_file.positions)
        self.system = None

    def setup_system(self):
        """
        Sets up the simulation system by adding a water box and ions.
        """
        self.modeller.addSolvent(self.forcefield, model='tip3p', padding=1.0*nanometers)
        self.system = self.forcefield.createSystem(self.modeller.topology, nonbondedMethod=PME,
                                                   nonbondedCutoff=1.0*nanometers, constraints=HBonds,
                                                   ewaldErrorTolerance=0.0005)
        print("System setup complete.")
        return self.system, self.modeller.topology, self.modeller.positions