from openmm import LangevinIntegrator, Platform
from openmm.app import Simulation, DCDReporter, StateDataReporter
from openmm.unit import kelvin, picoseconds, nanometer

class SimulationRunner:
    def __init__(self, topology, system, positions, production_steps):
        self.topology = topology
        self.system = system
        self.positions = positions
        self.production_steps = production_steps
        self.integrator = LangevinIntegrator(300*kelvin, 1/picoseconds, 0.002*picoseconds)
        self.simulation = None

    def run_simulation(self):
        """
        Runs the simulation, including minimization, equilibration, and production.
        """
        platform = Platform.getPlatformByName('CPU')
        self.simulation = Simulation(self.topology, self.system, self.integrator, platform)
        self.simulation.context.setPositions(self.positions)

        # Minimization
        print("Minimizing energy...")
        self.simulation.minimizeEnergy()

        # Equilibration
        print("Running equilibration...")
        self.simulation.context.setVelocitiesToTemperature(300*kelvin)
        self.simulation.step(1000)

        # Production
        print("Running production...")
        self.simulation.reporters.append(DCDReporter('production.dcd', 1000))
        self.simulation.reporters.append(StateDataReporter('production.log', 1000, step=True,
                                                           potentialEnergy=True, temperature=True, density=True))
        self.simulation.step(self.production_steps)

        print("Production run complete.")