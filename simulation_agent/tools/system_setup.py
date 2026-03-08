"""OpenMM system setup utilities for solvated protein systems."""

from openmm.app import (
    CutoffNonPeriodic,
    CutoffPeriodic,
    ForceField,
    HBonds,
    Modeller,
    NoCutoff,
    PME,
)
from openmm.unit import nanometers


NONBONDED_METHODS = {
    "PME": PME,
    "NoCutoff": NoCutoff,
    "CutoffNonPeriodic": CutoffNonPeriodic,
    "CutoffPeriodic": CutoffPeriodic,
}

CONSTRAINT_METHODS = {
    "HBonds": HBonds,
    "None": None,
}


class SystemSetup:
    """Build a solvated OpenMM system from an input PDB structure."""

    def __init__(self, pdb_file, options=None):
        """Initialize force field and setup options for system construction."""
        options = options or {}
        self.pdb_file = pdb_file
        forcefield_files = options.get("forcefield_files", ["amber14-all.xml", "amber14/tip3p.xml"])
        if not isinstance(forcefield_files, (list, tuple)) or not forcefield_files:
            raise ValueError("'forcefield_files' must be a non-empty list.")

        self.water_model = options.get("water_model", "tip3p")
        self.padding_nm = float(options.get("padding_nm", 1.0))
        self.nonbonded_cutoff_nm = float(options.get("nonbonded_cutoff_nm", 1.0))
        self.ewald_error_tolerance = float(options.get("ewald_error_tolerance", 0.0005))
        self.nonbonded_method_name = options.get("nonbonded_method", "PME")
        self.constraints_name = options.get("constraints", "HBonds")

        if self.nonbonded_method_name not in NONBONDED_METHODS:
            supported_methods = ", ".join(sorted(NONBONDED_METHODS))
            raise ValueError(
                f"Unsupported nonbonded_method '{self.nonbonded_method_name}'. "
                f"Supported values: {supported_methods}."
            )
        if self.constraints_name not in CONSTRAINT_METHODS:
            supported_constraints = ", ".join(sorted(CONSTRAINT_METHODS))
            raise ValueError(
                f"Unsupported constraints '{self.constraints_name}'. "
                f"Supported values: {supported_constraints}."
            )

        self.forcefield = ForceField(*forcefield_files)
        self.modeller = Modeller(self.pdb_file.topology, self.pdb_file.positions)
        self.system = None

    def setup_system(self):
        """
        Add solvent and build an OpenMM System object.
        """
        self.modeller.addSolvent(
            self.forcefield,
            model=self.water_model,
            padding=self.padding_nm * nanometers,
        )
        self.system = self.forcefield.createSystem(
            self.modeller.topology,
            nonbondedMethod=NONBONDED_METHODS[self.nonbonded_method_name],
            nonbondedCutoff=self.nonbonded_cutoff_nm * nanometers,
            constraints=CONSTRAINT_METHODS[self.constraints_name],
            ewaldErrorTolerance=self.ewald_error_tolerance,
        )
        print("System setup complete.")
        return self.system, self.modeller.topology, self.modeller.positions
