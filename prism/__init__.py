# prism/__init__.py
"""
PRISM - Protein Receptor Interaction Simulation Modeler

A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
"""

__version__ = "1.0.0"
__author__ = "PRISM Development Team"

from .builder import PRISMBuilder

__all__ = ["PRISMBuilder"]

# ===========================
# prism/forcefield/__init__.py
"""
Force field generators for PRISM
"""

from forcefield.base import ForceFieldGeneratorBase
from forcefield.gaff import GAFFForceFieldGenerator
from forcefield.openff import OpenFFForceFieldGenerator

__all__ = ["ForceFieldGeneratorBase", "GAFFForceFieldGenerator", "OpenFFForceFieldGenerator"]

# ===========================
# prism/utils/__init__.py
"""
Utility functions for PRISM
"""

# Future utilities can be added here