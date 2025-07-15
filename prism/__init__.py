# # # prism/__init__.py
# # """
# # PRISM - Protein Receptor Interaction Simulation Modeler

# # A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
# # """

# # __version__ = "1.0.0"
# # __author__ = "PRISM Development Team"

# # from .builder import PRISMBuilder

# # __all__ = ["PRISMBuilder"]

# # # ===========================
# # # prism/forcefield/__init__.py
# # """
# # Force field generators for PRISM
# # """

# # from forcefield.base import ForceFieldGeneratorBase
# # from forcefield.gaff import GAFFForceFieldGenerator
# # from forcefield.openff import OpenFFForceFieldGenerator

# # __all__ = ["ForceFieldGeneratorBase", "GAFFForceFieldGenerator", "OpenFFForceFieldGenerator"]

# # # ===========================
# # # prism/utils/__init__.py
# # """
# # Utility functions for PRISM
# # """

# # # Future utilities can be added here


# #!/usr/bin/env python3
# # -*- coding: utf-8 -*-

# """
# PRISM - Protein Receptor Interaction Simulation Modeler

# A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
# """

# __version__ = "1.0.0"
# __author__ = "PRISM Development Team"

# from .builder import PRISMBuilder
# from .core import PRISMSystem

# # High-level API functions
# def system(protein_path, ligand_path, config=None, **kwargs):
#     """
#     Create a protein-ligand system for MD simulation.
    
#     Parameters:
#     -----------
#     protein_path : str
#         Path to protein PDB file
#     ligand_path : str  
#         Path to ligand file (MOL2/SDF)
#     config : str or dict, optional
#         Configuration file path or configuration dictionary
#     **kwargs : optional
#         Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)
    
#     Returns:
#     --------
#     PRISMSystem
#         A system object with methods to build and run simulations
        
#     Examples:
#     ---------
#     >>> import prism as pm
#     >>> system = pm.system("protein.pdb", "ligand.mol2")
#     >>> system.build()
#     >>> system.generate_mdp_files()
    
#     >>> # With custom configuration
#     >>> system = pm.system("protein.pdb", "ligand.sdf", 
#     ...                    config="config.yaml",
#     ...                    ligand_forcefield="openff")
#     >>> system.build()
#     """
#     return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)

# def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
#     """
#     Build a complete protein-ligand system (one-step function).
    
#     Parameters:
#     -----------
#     protein_path : str
#         Path to protein PDB file
#     ligand_path : str
#         Path to ligand file (MOL2/SDF)  
#     output_dir : str
#         Output directory for generated files
#     **kwargs : optional
#         Additional parameters
        
#     Returns:
#     --------
#     str
#         Path to the output directory
        
#     Examples:
#     ---------
#     >>> import prism as pm
#     >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")
#     """
#     system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
#     return system_obj.build()

# # Export main classes and functions
# __all__ = [
#     "PRISMBuilder", 
#     "PRISMSystem",
#     "system", 
#     "build_system",
#     "__version__"
# ]

#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
PRISM - Protein Receptor Interaction Simulation Modeler

A comprehensive tool for building protein-ligand systems for molecular dynamics simulations.
"""

__version__ = "1.0.0"
__author__ = "PRISM Development Team"

from .builder import PRISMBuilder
from .core import PRISMSystem
from .sim import model

# High-level API functions
def system(protein_path, ligand_path, config=None, **kwargs):
    """
    Create a protein-ligand system for MD simulation.
    
    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str  
        Path to ligand file (MOL2/SDF)
    config : str or dict, optional
        Configuration file path or configuration dictionary
    **kwargs : optional
        Additional parameters (output_dir, ligand_forcefield, forcefield, etc.)
    
    Returns:
    --------
    PRISMSystem
        A system object with methods to build and run simulations
        
    Examples:
    ---------
    >>> import prism as pm
    >>> system = pm.system("protein.pdb", "ligand.mol2")
    >>> system.build()
    >>> system.generate_mdp_files()
    
    >>> # With custom configuration
    >>> system = pm.system("protein.pdb", "ligand.sdf", 
    ...                    config="config.yaml",
    ...                    ligand_forcefield="openff")
    >>> system.build()
    """
    return PRISMSystem(protein_path, ligand_path, config=config, **kwargs)

def build_system(protein_path, ligand_path, output_dir="prism_output", **kwargs):
    """
    Build a complete protein-ligand system (one-step function).
    
    Parameters:
    -----------
    protein_path : str
        Path to protein PDB file
    ligand_path : str
        Path to ligand file (MOL2/SDF)  
    output_dir : str
        Output directory for generated files
    **kwargs : optional
        Additional parameters
        
    Returns:
    --------
    str
        Path to the output directory
        
    Examples:
    ---------
    >>> import prism as pm
    >>> output_path = pm.build_system("protein.pdb", "ligand.mol2")
    """
    system_obj = PRISMSystem(protein_path, ligand_path, output_dir=output_dir, **kwargs)
    return system_obj.build()

# Export main classes and functions
__all__ = [
    "PRISMBuilder", 
    "PRISMSystem",
    "system", 
    "build_system",
    "model",
    "__version__"
]