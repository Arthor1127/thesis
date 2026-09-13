"""gds_stackup_pipeline - runs SQDMetal PALACE eigenmode simulations from a
raw GDS layout + an XML material stackup, via SQDMetal's native
gds_design= input branch. Requires SQDMetal to already be installed in the
active env (not declared as a pip dependency -- it's a local editable
install, see /home/ruiz/SQDmetal)."""
from .stackup import StackupInfo, resolve_stackup
from .build import run_gds_eigenmode

__all__ = ["StackupInfo", "resolve_stackup", "run_gds_eigenmode"]
