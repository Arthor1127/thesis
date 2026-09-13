"""sq_engine -- compatibility shim.

The SecondQuantizationEngine and the SSH PHOTONICChain now live together in the
single module ``ssh_chain_sq`` (installed in the qcg-quantum-design env). This
shim re-exports the engine half so ``from sq_engine import ...`` keeps working.

Note: inside ``notebooks/exploration/unitary/`` the original two local files
(``sq_engine.py`` / ``ssh_chain_sq.py``) still shadow this, so nothing there
changes.
"""
from ssh_chain_sq import (  # noqa: F401
    TensorTerm,
    Superposition,
    QuantumOperator,
    SecondQuantizationEngine,
    CompositeEngine,
    _CUPY_AVAILABLE,
)
