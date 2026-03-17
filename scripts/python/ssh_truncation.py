import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.special import genlaguerre, gammaln
import sys
import scienceplots
plt.style.use(['science', 'notebook'])


def compute_position_exponential(N, u): 
  matrix = np.zeros(shape=(N, N), dtype=complex)
  for m in range(N):
    for n in range(N):
      if m >= n:
        matrix[m][n] = (
          np.exp(-u**2/2)
          * np.exp(0.5 * (gammaln(n+1) - gammaln(m+1)))
          * (1.0j*u)**(m-n)
          * genlaguerre(n, m-n)(u**2)
        )
      else:
        matrix[m][n] = (
          np.exp(-u**2/2)
          * np.exp(0.5 * (gammaln(m+1) - gammaln(n+1)))
          * (1.0j*u)**(n-m)
          * genlaguerre(m, n-m)(u**2)
        )

  return qt.Qobj(matrix)
# hamiltoniano exacto con dimensión truncada en N estados de Fock
def exact_transmon_hamiltonian(N, E_c, E_j):
    charge_number = 0.5j * (E_j/(2*E_c))**0.25 * (qt.destroy(N).dag() - qt.destroy(N))
    phase_prefactor = (2 * E_c / E_j)**0.25
    cos_delta = 0.5 * (compute_position_exponential(N, phase_prefactor) + compute_position_exponential(N, -phase_prefactor))
    return 4 * E_c * charge_number**2 - E_j * cos_delta
# hamiltoniano para una cadena (anillo) de transmon
def exact_transmon_ssh(N, n, E_c, E_J, E_g1, E_g2):
    op_list = [qt.qeye(N) for _ in range(n)]
    full_hamiltonian = qt.qzero_like(qt.tensor(op_list))
    H_j = exact_transmon_hamiltonian(N, E_c, E_J)
    
    # suma de hamiltonianos actuando sobre cada uno de los transmon por separado
    for j in range(n):
        op_list[j] = H_j
        full_hamiltonian += qt.tensor(op_list)
        op_list[j] = qt.qeye(N)

    hbar_g1 = (E_J * E_c**3)**0.5 / (np.sqrt(2) * E_g1)
    hbar_g2 = (E_J * E_c**3)**0.5 / (np.sqrt(2) * E_g2)
    # sumamos términos de acoplamiento
    a = qt.destroy(N)
    charge_op = (a.dag() - a)

    for j in range(n):
        op_list[j] = charge_op
        op_list[(j+1) % n] = charge_op
        if j%2 == 0:
            full_hamiltonian -= hbar_g1 * qt.tensor(op_list)
        else:
            full_hamiltonian -= hbar_g2 * qt.tensor(op_list)

        op_list[j] = qt.qeye(N)
        op_list[(j+1) % n] = qt.qeye(N)
    
    return full_hamiltonian

# argumentos: numero total de pasos, indice

if len(sys.argv) != 2:
    print("usage: ./program_name total_steps index")


N_values = [5, 8]
n = 4
E_c = 1.0
E_j = 50.0
omega_q = np.sqrt(8 * E_c * E_j) - E_c
U = E_c / 2

nu_w_pairs = [(0.40, 0.60), (0.1, 0.9), (0.7, 0.3), (0.5, 0.2), (0.3, 0.4), (0.8, 0.1)]
E_g_values = 0.5 * (omega_q + 2 * U) * np.array(nu_w_pairs)

energy_matrix = []
for N in N_values:
    energies = []
    for E_g in E_g_values:
        hamiltonian = exact_transmon_ssh(N, n, E_c, E_j, *E_g)
        energies.append(np.sort(hamiltonian.eigenenergies()))
    energy_matrix.append(np.array(energies))