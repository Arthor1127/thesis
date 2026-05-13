import numpy as np
import matplotlib.pyplot as plt
import qutip as qt
from scipy.special import genlaguerre, gammaln
import itertools
import scienceplots
plt.style.use(['science', 'notebook'])

import scipy.sparse as sp
from scipy.sparse.linalg import eigsh


def ssh_energies(band, omega_q, nu, w, U, k):
    delta = np.sqrt(nu**2 + w**2 + 2 * nu * w * np.cos(k))
    if band == 0:
        energy = np.sqrt((omega_q)**2 - 4*U**2) * np.sqrt(1 + 2 * delta / (omega_q - 2 * U))
    elif band == 1:
        energy = np.sqrt((omega_q)**2 - 4*U**2) * np.sqrt(1 - 2 * delta / (omega_q - 2 * U))
        
    return energy

def ssh_ground_energy(n, omega_q, nu, w, U):
    k_values = np.linspace(-np.pi, np.pi, n, endpoint= False)
    energy_sum = np.sum([ssh_energies(0, omega_q, nu, w, U, k) + ssh_energies(1, omega_q, nu, w, U, k) for k in k_values])
    return 0.5 * energy_sum - n * omega_q

def embed(ops_at_indices: dict, n_total: int, fock_dim: int) -> qt.Qobj:
    """
    Embed a set of local operators into the full tensor product space.
    
    ops_at_indices: {index: qobj, ...}  — sites where non-identity ops go
    n_total: total number of sites (2*N for SSH, A and B interleaved)
    fock_dim: local Hilbert space dimension
    """
    tensor_aid = [qt.qeye(fock_dim)] * n_total
    for idx, op in ops_at_indices.items():
        tensor_aid[idx] = op
    return qt.tensor(tensor_aid)

def capacitive_ssh(N, fock_transmon, E_c, E_J, Ec1, Ec2):
    n_transmon = qt.num(fock_transmon)
    Lambda_0 = E_c
    capacitive_factor = np.sqrt(E_J * (E_c)**3)**(0.25)/np.sqrt(2)
    gc1 = capacitive_factor / Ec1
    gc2 = capacitive_factor / Ec2
    omega_q = np.sqrt(8*E_c*E_J) - E_c
    k_values = np.linspace(-np.pi, np.pi, N, endpoint= False)
    delta_ck = np.array([gc1 + np.exp(-1.0j * k) * gc2 for k in k_values])
    delta_ckA = Lambda_0 - np.abs(delta_ck)
    delta_ckB = Lambda_0 + np.abs(delta_ck)
    epsilon_ckA = omega_q + np.abs(delta_ck)
    epsilon_ckB = omega_q - np.abs(delta_ck)
    
    E_a = np.sqrt(epsilon_ckA**2 - delta_ckA**2)
    E_b = np.sqrt(epsilon_ckB**2 - delta_ckB**2)
    ground_state_energy = 0.5 * np.sum(E_a + E_b) - N * omega_q
    identity = qt.tensor([qt.qeye(fock_transmon)] * (2 * N))
    def make_n_kA():
        ops = []
        for j in range(N):
            ops.append(embed({2*j:   n_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    def make_n_kB():
        ops = []
        for j in range(N):
            ops.append(embed({2*j+1: n_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    return np.sum(E_a * make_n_kA() + E_b * make_n_kB()) + ground_state_energy * identity
    

def finally_transformed_squid_ssh(N, fock_transmon, E_c, E_J, Ec1, Ec2, EJ1, EJ2):
    omega_q = np.sqrt(8 * E_c * E_J) - E_c
    squid_factor = 2*np.sqrt(2 * E_c / E_J)
    gj1 = squid_factor * EJ1
    gj2 = squid_factor * EJ2
    k_values = np.linspace(-np.pi, np.pi, N, endpoint= False)
    capacitive_factor = np.sqrt(E_J * (E_c)**3)**(0.25)/np.sqrt(2)
    gc1 = capacitive_factor / Ec1
    gc2 = capacitive_factor / Ec2
    delta_ck = np.array([gc1 + np.exp(-1.0j * k) * gc2 for k in k_values])
    phase_correction = np.angle(delta_ck)
    delta_k = np.array([gj1 + np.exp(-1.0j * k) * gj2 for k in k_values]) * np.exp(-1.0j * phase_correction)
    
    diag_coeff_A_0 = (gj1 + gj2) - np.real(delta_k)
    diag_coeff_B_0 = (gj1 + gj2) + np.real(delta_k)
    self_pairing_coeff_A_0 = -0.5 * diag_coeff_B_0
    self_pairing_coeff_B_0 = -0.5 * diag_coeff_A_0
    hopping_coeff_0 = 1.0j * np.imag(delta_k)
    intersite_pairing_coeff_0 = hopping_coeff_0
    
    Lambda_0 = E_c
    
    delta_ckA = Lambda_0 - np.abs(delta_ck)
    delta_ckB = Lambda_0 + np.abs(delta_ck)
    epsilon_ckA = omega_q + np.abs(delta_ck)
    epsilon_ckB = omega_q - np.abs(delta_ck)
    
    E_a = np.sqrt(epsilon_ckA**2 - delta_ckA**2)
    E_b = np.sqrt(epsilon_ckB**2 - delta_ckB**2)
    
    cosh_2ra = epsilon_ckA / E_a
    cosh_2rb = epsilon_ckB / E_b
    sinh_2ra = delta_ckA / E_a
    sinh_2rb = delta_ckB / E_b
    
    ra = 0.5 * np.asinh(sinh_2ra)
    rb = 0.5 * np.asinh(sinh_2rb)
    
    n_transmon = qt.num(fock_transmon)
    b_transmon    = qt.destroy(fock_transmon)
    bdag_transmon = b_transmon.dag()
    
    
    def make_n_kA():
        ops = []
        for j in range(N):
            ops.append(embed({2*j:   n_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    def make_n_kB():
        ops = []
        for j in range(N):
            ops.append(embed({2*j+1: n_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    def make_self_pairing_kA():
        ops = []
        for j in range(N):
            mj = (N - j) % N
            
            if mj == j:
                ops.append(embed({2*j:   bdag_transmon * bdag_transmon}, 2*N, fock_transmon))
            else:
                ops.append(embed({2*j:   bdag_transmon, 2*mj:   bdag_transmon}, 2*N, fock_transmon))
                
        return np.array(ops)
    
    def make_self_pairing_kB():
        ops = []
        for j in range(N):
            mj = (N - j) % N
            
            if mj == j:
                ops.append(embed({2*j+1: bdag_transmon * bdag_transmon}, 2*N, fock_transmon))
            else:
                ops.append(embed({2*j+1: bdag_transmon, 2*mj+1: bdag_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    def make_hopping_k():
        ops  = []
        for j in range(N):
            ops.append(embed({2*j: bdag_transmon, 2*j+1: b_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    def make_intersite_pairing_k():
        ops  = []
        for j in range(N):
            mj = (N - j) % N
            ops.append(embed({2*j: bdag_transmon, 2*mj+1: bdag_transmon}, 2*N, fock_transmon))
        return np.array(ops)
    
    diag_coeff_A = diag_coeff_A_0 * cosh_2ra + 2.0 * self_pairing_coeff_A_0 * sinh_2ra
    diag_coeff_B = diag_coeff_B_0 * cosh_2rb + 2.0 * self_pairing_coeff_B_0 * sinh_2rb
    
    self_pairing_coeff_A =  self_pairing_coeff_A_0 * cosh_2ra + 0.5 * diag_coeff_A_0 * sinh_2ra
    self_pairing_coeff_B =  self_pairing_coeff_B_0 * cosh_2rb + 0.5 * diag_coeff_B_0 * sinh_2rb
    
    hopping_coeff = hopping_coeff_0 * np.exp(ra + rb)
    intersite_pairing_coeff = hopping_coeff
    
    ground_state_energy = 0.5 * np.sum(diag_coeff_A + diag_coeff_B - diag_coeff_A_0 - diag_coeff_B_0) + N * (0.5 * (gj1 + gj2) - 2.0 * (EJ1 + EJ2))
    H_diag = np.sum(diag_coeff_A * make_n_kA() + diag_coeff_B * make_n_kB())
    H_self_pairing = np.sum(self_pairing_coeff_A * make_self_pairing_kA() + self_pairing_coeff_B * make_self_pairing_kB())
    H_hop = np.sum(hopping_coeff * make_hopping_k())
    H_intersite_pairing = np.sum(intersite_pairing_coeff * make_intersite_pairing_k())
    identity = qt.tensor([qt.qeye(fock_transmon)] * (2 * N))
    
    return H_diag + H_self_pairing + H_self_pairing.dag() + H_hop + H_hop.dag() + H_intersite_pairing + H_intersite_pairing.dag() + ground_state_energy * identity

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

def ssh_squid_chain(N, fock_photon, fock_transmon, E_c, E_J, omega_c, gamma, Ec1, Ec2, E_J1, E_J2):
  # hamiltonianos de la cadena
  
  H0 = capacitive_ssh(N, fock_transmon, E_c, E_J, Ec1, Ec2)
  H1 = finally_transformed_squid_ssh(N, fock_transmon, E_c, E_J, Ec1, Ec2, E_J1, E_J2)
  eye_transmon = qt.tensor([qt.qeye(fock_transmon) for _ in range(2*N)])
  # hamiltonianos fotonicos
  eye_photon = qt.qeye(fock_photon)
  H_ph = omega_c * qt.num(fock_photon)
  sin_term = -0.5j * (compute_position_exponential(fock_photon, gamma) - compute_position_exponential(fock_photon, -gamma))
  
  return qt.tensor(H_ph, eye_transmon) + qt.tensor(eye_photon, H0) - qt.tensor(sin_term, H1)

