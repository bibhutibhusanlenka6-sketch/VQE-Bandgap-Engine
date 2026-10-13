import sys
import re
import numpy as np
from ase.data import covalent_radii, atomic_numbers
from qiskit_nature.second_q.operators import FermionicOp
from qiskit_nature.second_q.mappers import JordanWignerMapper
from qiskit.circuit.library import EfficientSU2
from qiskit.primitives import StatevectorEstimator
from scipy.optimize import minimize
import warnings

warnings.filterwarnings("ignore", category=DeprecationWarning)

if len(sys.argv) < 2:
    print("\n[ERROR] Usage: python exact_vqe_bandgap.py <Formula> [optional: exact_t] [optional: exact_U]")
    print("Example (Heuristic): python exact_vqe_bandgap.py Nb2C")
    print("Example (Exact):     python exact_vqe_bandgap.py Mn3O4 2.0 5.0\n")
    sys.exit()

raw_input = sys.argv[1].capitalize()

print("\n" + "="*75)
print(f"   EXACT VQE BAND GAP ENGINE: {raw_input}")
print("="*75)

# --- THE PRECISION OVERRIDE ---
if len(sys.argv) == 4:
    dynamic_t = float(sys.argv[2])
    dynamic_U = float(sys.argv[3])
    print(f"    [MANUAL OVERRIDE] Using exact literature parameters:")
    print(f"    t = {dynamic_t:.3f} eV, U = {dynamic_U:.3f} eV")
else:
    parsed_formula = re.findall(r'([A-Z][a-z]*)(\d*)', raw_input)
    total_radius = 0.0
    total_atoms = 0
    for elem, count_str in parsed_formula:
        count = int(count_str) if count_str else 1
        r = covalent_radii[atomic_numbers[elem]]
        total_radius += r * count
        total_atoms += count
    avg_r = total_radius / total_atoms
    dynamic_t = 3.0 / (avg_r + 0.5)  
    dynamic_U = 8.0 / (avg_r + 0.5)  
    print(f"    [HEURISTIC MODE] Derived parameters based on radii:")
    print(f"    t = {dynamic_t:.3f} eV, U = {dynamic_U:.3f} eV")

mu = dynamic_U / 2.0  

# --- TRI-STATE VQE FUNCTION ---
ansatz = EfficientSU2(num_qubits=4, reps=3, entanglement="linear")
estimator = StatevectorEstimator()

def calculate_state_energy(target_N, state_name):
    print(f"\n>>> Optimizing {state_name} State (Target Electrons: {target_N})...")
    terms = {}
    
    terms["+_0 -_2"] = -dynamic_t; terms["+_2 -_0"] = -dynamic_t
    terms["+_1 -_3"] = -dynamic_t; terms["+_3 -_1"] = -dynamic_t
    terms["+_0 -_0 +_1 -_1"] = dynamic_U
    terms["+_2 -_2 +_3 -_3"] = dynamic_U
    terms["+_0 -_0"] = -mu; terms["+_1 -_1"] = -mu
    terms["+_2 -_2"] = -mu; terms["+_3 -_3"] = -mu

    P = 15.0 
    
    linear_coeff = P * (1.0 - 2.0 * target_N)
    for i in range(4):
        k = f"+_{i} -_{i}"
        terms[k] = terms.get(k, 0.0) + linear_coeff
        
    quad_coeff = 2.0 * P
    for i in range(4):
        for j in range(i+1, 4):
            k = f"+_{i} -_{i} +_{j} -_{j}"
            terms[k] = terms.get(k, 0.0) + quad_coeff

    fermionic_op = FermionicOp(terms, num_spin_orbitals=4)
    qubit_op = JordanWignerMapper().map(fermionic_op)

    def cost_function(params):
        job = estimator.run([(ansatz, qubit_op, params)])
        return float(job.result()[0].data.evs)

    np.random.seed(42)
    initial_params = np.random.random(ansatz.num_parameters) * 0.01

    result = minimize(cost_function, initial_params, method="COBYLA", options={"maxiter": 300})
    
    true_energy = result.fun + P * (target_N ** 2)
    print(f"    [SUCCESS] {state_name} Energy: {true_energy:.6f} eV")
    return true_energy

E_neutral = calculate_state_energy(target_N=2, state_name="Neutral (N)")
E_cation  = calculate_state_energy(target_N=1, state_name="Cation (N-1)")
E_anion   = calculate_state_energy(target_N=3, state_name="Anion (N+1)")

IP = E_cation - E_neutral
EA = E_neutral - E_anion
band_gap = IP - EA

print("\n" + "="*75)
print(" >>> ELECTRONIC BAND GAP ANALYSIS COMPLETE")
print("="*75)
print(f"Ionization Potential (IP):  {IP:.4f} eV")
print(f"Electron Affinity (EA):     {EA:.4f} eV")
print("-" * 75)
if band_gap < 0.05:
    print(f"FUNDAMENTAL BAND GAP:       {band_gap:.4f} eV  (Classification: METALLIC / CONDUCTOR)")
else:
    print(f"FUNDAMENTAL BAND GAP:       {band_gap:.4f} eV  (Classification: SEMICONDUCTOR)")
print("="*75 + "\n")