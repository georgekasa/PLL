#https://www.ti.com/lit/ml/snaa106c/snaa106c.pdf PAGE 350
import numpy as np
from scipy.optimize import fsolve

# Define the equation for T1
def equation(T1):
    term1 = np.arctan(gamma / (wc * T1 * (1 + T31)))
    term2 = np.arctan(wc * T1)
    term3 = np.arctan(wc * T1 * T31)
    return term1 - term2 - term3 - phi_rad



def sanity_check(T1, tolerance=1e-6):
    term1 = np.arctan(gamma / (wc * T1 * (1 + T31)))
    term2 = np.arctan(wc * T1)
    term3 = np.arctan(wc * T1 * T31)
    # Check if the equation is close to phi_rad
    return np.abs(term1 - term2 - term3 - phi_rad) < tolerance  # Adjust tolerance if needed

BW = 2e3
phase_margin = 47.1
gamma = 1.136
K_cp = 4e-3
Kvco = 30e6
fvco = 1.392e9
fref = 60e3
T31 = 0.6

N = fvco/fref
wc = 2*np.pi*BW

# Convert phi to radians, TAN & arctan are in radians!!!!
phi_rad = np.radians(phase_margin)


# Initial guess for T1
T1_initial_guess = 1e-8

# Solve for T1 using fsolve
T1_solution = fsolve(equation, T1_initial_guess)
T1 = T1_solution[0]
term1 = np.arctan(gamma / (wc * T1 * (1 + T31)))
term2 = np.arctan(wc * T1)
term3 = np.arctan(wc * T1 * T31)

print("sanity check: ", term1 - term2 - term3, phi_rad)
# Output the solution
print("Solution for T1:", T1_solution[0])
# Iterative solution with sanity check
T1_initial_guess = 1e-6 # Initial guess
max_iterations = 10    # Maximum number of iterations
tolerance = 1e-7        # Tolerance for the solution

for i in range(max_iterations):
    # Solve the equation
    T1_solution = fsolve(equation, T1_initial_guess)
    T1 = T1_solution[0]
    
    # Perform sanity check
    if sanity_check(T1):
        print(f"Solution found after {i+1} iterations: T1 = {T1}")
        break
    else:
        print(f"Iteration {i+1}: Solution {T1} did not pass sanity check. Trying a new guess.")
        T1_initial_guess *= 0.1  # Adjust guess (you can modify the step size)

else:
    print("Failed to find a valid solution within the maximum number of iterations.")




T3 = T1 * T31
T2 = gamma/((wc**2) * (T1 + T3))
A0_a = K_cp*Kvco/(N*wc**2)
A0_b_num = np.sqrt(1 + (wc*T2)**2)
A0_b_den = np.sqrt((1 + (wc*T1)**2)*(1 + (wc*T3)**2))
A0 = A0_a*A0_b_num/A0_b_den
A1 = A0*(T1+T3)
A2 = A0*T1*T3

#solve for components
C1 = (A2/(T2**2))*(1+np.sqrt(1+(T2/A2)*(T2*A0-A1)))
C3_num = -(T2*C1)**2  + T2*A1*C1 - A0*A2
C3 = C3_num/(C1*T2*T2 - A2)
C2 = A0 - C1 - C3
R2 = T2/C2
R3 = A2/(C1*C3*T2)

#print components
print("C1: ", C1)
print("C2: ", C2)
print("C3: ", C3)
print("R2: ", R2)
print("R3: ", R3) 


