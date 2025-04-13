#https://www.ti.com/lit/ml/snaa106c/snaa106c.pdf PAGE 350
import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
from scipy import signal
import sympy as sp

def gain_brute_sanity( A0, T1, T2, T3, Icp, Fref, Fvco,  Kvco):
    # Calculate the transient response of the PLL system using the given parameters
    # This is a placeholder function. You can implement the actual calculation here.
    f = np.logspace(1, 10, 1000)  # Frequency range from 10 Hz to 10 GHz
    s = 1j * 2 * np.pi * f       # Complex frequency variable
    print(A0, T1, T2, T3, Icp, Fref, Fvco,  Kvco)
    Gain = A0
    Hol_numerator_filter =Gain* (s*T2+1)
    Hol_denominator_filter = s*A0*(1+s*T1)*(1+s*T3)
    Hol_tf_filter = Hol_numerator_filter/Hol_denominator_filter

    G_pll = Hol_tf_filter*Icp*Kvco/(s)
    N = Fvco/Fref
    H = 1.0/N
    Hol_num_pll = G_pll/(1+G_pll*H)
    # Plotting magnitude (gain)
    plt.figure(figsize=(8, 5))
    plt.semilogx(f, 20 * np.log10(np.abs(Hol_num_pll)), label='Magnitude (Gain)')
    plt.title("Frequency Response |H(jω)| (in dB)")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Gain (dB)")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()




def gain_brute_force( R2, R3, C1, C2, C3, Icp, Fref, Fvco,  Kvco):
    # Calculate the transient response of the PLL system using the given parameters
    # This is a placeholder function. You can implement the actual calculation here.
    f = np.logspace(1, 10, 1000)  # Frequency range from 10 Hz to 10 GHz
    s = 1j * 2 * np.pi * f       # Complex frequency variable
    den_s = ((C1+C2+C3)/(C1*C2*C3*R2*R3))
    dummy = C1*C2*R2+C1*C3*R3+C2*C3*R2+C2*R3*C3
    den_s2 = dummy/(C1*C2*C3*R2*R3)#s2
    den_s3 = 1
    den_s0 = 0
    Gain = 1.0/(C1*C3*R3)
    Hol_numerator_filter =Gain* (s+1.0/(R2*C2))
    Hol_denominator_filter = den_s3*s**3 + s**2 * (den_s2) + s * (den_s) + den_s0
    Hol_tf_filter = Hol_numerator_filter/Hol_denominator_filter

    G_pll = Hol_tf_filter*Icp*Kvco/(s)
    N = Fvco/Fref
    H = 1.0/N
    Hol_num_pll = G_pll/(1+G_pll*H)
    # Plotting magnitude (gain)
    plt.figure(figsize=(8, 5))
    plt.semilogx(f, 20 * np.log10(np.abs(Hol_num_pll)), label='Magnitude (Gain)')
    plt.title("Frequency Response |H(jω)| (in dB)")
    plt.xlabel("Frequency (Hz)")
    plt.ylabel("Gain (dB)")
    plt.grid(True, which='both', linestyle='--', linewidth=0.5)
    plt.tight_layout()
    plt.show()





def calculate_system_char( R2, R3, C1, C2, C3, Icp, Fref, Fvco,  Kvco):
    # Calculate the transient response of the PLL system using the given parameters
    # This is a placeholder function. You can implement the actual calculation here.

    s = sp.symbols('s')
    den_s = ((C1+C2+C3)/(C1*C2*C3*R2*R3))
    dummy = C1*C2*R2+C1*C3*R3+C2*C3*R2+C2*R3*C3
    den_s2 = dummy/(C1*C2*C3*R2*R3)#s2
    den_s3 = 1
    den_s0 = 0
    Gain = 1.0/(C1*C3*R3)
    Hol_numerator_filter =Gain* (s+1.0/(R2*C2))
    Hol_denominator_filter = den_s3*s**3 + s**2 * (den_s2) + s * (den_s) + den_s0
    Hol_tf_filter = Hol_numerator_filter/Hol_denominator_filter

    G_pll = Hol_tf_filter*Icp*Kvco/(s)
    N = Fvco/Fref
    H = 1.0/N


    H_cl = G_pll / (1 + G_pll * H)
    H_cl = sp.simplify(H_cl)

    # Get numerator and denominator as polynomials in 's'
    num, den = sp.fraction(H_cl)

    # Extract polynomial coefficients in descending powers of s
    num_coeffs = sp.Poly(num, s).all_coeffs()
    den_coeffs = sp.Poly(den, s).all_coeffs()

    # Convert symbolic coefficients to float for numerical analysis
    num_coeffs_f = [float(c.evalf()) for c in num_coeffs]
    den_coeffs_f = [float(c.evalf()) for c in den_coeffs]

    # Create LTI system
    sys = signal.TransferFunction(num_coeffs_f, den_coeffs_f)

    w, mag, phase = signal.bode(sys)
    # Find bandwidth frequency where mag drops below -3 dB from 0 Hz gain
    mag_linear = 10**(mag / 20)
    idx_bw = np.where(mag <= mag[0] - 3)[0]
    bw_hz = w[idx_bw[0]] / (2 * np.pi) if len(idx_bw) > 0 else None

    gcross_idx = np.where(np.diff(np.sign(mag)))[0]
    if len(gcross_idx) > 0:
        # Interpolate to get more precise frequency
        w_gc = np.interp(0, [mag[gcross_idx[0]], mag[gcross_idx[0] + 1]],
                            [w[gcross_idx[0]], w[gcross_idx[0] + 1]])
        phase_gc = np.interp(w_gc, w, phase)
        phase_margin = 180 + phase_gc
        print(f"Phase Margin: {phase_margin:.2f}° at {w_gc/(2*np.pi):.2e} Hz")
    else:
        print("No gain crossover found.")



    print(f"Estimated Bandwidth: {bw_hz:.2e} Hz")

    plt.figure(figsize=(10, 6))

    plt.subplot(2, 1, 1)
    plt.semilogx(w / (2*np.pi), mag)
    plt.title('Bode Plot')
    plt.ylabel('Magnitude (dB)')
    plt.grid(True, which='both')

    plt.subplot(2, 1, 2)
    plt.semilogx(w / (2*np.pi), phase)
    plt.xlabel('Frequency (Hz)')
    plt.ylabel('Phase (degrees)')
    plt.grid(True, which='both')

    plt.tight_layout()
    plt.show()






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
K_cp = 1e-3
Kvco = 30e6
fvco = 1.392e9
fref = 60e3
N = fvco/fref
K = K_cp * Kvco / N
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


calculate_system_char( R2, R3, C1, C2, C3, K_cp, fref, fvco,  Kvco)

#gain_brute_sanity( A0, T1, T2, T3, K_cp, fref, fvco,  Kvco)
#gain_brute_force( R2, R3, C1, C2, C3, K_cp, fref, fvco,  Kvco)
