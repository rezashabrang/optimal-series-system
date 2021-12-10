from functions import (
    R,
    set_globals,
    W,
    calculate_equipment_cost
)
# ---------------------- Globals ----------------------
T = 210  # Time units of mission
J = 2  # N operation systems
H = 2  # N rescue systems
m = 60  # N intervals
pi = 50  # FDSM paramter
CF = 10000  # Failure Cost
CL = 300000  # System Loss Cost
CE = calculate_equipment_cost()  # Equipment Cost
set_globals(T, m, pi)

# ----------------- Calculating Probs -----------------
R_val = R()
W_val = W()
SS = R_val + W_val

# ----------------- Calculating total cost -----------------
C_tot = CF * (1 - R_val) + CL * (1 - SS) + CE

# ----------------- Printing results -----------------
print(f'R = {R_val}')
print(f'SS = {SS}')
print(f'C_E = {CE}')
print(f'C_tot = {C_tot}')
