import math
from SYSTEM import SPECS, PROPOSED_SYSTEM
from tqdm import tqdm

# -------------- Constants (globals) --------------

# Duration of simulation
T = 0

# Number of intervals
m = 1

# Delta
delta = 1

# FDSM Param
pi_fdsm = 99500


def set_globals(t_val, m_val, pi):
    global T
    global m
    global delta
    global pi_fdsm
    T = t_val
    m = m_val
    delta = T / m
    pi_fdsm = pi


def calculate_equipment_cost():
    """
    Calculating equipment Cost
    """
    total_cost = 0
    for key, val in PROPOSED_SYSTEM.items():
        spec = SPECS[key]
        for i in val:
            total_cost += spec[i]['c']
    return total_cost


def torque(t):
    """
    duration of successful rescue procedure when it is activated at time t
    args:
        t: Entire time from the entire mission begining
    """
    if t % 70 == 0:
        return 5
    t_c = t - 70 * math.floor(t / 70)
    res = 5 + 45 / (1 + math.pow(0.05 * t_c, -4.2))
    return res


def mu(v):
    """
    """
    global delta
    return torque(v * delta) / delta


def FDSM(t):
    """
    failure detection and switching mechanism probability
    args:
        t: Entire time from the entire mission begining
    """
    global pi_fdsm
    res = 1 - math.exp(-t / pi_fdsm)
    return res


def Z(i):
    """

    """
    global delta
    res_Z = FDSM(delta * i)
    if res_Z < 0:
        return 0
    return res_Z


def F_Weibull(t, j, k):
    """
    Weibull time-to failure distributions
    args:
        t: Entire time from the entire mission begining
        j: Subsystem type
        k: Component type
    """
    # Defining Constants
    n_j = SPECS[j][k]['n']
    b_j = SPECS[j][k]['b']
    res = 1 - math.exp((-1) * math.pow(t / n_j, b_j))
    return res


def p_hat(i_s, i_o, j, k):
    """
    probability that component of type component_type that should be activated
    in interval i_s fails before functioning for i_o time intervals
    args:
        i_s: Standby time interval,
        i_o: Opertaion time interval,
        j: System Type,
        k: Type of component
    """
    # Defining Constant
    d_j = SPECS[j][k]['d']
    return F_Weibull(delta * (d_j * i_s + i_o), j, k)


def p(i_s, i_o, j, k):
    """
    probability that component of type k that should be activated
    in interval i_s fails after functioning for exactly i_o time interval
    args:
        i_s: Standby time interval,
        i_o: Opertaion time interval,
        j: System Type,
        k: Type of component
    """
    return p_hat(i_s, i_o + 1, j, k) - p_hat(i_s, i_o, j, k)


def Q(j, k_j, x, k):
    """
    pmf of the random variable X(j,k) which represents the time interval that component k
    from subsystem j fails.
    args:
        j: Subsystem index
        k_j: Component index
        x: Time interval index
        k: Type of component
    """
    # First defining the inital values
    if k_j == 0:
        if x == 0:
            return 1
        else:
            return 0

    # Defining previous component type
    k = PROPOSED_SYSTEM[j][k_j - 1]

    # Calculating terms
    first_term = Q(j, k_j - 1, x, k) * p_hat(x, 0, j, k)
    second_term = 0
    for y in range(x):
        second_term += Q(j, k_j - 1, y, k) * p(y, x - y, j, k)

    return first_term + second_term


def U(j, i):
    """
    The probability that subsystem j fails in no later than time
    interval i.
    args:
        j: Subsystem index
        i: Time interval
        K_j: Last component of subsystem j
        x: Time interval
        k: Component Type
    """
    # Getting the last component specs of the subsystem j
    K_j = len(PROPOSED_SYSTEM[j])  # Getting the index (MAX should be four)
    k = PROPOSED_SYSTEM[j][-1]
    res = 0
    for x in range(i):
        res += Q(j, K_j, x, k)
    return res


def U_OS(i):
    """
    probability that the entire operation system fails in no later than time interval i
    args:
        i: Time interval index
        J: Operation subsystems
    """
    J = 2
    res = 1
    for j in range(1, J + 1):
        res *= (1 - U(j, i))
    return 1 - res


def Q_OS(i):
    """
    The pmf of the discrete time-to-failure of the entire operation system.
    args:
        i: Time interval index
    """
    return U_OS(i) - U_OS(i - 1)


# -------------------------------------------- MSP --------------------------------------------
def PR_e_j_k(j, k_j):
    """
    Probability of event e_j,k that component s_j(k) successfully accomplishes the mission task.
    """
    k = PROPOSED_SYSTEM[j][k_j - 1]
    res = 0
    for x in range(m):
        res += Q(j, k_j - 1, x, k) * (1 - p_hat(x, m - x, j, k))

    return res


def r(j):
    """
    The probability that subsystem j can successfully accomplish the mission.
    Args:
        j: Subsystem index
    """
    # Defining K_j
    K_j = len(PROPOSED_SYSTEM[j])

    res = 0
    for k_j in range(1, K_j + 1):
        res += PR_e_j_k(j, k_j)
    return res


def R():
    """
    MSP
    """
    J = 2
    res = 1
    for j in range(1, J + 1):
        res *= r(j)

    return res


# -------------------------------------------- SS --------------------------------------------

def Q_SS(j, k_j, v, x, k=0):
    """
    Q for SS
    args:
        j: Subsystem index
        k_j: Component index
        x: Time interval index
        k: Type of component
    """
    # First defining the inital values
    if k_j == 0:
        if x == v:
            return 1
        else:
            return 0
    else:
        if x < v:
            return 0
        elif v <= x:
            # Defining previous component type
            if k_j - 1 <= 0:
                curr_k = PROPOSED_SYSTEM[j][k_j - 1]
                prev_k = curr_k
            else:
                prev_k = PROPOSED_SYSTEM[j][k_j - 2]
                curr_k = PROPOSED_SYSTEM[j][k_j - 1]

            # Calculating terms
            # changed here
            first_term = Q_SS(j, k_j - 1, v, x, prev_k) * p_hat(v, 0, j, curr_k)
            second_term = 0
            for y in range(v, x + 1):
                second_term += Q_SS(j, k_j - 1, v, y, prev_k) * p(y, x - y, j, curr_k)

            return first_term + second_term
        else:
            return 1


def PR_e_j_k_SS(j, k_j, v):
    """
    SS PR_e_j_k
    Args:
        k_j: Index of component,
        v: Time interval
    """
    k = PROPOSED_SYSTEM[j][k_j - 1]
    res = 0
    for x in range(v, m):
        # changed here (keep in mind)
        i_o = v + mu(v) - x
        if i_o < 0:
            i_o = 0
        res += Q_SS(j, k_j - 1, v, x) * (1 - p_hat(x, i_o, j, k))

    return res


def w(j, v):
    """
    Probability that rescue system activated in time interval v
    completes the rescue procedure.
    Args:
        j: Subsystem index
        v: Time interval
    """
    subsystem = PROPOSED_SYSTEM[j]
    K_j = subsystem.index(subsystem[-1]) + 1
    res = 0
    for k_j in range(1, K_j + 1):
        res += PR_e_j_k_SS(j, k_j, v)
    return res


def W():
    """
    Overall success probability of the rescue procedure
    """
    global m
    res = 0
    # if no rescue system is provided then return 0 for this part
    if PROPOSED_SYSTEM[3] == [] and PROPOSED_SYSTEM[4] == []:
        return 0
    for v in tqdm(range(m)):
        w_res = 1
        for j in [3, 4]:
            # If there is no rescue system then skip it
            if PROPOSED_SYSTEM[j] == []:
                continue
            w_res *= w(j, v)
        res += Q_OS(v) * (1 - Z(v)) * w_res

    return res


def SS():
    """
    System survivability
    """
    return R() + W()
