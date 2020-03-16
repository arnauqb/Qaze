export G, C, M_P, K_B, SIGMA_T, SIGMA_SB, M_SUN, YEAR_TO_SEC, SIGMA_E, K_INTERP_K_VALUES, K_INTERP_XI_VALUES, ETAMAX_INTERP_ETAMAX_VALUES, ETAMAX_INTERP_XI_VALUES
const G = 6.6743e-8
const C = 29979245800.
const M_P = 1.67262192369e-24
const K_B = 1.380649e-16
const SIGMA_T = 6.6524587321000005e-25
const SIGMA_E = SIGMA_T / M_P
const SIGMA_SB = 5.6703744191844314e-05
const M_SUN = 1.988409870698051e+33
const YEAR_TO_SEC = 31557600.0

# interpolation values for force multiplier #
const K_INTERP_XI_VALUES = [-4, -3, -2.26, -2.00, -1.50, -1.00,
                            -0.42, 0.00, 0.22, 0.50, 1.0,
                            1.5, 1.8, 2.0, 2.18, 2.39,
                            2.76, 3.0, 3.29, 3.51, 3.68, 4.0]
const K_INTERP_K_VALUES = [0.411, 0.411, 0.400, 0.395, 0.363, 0.300,
                           0.200, 0.132, 0.100, 0.068, 0.042,
                           0.034, 0.033, 0.021, 0.013, 0.048,
                           0.046, 0.042, 0.044, 0.045, 0.032,
                            0.013]
const ETAMAX_INTERP_XI_VALUES = [-3, -2.5, -2.00, -1.50, -1.00,
                           -0.5, -0.23, 0.0, 0.32, 0.50,
                           1.0, 1.18, 1.50, 1.68, 2.0,
                           2.02, 2.16, 2.25, 2.39, 2.79,
                           3.0, 3.32, 3.50, 3.75, 4.00]

const ETAMAX_INTERP_ETAMAX_VALUES = [6.95, 6.95, 6.98, 7.05, 7.26,
                               7.56, 7.84, 8.00, 8.55, 8.95,
                               8.47, 8.00, 6.84, 6.00, 4.32,
                               4.00, 3.05, 2.74, 3.00, 3.10,
                               2.73, 2.00, 1.58, 1.20, 0.78]