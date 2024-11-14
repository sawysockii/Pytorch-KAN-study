"""
Модуль с константами, используемыми в расчетных модулях
"""
from scipy.interpolate import interp1d

# Таблица корректировочных коэффициентов по книге Такаса NU_SSU_POINTS = [50, 80, 100, 150, 200,
# 300, 400, 500, 600, 700, 800, 900, 1000, 1500, 2000, 2500, 3000, 4000, 5000] Q_COEFF_POINTS = [
# 1, 0.980, 0.970, 0.947, 0.924, 0.886, 0.847, 0.819, 0.792, 0.766, 0.745, 0.727, 0.708, 0.659,
# 0.621, 0.590, 0.562, 0.518, 0.479] HEAD_COEFFS_POINTS = [1, 0.990, 0.985, 0.970, 0.958, 0.933,
# 0.909, 0.897, 0.883, 0.868, 0.858, 0.846, 0.833, 0.799, 0.771, 0.750, 0.733, 0.702,
# 0.677] EFF_COEFF_POINTS = [0.945, 0.870, 0.825, 0.736, 0.674, 0.566, 0.497, 0.462, 0.434,
# 0.410, 0.390, 0.368, 0.349, 0.307, 0.272, 0.245, 0.218, 0.178, 0.149] POWER_COEFF_POINTS = [
# 1.058, 1.115, 1.158, 1.248, 1.341, 1.460, 1.549, 1.590, 1.611, 1.622, 1.639, 1.671, 1.690,
# 1.715, 1.760, 1.806, 1.890, 2.043, 2.176]

# Таблица корректировочных коэффициентов вязкости по Pipesim
NU_SSU_POINTS = [
    46.96,
    51.56,
    63.33,
    73.82,
    82.81,
    114.60,
    166.93,
    254.13,
    278.56,
    404.03,
    628.13,
    702.33,
    938.15,
    1472.40,
    2500.58,
    4925.41,
]
POWER_COEFF_POINTS = [
    1.053,
    1.053,
    1.089,
    1.121,
    1.144,
    1.212,
    1.292,
    1.384,
    1.404,
    1.487,
    1.586,
    1.611,
    1.677,
    1.779,
    1.898,
    2.046,
]
Q_COEFF_POINTS = [
    1.000,
    1.000,
    0.994,
    0.988,
    0.983,
    0.964,
    0.935,
    0.894,
    0.884,
    0.841,
    0.784,
    0.769,
    0.728,
    0.663,
    0.584,
    0.485,
]
EFF_COEFF_POINTS = [
    0.950,
    0.950,
    0.910,
    0.876,
    0.852,
    0.780,
    0.698,
    0.609,
    0.590,
    0.516,
    0.434,
    0.415,
    0.367,
    0.299,
    0.231,
    0.161,
]
HEAD_COEFFS_POINTS = [
    1.000,
    1.000,
    0.997,
    0.994,
    0.991,
    0.980,
    0.964,
    0.942,
    0.936,
    0.911,
    0.878,
    0.869,
    0.844,
    0.803,
    0.751,
    0.680,
]

# Список переменных для вывода распределений в скважинах
DISTRS_PVT = [
    "rs",
    "pb",
    "muo",
    "mug",
    "muw",
    "mu_liq",
    "mu_mix",
    "z",
    "bo",
    "bg",
    "bw",
    "rho_oil_rc",
    "rho_gas_rc",
    "rho_wat_rc",
    "rho_liq_rc",
    "rho_mix_rc",
    "compro",
    "q_oil_rc",
    "q_gas_rc",
    "q_wat_rc",
    "q_liq_rc",
    "q_mix_rc",
    "gas_fraction",
    "st_oil_gas",
    "st_wat_gas",
    "st_liq_gas",
]
DISTRS_PIPE = [
    "dp_dl",
    "dp_dl_fric",
    "dp_dl_grav",
    "dp_dl_acc",
    "liquid_holdup",
    "friction_factor",
    "vsl",
    "vsg",
    "vsm",
    "flow_pattern",
    "lambda_l",
    "n_re",
    "angle",
    "vl",
    "vg",
    # Минимальная скорость выноса жидкости с забоя
    "v_mix",
    "v_mix_krit",
]

DISTRS = DISTRS_PVT + DISTRS_PIPE

DISTRS_NONE = [
    "rs",
    "pb",
    "muo",
    "muw",
    "bo",
    "bw",
    "compro",
    "flow_pattern",
    "v_mix",
    "v_mix_krit"
    ]

# Список PVT-корреляций по умолчанию
OIL_CORRS = {
    "pb": "standing",
    "rs": "standing",
    "b": "standing",
    "mud": "beggs",
    "mus": "beggs",
    "compr": "vasquez",
    "hc": "wright",
    "st_oil_gas": "baker",
}
WAT_CORRS = {
    "b": "mccain",
    "mu": "mccain",
    "rho": "standing",
    "hc": "const",
    "st_wat_gas": "katz",
}
GAS_CORRS = {
    "ppc": "standing",
    "tpc": "standing",
    "z": "standing",
    "mu": "lee",
    "hc": "mahmood",
}

# Словарь свойств и функций для BlackOilModel
ALL_PROPERTIES = {
    "bw": ["wat_corrs", "calc_water_fvf"],
    "rho_wat": ["wat_corrs", "calc_water_density"],
    "muw": ["wat_corrs", "calc_water_viscosity"],
    "hc_wat": ["wat_corrs", "calc_heat_capacity"],
    "salinity": ["wat_corrs", "calc_salinity"],
    "z": ["gas_corrs", "calc_z"],
    "bg": ["gas_corrs", "calc_gas_fvf"],
    "rho_gas": ["gas_corrs", "calc_gas_density"],
    "mug": ["gas_corrs", "calc_gas_viscosity"],
    "hc_gas": ["gas_corrs", "calc_heat_capacity"],
    "pb": ["oil_corrs", "calc_pb"],
    "rs": ["oil_corrs", "calc_rs"],
    "compro": ["oil_corrs", "calc_oil_compressibility"],
    "bo": ["oil_corrs", "calc_oil_fvf"],
    "rho_oil": ["oil_corrs", "calc_oil_density"],
    "muo": ["oil_corrs", "calc_oil_viscosity"],
    "hc_oil": ["oil_corrs", "calc_heat_capacity"],
    "st_wat_gas": ["wat_corrs", "calc_st_wat_gas"],
    "st_oil_gas": ["oil_corrs", "calc_st_oil_gas"],
}

# Максимальное значение давления насыщения = 20 000 psi
PB_MAX = 137895145.863367

# 1 Атмосфера = 101325 Па
ATM = 101325

# Стандартные корреляции для расчета штуцера
CHOKE_CORRS = {"subcritical": "mechanistic", "critical": "mechanistic"}

# Показатель адиабаты
CP_CV = 1.26

# Фактор критического давления (2 / (CP_CV + 1)) ** (CP_CV / (CP_CV - 1))
CPR = 0.5530618435484577

# Коэффициент перевода psi в Па
PSI = 0.00014503773773020924

# Верхняя граница оптимизатора для давления = 700 атм
P_UP_LIMIT = 70927500

# Наборы констант для расчета корреляции Hagedor&Brown
PHI_HAGEDORN_BROWN_CONSTANTS = [1, 1.11, 1.4, 1.6, 1.69, 1.749, 1.785, 1.81, 1.84]
CB_HAGEDORN_BROWN_CONSTANTS = [0.012, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.089]

# Безразмерное число n_l в корреляции Hagedor&Brown
N_L_CONSTANTS = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
N_LC_CONSTANTS = [0.002, 0.0022, 0.0025, 0.0029, 0.0041, 0.0066, 0.009, 0.0145]

C_A_CONSTANTS = [
    0.000002,
    0.000005,
    0.00001,
    0.000015,
    0.00002,
    0.00005,
    0.0001,
    0.0002,
    0.0005,
    0.001,
    0.002,
    0.005,
    0.01,
]
HL_PHI_CONSTANTS = [
    0.04,
    0.09,
    0.145,
    0.17,
    0.18,
    0.25,
    0.34,
    0.44,
    0.65,
    0.82,
    0.92,
    0.96,
    1,
]

# Наборы констант для расчета корреляции Hagedor&Brown
PHI_HAGEDORN_BROWN_CONSTANTS = [1, 1.11, 1.4, 1.6, 1.69, 1.749, 1.785, 1.81, 1.84]
CB_HAGEDORN_BROWN_CONSTANTS = [0.012, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.089]

# Безразмерное число n_l в корреляции Hagedor&Brown
N_L_CONSTANTS = [0.002, 0.005, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5]
N_LC_CONSTANTS = [0.002, 0.0022, 0.0025, 0.0029, 0.0041, 0.0066, 0.009, 0.0145]

C_A_CONSTANTS = [
    0.000002,
    0.000005,
    0.00001,
    0.000015,
    0.00002,
    0.00005,
    0.0001,
    0.0002,
    0.0005,
    0.001,
    0.002,
    0.005,
    0.01,
]

HL_PHI_CONSTANTS = [
    0.04,
    0.09,
    0.145,
    0.17,
    0.18,
    0.25,
    0.34,
    0.44,
    0.65,
    0.82,
    0.92,
    0.96,
    1,
]

# Массив корректировочный коэффициентов влияния количества свободного газа на приеме насоса по SPE 206468
GAS_COEF_CORR = [
    [
        1,
        0.94488189,
        0.535433071,
        0.322834646,
        0.299212598,
        0.291338583,
        0.267716535,
        0.251968504,
        0.220472441,
        0.196850394,
        0.173228346,
        0.141732283,
    ],
    [
        1,
        0.952380952,
        0.634920635,
        0.428571429,
        0.396825397,
        0.380952381,
        0.357142857,
        0.333333333,
        0.293650794,
        0.261904762,
        0.23015873,
        0.19047619,
    ],
    [1, 0.96, 0.736, 0.536, 0.496, 0.472, 0.448, 0.416, 0.368, 0.328, 0.288, 0.24],
    [
        1,
        0.97199341,
        0.827018122,
        0.658978583,
        0.601317957,
        0.57660626,
        0.551894563,
        0.510708402,
        0.461285008,
        0.39538715,
        0.345963756,
        0.28830313,
    ],
    [
        1,
        0.979831933,
        0.87394958,
        0.722689076,
        0.655462185,
        0.621848739,
        0.596638655,
        0.554621849,
        0.487394958,
        0.428571429,
        0.369747899,
        0.305882353,
    ],
    [
        1,
        0.982905983,
        0.907692308,
        0.774358974,
        0.700854701,
        0.666666667,
        0.632478632,
        0.581196581,
        0.514529915,
        0.444444444,
        0.384615385,
        0.290598291,
    ],
    [
        1,
        0.98540146,
        0.948905109,
        0.857664234,
        0.775547445,
        0.708029197,
        0.638686131,
        0.565693431,
        0.501824818,
        0.419708029,
        0.328467153,
        0.200729927,
    ],
    [
        1,
        0.983935743,
        0.963855422,
        0.893574297,
        0.799196787,
        0.698795181,
        0.592369478,
        0.512048193,
        0.441767068,
        0.335341365,
        0.200803213,
        -1.11468 * 10 ** (-17),
    ],
    [
        1,
        0.985324948,
        0.964360587,
        0.903563941,
        0.807127883,
        0.693920335,
        0.58490566,
        0.503144654,
        0.419287212,
        0.299790356,
        0.146750524,
        -0.1,
    ],
    [
        1,
        0.983796296,
        0.965277778,
        0.909722222,
        0.805555556,
        0.659722222,
        0.546296296,
        0.428240741,
        0.300925926,
        0.122685185,
        -0.162037037,
        -0.6,
    ],
    [
        1,
        0.987951807,
        0.963855422,
        0.910843373,
        0.8,
        0.648192771,
        0.530120482,
        0.397590361,
        0.240963855,
        0.024096386,
        -0.361445783,
        -1.1,
    ],
    [
        1,
        0.986842105,
        0.960526316,
        0.907894737,
        0.784210526,
        0.626315789,
        0.497368421,
        0.342105263,
        0.131578947,
        -0.263157895,
        -0.789473684,
        -1.789473684,
    ],
    [
        1,
        0.984615385,
        0.950769231,
        0.892307692,
        0.753846154,
        0.584615385,
        0.415384615,
        0.215384615,
        -0.307692308,
        -1,
        -2,
        -3,
    ],
    [
        1,
        0.98245614,
        0.936842105,
        0.870175439,
        0.698245614,
        0.526315789,
        0.333333333,
        -0.070175439,
        -1,
        -2,
        -3,
        -4,
    ],
    [
        1,
        0.980392157,
        0.921568627,
        0.843137255,
        0.654901961,
        0.470588235,
        0.254901961,
        -0.392156863,
        -2,
        -3,
        -4,
        -5,
    ],
    [
        1,
        0.980392157,
        0.892156863,
        0.794117647,
        0.56372549,
        0.367647059,
        0,
        -1.5,
        -5,
        -6,
        -7,
        -8,
    ],
    [
        1,
        0.981818182,
        0.848484848,
        0.727272727,
        0.478787879,
        0.181818182,
        -0.363636364,
        -3,
        -7,
        -8,
        -9,
        -10,
    ],
    [
        1,
        0.954545455,
        0.718181818,
        0.536363636,
        0.227272727,
        -0.454545455,
        -1.818181818,
        -6,
        -10,
        -11,
        -12,
        -13,
    ],
    [1, 0.9, 0.2, 0.05, -0.5, -2, -5, -11, -15, -16, -17, -18],
    [1, 0.7, -0.5, -0.7, -1.5, -5, -10, -17, -20, -21, -22, -24],
]
# Диапазон газосодержания, д.ед
GAS_COEF_CORR_X = [
    0,
    0.01,
    0.02,
    0.03,
    0.05,
    0.07,
    0.09,
    0.12,
    0.15,
    0.18,
    0.21,
    0.25,
]
# Диапазон коэффицинта подачи, д.ед
GAS_COEF_CORR_Y = [
    0,
    0.11,
    0.21,
    0.32,
    0.37,
    0.42,
    0.53,
    0.63,
    0.65,
    0.72,
    0.74,
    0.77,
    0.81,
    0.84,
    0.86,
    0.89,
    0.92,
    0.95,
    1.0,
    1.05,
]

# Наборы констант для расчета корреляции Duns & Ros

# 1) Для расчета границ между пузырьковым и пробковым режимами
L1_x = [0, 16, 20, 30, 40, 50, 70, 275]
L1_y = [2, 2, 1.9, 1.6, 1.25, 1.1, 1, 1]

L2_x = [7.5, 20, 30, 40, 50, 60, 100, 275]
L2_y = [0.5, 0.75, 0.9, 1, 1.1, 1.1, 1.11, 1.11]


# 2) Для расчета hl
N_l_x = [0.002, 0.006, 0.007, 0.01, 0.02, 0.04, 0.05, 0.1, 0.2, 0.4, 2]
F1_y = [1.3, 1.3, 1.3, 1.3, 1.3, 1.5, 1.65, 2, 2, 1.8, 0.9]
F2_y = [0.25, 0.25, 0.25, 0.25, 0.28, 0.45, 0.6, 0.95, 1, 0.98, 0.7]
F3_y = [0.8, 0.9, 1, 1.3, 1.9, 2.5, 3, 3.25, 3.5, 3.75, 4]
F4_y = [-20, 5, 10, 25, 38, 50, 52, 55, 55, 55, 55]
F5_y = [0.22, 0.2, 0.19, 0.18, 0.17, 0.15, 0.13, 0.065, 0.048, 0.06, 0.11]
F6_y = [0.8, 0.05, 0, -0.1, -0.1, 0.65, 1, 2.1, 1.9, 1.8, 1.75]
F7_y = [0.14, 0.101, 0.099, 0.09, 0.07, 0.056, 0.052, 0.04, 0.034, 0.03, 0.025]

# 3) Для расчета ff_l (коэффициент трения)
coef_x = [0.001, 0.4, 0.7, 1, 2, 3, 6, 10, 20, 40, 100]
f2_y = [1, 1, 0.9, 0.75, 0.6, 0.5, 0.4, 0.35, 0.3, 0.25, 0.2]

f1_func = interp1d(
    x=N_l_x,
    y=F1_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f2_func = interp1d(
    x=N_l_x,
    y=F2_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f3_func = interp1d(
    x=N_l_x,
    y=F3_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f4_func = interp1d(
    x=N_l_x,
    y=F4_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f5_func = interp1d(
    x=N_l_x,
    y=F5_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f6_func = interp1d(
    x=N_l_x,
    y=F6_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f7_func = interp1d(
    x=N_l_x,
    y=F7_y,
    fill_value="extrapolate",
    kind="quadratic",
)
f2_y_func = interp1d(
    x=coef_x,
    y=f2_y,
    fill_value="extrapolate",
    kind="linear",
)
l1_func = interp1d(
    x=L1_x,
    y=L1_y,
    fill_value="extrapolate",
    kind="quadratic",
)
l2_func = interp1d(
    x=L2_x,
    y=L2_y,
    fill_value="extrapolate",
    kind="quadratic",
)

DR_const = {
    "f1": f1_func,
    "f2": f2_func,
    "f3": f3_func,
    "f4": f4_func,
    "f5": f5_func,
    "f6": f6_func,
    "f7": f7_func,
    "f2_y": f2_y_func,
    "l1": l1_func,
    "l2": l2_func,
}

# Константы, используемые в расчете температурной корреляции
TEMP_CORR = {
    "time_sec": 1209600,  # Время, которое отработала скважина на момент расчета
    "cp_earth": 837.4,  # Теплоемкость грунта
    "cp_an": 1004.81,  # Теплоемкость флюида в затрубном пространстве
    "thermal_conduct_tube": 60.55,  # Теплопроводность металла НКТ
    "thermal_conduct_cem": 0.779,  # Теплопроводность цемента вокруг скважины
    "thermal_conduct_earth": 2.422,  # Теплопроводность грунта
    "thermal_conduct_an": 0.865,  # Теплопроводность в затрубном пространстве
    "thermal_conduct_l": 0.1384,  # Теплопроводность жидкости
    "thermal_conduct_g": 0.0346,  # Теплопроводность газа
    "rho_earth": 2242,  # Плотность породы
    "mu_an": 0.0001,  # Вязкость флюида в затрубном пространстве
    "rho_an": 36.92,  # Плотность флюида в затрубном пространстве
    "u": 28.39132,  # Коэффициент теплопередачи системы скважинный флюид/скважина/горная порода
}