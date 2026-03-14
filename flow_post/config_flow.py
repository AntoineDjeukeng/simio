from pathlib import Path

ROOT = Path("/data/antoine/Flow_CDI")
RESULTS = ROOT / "RESULTS"

# Set this once (constant across all simulations)
LY_NM = 5.112  # TODO: fill in your Ly in nm

FIELD_TO_E = {
    "FIELD_00": 0.0,
    "FIELD_01": 0.004,
    "FIELD_02": 0.01,
    "FIELD_03": 0.03,
}

FIELDS = ["FIELD_00", "FIELD_01", "FIELD_02", "FIELD_03"]
CHARGES = ["neutral", "pos", "neg", "positive", "negative"]
E_ORDER = [0.0, 0.004, 0.01, 0.03]

# Fit rule (methodology promise: long-time slope)
FIT_LAST_NS_DEFAULT = 45.0

# Linear-regime selection for conductivity and Ohmic fits.
# Set to None to disable threshold filtering.
SIGMA_E_MAX_V_NM = 0.01
OHMIC_E_MAX_V_NM = 0.01
OHMIC_E_LIST = None  # Example: [0.004, 0.01]

# Ohmic fit options
OHMIC_FIT_MODE = "WLS"   # "WLS" or "OLS"
OHMIC_WEIGHT = "I_sem"   # "I_sem" or "none"

# Replicas
I_START, I_END = 1, 20
PATTERN = "count_{i:02d}.dat"
