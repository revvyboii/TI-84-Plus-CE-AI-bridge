/*
 * CHEM - Chemistry 1 problem solver for the TI-84 Plus CE (offline)
 *
 * Same engine and screens as PHYS: pick a topic, type the numbers the word
 * problem gives you, and every other quantity fills in by itself.
 * ENTER on a solved value shows the steps.
 *
 * 23 problem types: density, % error, temperature, moles, % by mass,
 * average atomic mass, empirical formula, molarity, dilution, ideal gas,
 * combined gas law, STP, partial pressures, stoichiometry, limiting
 * reactant, % yield, q = mc dT, calorimetry, heat of reaction, light,
 * hydrogen atom, pH, titration.
 *
 * Keys (form): digits . (-) type a value   ^ or , = E (6.02E23)
 *   LEFT = negative, RIGHT = positive   UP/DOWN field   DEL backspace
 *   CLEAR erase / go back   ENTER steps for a solved value   MODE theme
 *   Y= menu  WINDOW reset  ZOOM tips  TRACE steps  GRAPH periodic table
 * Periodic table (GRAPH in a topic, WINDOW at home): arrows move,
 *   ALPHA + letter jumps (ALPHA C = carbon, again = Cl, Ca...),
 *   ENTER adds one atom, digits set how many, DEL removes,
 *   TRACE = % composition, GRAPH = put the molar mass in the box.
 * Home or topic list: ZOOM = Help me pick (asks a few questions, then
 *   walks you box by box: ENTER = next box, UP = back a box).
 */

#ifndef USE_CLAUDE
#define USE_CLAUDE 0     /* offline only */
#endif

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

/* ================================================================== */
/* Solver engine (plain C, also compiled on a PC for testing)          */
/* ================================================================== */

#define MAXV    16
#define GRAV    9.80
#define PI      3.14159265358979
#define DEG     (PI / 180.0)
#define BIGG    6.674e-11
#define UNK     (0.0 / 0.0)     /* "no default" marker */
#define B(x)    (1u << (x))

static bool ok(double x) { return x == x && fabs(x) < 1e37; }

typedef struct {
    const char *sym;    /* short symbol shown in the form   */
    const char *unit;
    const char *desc;   /* help line + sent to Claude       */
    double def;         /* value assumed if left blank, or UNK */
} var_t;

typedef double (*fn_t)(const double *v);
/* rule flags: either sign of a square root / skip the agreement check / has a second root */
enum { R_ABS = 1, R_NOCHK = 2, R_ALT = 4 };
typedef struct {
    uint8_t target;
    uint16_t need;
    fn_t fn;
    const char *eq;
    uint8_t flags;
} rule_t;

typedef struct {
    const char *name;
    const var_t *vars;
    uint8_t nvars;
    const rule_t *rules;
    uint8_t nrules;
    const char *tips;   /* word-problem phrases -> what to type */
} topic_t;

static uint16_t flip;   /* quantities to take with the other root (see solve) */
static double clamp1(double x) { return x > 1 && x < 1.00001 ? 1 : x < -1 && x > -1.00001 ? -1 : x; }
static double acosz(double x) { return acos(clamp1(x)); }
static double asinz(double x) { return asin(clamp1(x)); }

#define FN(name, expr) static double name(const double *v) { return (expr); }
#define SQ(x) ((x) * (x))
static double root(double x) { return x < 0 ? UNK : sqrt(x); }
/* square root of a sum whose terms add up to `scale` in size: a slightly
   negative result is rounding in what was typed (e.g. "ends at rest"), so 0 */
static double rootz(double x, double scale) {
    if (!ok(x) || !ok(scale)) return UNK;
    if (x < 0 && -x <= 2e-3 * scale) return 0;
    return root(x);
}

/* ================================================================== */
/* Chemistry 1 topics                                                  */
/* ================================================================== */
#define NA      6.022e23        /* Avogadro, per mol */
#define RGAS    0.08206         /* L atm / (mol K) */
#define HPL     6.626e-34       /* Planck, J s */
#define CLIGHT  2.998e8         /* m/s */
#define RYD     2.18e-18        /* hydrogen energy constant, J */
#define VSTP    22.4            /* L per mol at STP */
#define ZEROC   273.15
#define CWATER  4.184           /* J / (g C) */
static double min2(double a, double b) { return a < b ? a : b; }
/* a few shortcut rules apply only while some boxes are still blank */
static uint8_t st[MAXV];        /* box states (the form state below; 0 = blank) */
#define BLANK(i) (st[i] == 0)
static double cbrt_(double x) { return x > 0 ? pow(x, 1.0 / 3) : x < 0 ? -pow(-x, 1.0 / 3) : 0; }

/* ---------- density ------------------------------------------------- */
enum { DN_M, DN_V, DN_D, DN_VI, DN_VF, DN_L };
FN(dn_d,  v[DN_M] / v[DN_V])
FN(dn_m,  v[DN_D] * v[DN_V])
FN(dn_v,  v[DN_M] / v[DN_D])
FN(dn_vw, v[DN_VF] - v[DN_VI])
FN(dn_vc, v[DN_L] * v[DN_L] * v[DN_L])
FN(dn_vf, v[DN_VI] + v[DN_V])
FN(dn_vi, v[DN_VF] - v[DN_V])
FN(dn_l,  v[DN_V] > 0 ? cbrt_(v[DN_V]) : UNK)
static const var_t dens_vars[] = {
    {"m",  "g",    "m: mass", UNK},
    {"V",  "mL",   "V: volume (1 mL = 1 cm^3)", UNK},
    {"d",  "g/mL", "d: density", UNK},
    {"Vi", "mL",   "Vi: water level before", UNK},
    {"Vf", "mL",   "Vf: water level with object in", UNK},
    {"L",  "cm",   "L: side of a cube", UNK},
};
#define DN_EQ "d = m / V"
#define DW_EQ "V = Vf - Vi (water pushed up)"
static const rule_t dens_rules[] = {
    {DN_D,  B(DN_M)|B(DN_V), dn_d, DN_EQ, 0},
    {DN_V,  B(DN_VI)|B(DN_VF), dn_vw, DW_EQ, 0},
    {DN_V,  B(DN_L), dn_vc, "V = L^3 (cube)", 0},
    {DN_M,  B(DN_D)|B(DN_V), dn_m, DN_EQ, 0},
    {DN_V,  B(DN_M)|B(DN_D), dn_v, DN_EQ, 0},
    {DN_VF, B(DN_VI)|B(DN_V), dn_vf, DW_EQ, 0},
    {DN_VI, B(DN_VF)|B(DN_V), dn_vi, DW_EQ, 0},
    {DN_L,  B(DN_V), dn_l, "V = L^3 (cube)", 0},
};

/* ---------- percent error ------------------------------------------- */
enum { PE_X, PE_A, PE_E };
FN(pe_e, fabs(v[PE_X] - v[PE_A]) / fabs(v[PE_A]) * 100)
static const var_t perr_vars[] = {
    {"meas", "", "meas: what you measured", UNK},
    {"true", "", "true: accepted (book) value", UNK},
    {"err",  "%", "err: percent error", UNK},
};
static const rule_t perr_rules[] = {
    {PE_E, B(PE_X)|B(PE_A), pe_e, "err = |meas - true| / true x 100", 0},
};

/* ---------- temperature units --------------------------------------- */
enum { TM_C, TM_K, TM_F };
FN(tm_k,  v[TM_C] + ZEROC)
FN(tm_ck, v[TM_K] - ZEROC)
FN(tm_f,  1.8 * v[TM_C] + 32)
FN(tm_cf, (v[TM_F] - 32) / 1.8)
static const var_t temp_vars[] = {
    {"C", "deg C", "C: Celsius", UNK},
    {"K", "K",     "K: kelvin (never negative)", UNK},
    {"F", "deg F", "F: Fahrenheit", UNK},
};
static const rule_t temp_rules[] = {
    {TM_C, B(TM_K), tm_ck, "K = C + 273.15", 0},
    {TM_C, B(TM_F), tm_cf, "F = 1.8 C + 32", 0},
    {TM_K, B(TM_C), tm_k,  "K = C + 273.15", 0},
    {TM_F, B(TM_C), tm_f,  "F = 1.8 C + 32", 0},
};

/* ---------- grams, moles, particles --------------------------------- */
enum { MO_M, MO_MM, MO_N, MO_P, MO_K, MO_AT };
FN(mo_n,  v[MO_M] / v[MO_MM])
FN(mo_np, v[MO_P] / NA)
FN(mo_m,  v[MO_N] * v[MO_MM])
FN(mo_mm, v[MO_M] / v[MO_N])
FN(mo_p,  v[MO_N] * NA)
FN(mo_pa, v[MO_AT] / v[MO_K])
FN(mo_at, v[MO_P] * v[MO_K])
static const var_t mole_vars[] = {
    {"m",   "g",     "m: mass in grams", UNK},
    {"M",   "g/mol", "M: molar mass (GRAPH = table)", UNK},
    {"n",   "mol",   "n: moles", UNK},
    {"N",   "",      "N: molecules (or formula units)", UNK},
    {"k",   "",      "k: atoms of X in one (H2O: H=2)", UNK},
    {"Nat", "",      "Nat: atoms of that one kind", UNK},
};
static const rule_t mole_rules[] = {
    {MO_N,  B(MO_M)|B(MO_MM), mo_n, "n = m / M", 0},
    {MO_N,  B(MO_P), mo_np, "N = n x 6.022E23", 0},
    {MO_P,  B(MO_AT)|B(MO_K), mo_pa, "Nat = N x k", 0},
    {MO_M,  B(MO_N)|B(MO_MM), mo_m, "n = m / M", 0},
    {MO_MM, B(MO_M)|B(MO_N), mo_mm, "n = m / M", 0},
    {MO_P,  B(MO_N), mo_p, "N = n x 6.022E23", 0},
    {MO_AT, B(MO_P)|B(MO_K), mo_at, "Nat = N x k", 0},
};

/* ---------- percent by mass ----------------------------------------- */
enum { PM_PART, PM_WHOLE, PM_PCT };
FN(pm_pct,   100 * v[PM_PART] / v[PM_WHOLE])
FN(pm_part,  v[PM_PCT] * v[PM_WHOLE] / 100)
FN(pm_whole, 100 * v[PM_PART] / v[PM_PCT])
static const var_t pmass_vars[] = {
    {"part",  "g", "part: grams of the element/part", UNK},
    {"whole", "g", "whole: grams of the whole thing", UNK},
    {"pct",   "%", "pct: percent by mass", UNK},
};
#define PM_EQ "pct = part / whole x 100"
static const rule_t pmass_rules[] = {
    {PM_PCT,   B(PM_PART)|B(PM_WHOLE), pm_pct, PM_EQ, 0},
    {PM_PART,  B(PM_PCT)|B(PM_WHOLE), pm_part, PM_EQ, 0},
    {PM_WHOLE, B(PM_PART)|B(PM_PCT), pm_whole, PM_EQ, 0},
};

/* ---------- average atomic mass ------------------------------------- */
enum { AV_M1, AV_F1, AV_M2, AV_F2, AV_M3, AV_F3, AV_AVG };
FN(av_avg, (v[AV_M1] * v[AV_F1] + v[AV_M2] * v[AV_F2] + v[AV_M3] * v[AV_F3]) / 100)
FN(av_f2,  100 - v[AV_F1] - v[AV_F3])
FN(av_f1r, 100 - v[AV_F2] - v[AV_F3])
FN(av_f3,  100 - v[AV_F1] - v[AV_F2])
FN(av_f1,  (100 * (v[AV_AVG] - v[AV_M2]) + v[AV_F3] * (v[AV_M2] - v[AV_M3])) / (v[AV_M1] - v[AV_M2]))
FN(av_m2,  (100 * v[AV_AVG] - v[AV_M1] * v[AV_F1] - v[AV_M3] * v[AV_F3]) / v[AV_F2])
static const var_t avg_vars[] = {
    {"m1",  "amu", "m1: mass of isotope 1", UNK},
    {"f1",  "%",   "f1: its abundance (75.78%)", UNK},
    {"m2",  "amu", "m2: mass of isotope 2", UNK},
    {"f2",  "%",   "f2: its abundance (blank = rest)", UNK},
    {"m3",  "amu", "m3: isotope 3 (blank = none)", 0},
    {"f3",  "%",   "f3: its abundance", 0},
    {"avg", "amu", "avg: average atomic mass", UNK},
};
#define AV_EQ "avg = (m1 f1 + m2 f2 + m3 f3) / 100"
static const rule_t avg_rules[] = {
    {AV_F2,  B(AV_F1)|B(AV_F3), av_f2, "f1 + f2 + f3 = 100", 0},
    {AV_F1,  B(AV_F2)|B(AV_F3), av_f1r, "f1 + f2 + f3 = 100", 0},
    {AV_F3,  B(AV_F1)|B(AV_F2), av_f3, "f1 + f2 + f3 = 100", 0},
    {AV_AVG, B(AV_M1)|B(AV_F1)|B(AV_M2)|B(AV_F2)|B(AV_M3)|B(AV_F3), av_avg, AV_EQ, 0},
    {AV_F1,  B(AV_AVG)|B(AV_M1)|B(AV_M2)|B(AV_M3)|B(AV_F3), av_f1, AV_EQ, 0},
    {AV_M2,  B(AV_AVG)|B(AV_M1)|B(AV_F1)|B(AV_F2)|B(AV_M3)|B(AV_F3), av_m2, AV_EQ, 0},
};

/* ---------- empirical formula --------------------------------------- */
enum { EF_G1, EF_A1, EF_G2, EF_A2, EF_G3, EF_A3, EF_N1, EF_N2, EF_N3,
       EF_R1, EF_R2, EF_R3, EF_ME, EF_MM, EF_X };
FN(ef_n1, v[EF_G1] / v[EF_A1])
FN(ef_n2, v[EF_G2] / v[EF_A2])
FN(ef_n3, v[EF_G3] / v[EF_A3])
#define EF_MIN3 min2(min2(v[EF_N1], v[EF_N2]), v[EF_N3])
#define EF_MIN2 min2(v[EF_N1], v[EF_N2])
FN(ef_r1c, v[EF_N1] / EF_MIN3)
FN(ef_r2c, v[EF_N2] / EF_MIN3)
FN(ef_r3c, v[EF_N3] / EF_MIN3)
FN(ef_r1b, BLANK(EF_G3) ? v[EF_N1] / EF_MIN2 : UNK)
FN(ef_r2b, BLANK(EF_G3) ? v[EF_N2] / EF_MIN2 : UNK)
FN(ef_x,  v[EF_MM] / v[EF_ME])
static const var_t emp_vars[] = {
    {"g1",   "g",     "g1: grams (or %) of element 1", UNK},
    {"A1",   "g/mol", "A1: its atomic mass (GRAPH)", UNK},
    {"g2",   "g",     "g2: grams (or %) of element 2", UNK},
    {"A2",   "g/mol", "A2: its atomic mass (GRAPH)", UNK},
    {"g3",   "g",     "g3: element 3 (blank = none)", UNK},
    {"A3",   "g/mol", "A3: its atomic mass", UNK},
    {"n1",   "mol",   "n1: moles of element 1", UNK},
    {"n2",   "mol",   "n2: moles of element 2", UNK},
    {"n3",   "mol",   "n3: moles of element 3", UNK},
    {"r1",   "",      "r1: subscript 1 (round it)", UNK},
    {"r2",   "",      "r2: subscript 2 (round it)", UNK},
    {"r3",   "",      "r3: subscript 3 (round it)", UNK},
    {"Memp", "g/mol", "Memp: empirical formula mass", UNK},
    {"Mmol", "g/mol", "Mmol: molecular mass (if given)", UNK},
    {"x",    "",      "x: multiply subscripts by this", UNK},
};
#define EF_REQ "r = n / (smallest n)"
static const rule_t emp_rules[] = {
    {EF_N1, B(EF_G1)|B(EF_A1), ef_n1, "n = g / A", 0},
    {EF_N2, B(EF_G2)|B(EF_A2), ef_n2, "n = g / A", 0},
    {EF_N3, B(EF_G3)|B(EF_A3), ef_n3, "n = g / A", 0},
    {EF_R1, B(EF_N1)|B(EF_N2)|B(EF_N3), ef_r1c, EF_REQ, 0},
    {EF_R2, B(EF_N1)|B(EF_N2)|B(EF_N3), ef_r2c, EF_REQ, 0},
    {EF_R3, B(EF_N1)|B(EF_N2)|B(EF_N3), ef_r3c, EF_REQ, 0},
    {EF_R1, B(EF_N1)|B(EF_N2), ef_r1b, EF_REQ, R_NOCHK},     /* only 2 elements */
    {EF_R2, B(EF_N1)|B(EF_N2), ef_r2b, EF_REQ, R_NOCHK},
    {EF_X,  B(EF_MM)|B(EF_ME), ef_x, "x = Mmol / Memp", 0},
};

/* ---------- molarity ------------------------------------------------ */
enum { ML_M, ML_MM, ML_N, ML_ML, ML_V, ML_C };
FN(ml_n,  v[ML_M] / v[ML_MM])
FN(ml_m,  v[ML_N] * v[ML_MM])
FN(ml_mm, v[ML_M] / v[ML_N])
FN(ml_v,  v[ML_ML] / 1000)
FN(ml_ml, v[ML_V] * 1000)
FN(ml_c,  v[ML_N] / v[ML_V])
FN(ml_nc, v[ML_C] * v[ML_V])
FN(ml_vc, v[ML_N] / v[ML_C])
static const var_t mlr_vars[] = {
    {"m",  "g",     "m: grams of solute", UNK},
    {"M",  "g/mol", "M: molar mass (GRAPH = table)", UNK},
    {"n",  "mol",   "n: moles of solute", UNK},
    {"mL", "mL",    "mL: volume of SOLUTION", UNK},
    {"V",  "L",     "V: same volume in liters", UNK},
    {"c",  "mol/L", "c: molarity ('0.50 M')", UNK},
};
#define ML_EQ "c = n / V"
static const rule_t mlr_rules[] = {
    {ML_N,  B(ML_M)|B(ML_MM), ml_n, "n = m / M", 0},
    {ML_V,  B(ML_ML), ml_v, "V = mL / 1000", 0},
    {ML_ML, B(ML_V), ml_ml, "V = mL / 1000", 0},
    {ML_C,  B(ML_N)|B(ML_V), ml_c, ML_EQ, 0},
    {ML_N,  B(ML_C)|B(ML_V), ml_nc, ML_EQ, 0},
    {ML_V,  B(ML_N)|B(ML_C), ml_vc, ML_EQ, 0},
    {ML_M,  B(ML_N)|B(ML_MM), ml_m, "n = m / M", 0},
    {ML_MM, B(ML_M)|B(ML_N), ml_mm, "n = m / M", 0},
};

/* ---------- dilution ------------------------------------------------ */
enum { DL_C1, DL_V1, DL_C2, DL_V2, DL_W };
FN(dl_c1, v[DL_C2] * v[DL_V2] / v[DL_V1])
FN(dl_v1, v[DL_C2] * v[DL_V2] / v[DL_C1])
FN(dl_c2, v[DL_C1] * v[DL_V1] / v[DL_V2])
FN(dl_v2, v[DL_C1] * v[DL_V1] / v[DL_C2])
FN(dl_w,  v[DL_V2] - v[DL_V1])
FN(dl_v2w, v[DL_V1] + v[DL_W])
FN(dl_v1w, v[DL_V2] - v[DL_W])
static const var_t dil_vars[] = {
    {"C1",    "M",  "C1: stock (strong) molarity", UNK},
    {"V1",    "mL", "V1: volume of stock used", UNK},
    {"C2",    "M",  "C2: molarity after diluting", UNK},
    {"V2",    "mL", "V2: final TOTAL volume", UNK},
    {"water", "mL", "water: water added", UNK},
};
#define DL_EQ "C1 V1 = C2 V2"
#define DLW_EQ "V2 = V1 + water"
static const rule_t dil_rules[] = {
    {DL_V2, B(DL_V1)|B(DL_W), dl_v2w, DLW_EQ, 0},
    {DL_V1, B(DL_V2)|B(DL_W), dl_v1w, DLW_EQ, 0},
    {DL_C1, B(DL_C2)|B(DL_V2)|B(DL_V1), dl_c1, DL_EQ, 0},
    {DL_V1, B(DL_C2)|B(DL_V2)|B(DL_C1), dl_v1, DL_EQ, 0},
    {DL_C2, B(DL_C1)|B(DL_V1)|B(DL_V2), dl_c2, DL_EQ, 0},
    {DL_V2, B(DL_C1)|B(DL_V1)|B(DL_C2), dl_v2, DL_EQ, 0},
    {DL_W,  B(DL_V1)|B(DL_V2), dl_w, DLW_EQ, 0},
};

/* ---------- ideal gas ----------------------------------------------- */
enum { IG_P, IG_V, IG_N, IG_T, IG_TC, IG_M, IG_MM, IG_D, IG_KPA, IG_MMHG };
FN(ig_t,   v[IG_TC] + ZEROC)
FN(ig_tc,  v[IG_T] - ZEROC)
FN(ig_pk,  v[IG_KPA] / 101.325)
FN(ig_pm,  v[IG_MMHG] / 760)
FN(ig_kpa, v[IG_P] * 101.325)
FN(ig_mmhg, v[IG_P] * 760)
FN(ig_n,   v[IG_M] / v[IG_MM])
FN(ig_p,   v[IG_N] * RGAS * v[IG_T] / v[IG_V])
FN(ig_v,   v[IG_N] * RGAS * v[IG_T] / v[IG_P])
FN(ig_np,  v[IG_P] * v[IG_V] / (RGAS * v[IG_T]))
FN(ig_tt,  v[IG_P] * v[IG_V] / (v[IG_N] * RGAS))
FN(ig_m,   v[IG_N] * v[IG_MM])
FN(ig_mm,  v[IG_M] / v[IG_N])
FN(ig_d,   v[IG_M] / v[IG_V])
FN(ig_dp,  v[IG_P] * v[IG_MM] / (RGAS * v[IG_T]))
FN(ig_mmd, v[IG_D] * RGAS * v[IG_T] / v[IG_P])
static const var_t gas_vars[] = {
    {"P",    "atm",   "P: pressure in atm", UNK},
    {"V",    "L",     "V: volume in liters", UNK},
    {"n",    "mol",   "n: moles of gas", UNK},
    {"T",    "K",     "T: temperature in kelvin", UNK},
    {"Tc",   "deg C", "Tc: or type Celsius here", UNK},
    {"m",    "g",     "m: grams of gas", UNK},
    {"M",    "g/mol", "M: molar mass (GRAPH = table)", UNK},
    {"d",    "g/L",   "d: density of the gas", UNK},
    {"kPa",  "kPa",   "kPa: or pressure in kPa", UNK},
    {"mmHg", "mmHg",  "mmHg: or in mmHg (= torr)", UNK},
};
#define IG_EQ "P V = n R T   (R = 0.08206)"
static const rule_t gas_rules[] = {
    {IG_T,    B(IG_TC), ig_t, "T = Tc + 273.15", 0},
    {IG_P,    B(IG_KPA), ig_pk, "1 atm = 101.325 kPa", 0},
    {IG_P,    B(IG_MMHG), ig_pm, "1 atm = 760 mmHg", 0},
    {IG_N,    B(IG_M)|B(IG_MM), ig_n, "n = m / M", 0},
    {IG_P,    B(IG_N)|B(IG_T)|B(IG_V), ig_p, IG_EQ, 0},
    {IG_V,    B(IG_N)|B(IG_T)|B(IG_P), ig_v, IG_EQ, 0},
    {IG_N,    B(IG_P)|B(IG_V)|B(IG_T), ig_np, IG_EQ, 0},
    {IG_T,    B(IG_P)|B(IG_V)|B(IG_N), ig_tt, IG_EQ, 0},
    {IG_MM,   B(IG_D)|B(IG_T)|B(IG_P), ig_mmd, "d = P M / (R T)", 0},
    {IG_M,    B(IG_N)|B(IG_MM), ig_m, "n = m / M", 0},
    {IG_MM,   B(IG_M)|B(IG_N), ig_mm, "n = m / M", 0},
    {IG_D,    B(IG_M)|B(IG_V), ig_d, "d = m / V", 0},
    {IG_D,    B(IG_P)|B(IG_MM)|B(IG_T), ig_dp, "d = P M / (R T)", 0},
    {IG_TC,   B(IG_T), ig_tc, "T = Tc + 273.15", 0},
    {IG_KPA,  B(IG_P), ig_kpa, "1 atm = 101.325 kPa", 0},
    {IG_MMHG, B(IG_P), ig_mmhg, "1 atm = 760 mmHg", 0},
};

/* ---------- combined gas law ---------------------------------------- */
/* Before and after. A pair left blank is held constant (Boyle, Charles,
   Gay-Lussac): those shortcut rules skip the agreement check because they
   only apply when the third pair really is constant. */
enum { CG_P1, CG_V1, CG_T1, CG_P2, CG_V2, CG_T2, CG_C1, CG_C2 };
FN(cg_t1c, v[CG_C1] + ZEROC)
FN(cg_t2c, v[CG_C2] + ZEROC)
FN(cg_c1,  v[CG_T1] - ZEROC)
FN(cg_c2,  v[CG_T2] - ZEROC)
FN(cg_p2,  v[CG_P1] * v[CG_V1] * v[CG_T2] / (v[CG_T1] * v[CG_V2]))
FN(cg_v2,  v[CG_P1] * v[CG_V1] * v[CG_T2] / (v[CG_T1] * v[CG_P2]))
FN(cg_t2,  v[CG_P2] * v[CG_V2] * v[CG_T1] / (v[CG_P1] * v[CG_V1]))
FN(cg_p1,  v[CG_P2] * v[CG_V2] * v[CG_T1] / (v[CG_V1] * v[CG_T2]))
FN(cg_v1,  v[CG_P2] * v[CG_V2] * v[CG_T1] / (v[CG_P1] * v[CG_T2]))
FN(cg_t1,  v[CG_P1] * v[CG_V1] * v[CG_T2] / (v[CG_P2] * v[CG_V2]))
#define SAME_T (BLANK(CG_T1) && BLANK(CG_T2))
#define SAME_P (BLANK(CG_P1) && BLANK(CG_P2))
#define SAME_V (BLANK(CG_V1) && BLANK(CG_V2))
FN(by_p2,  SAME_T ? v[CG_P1] * v[CG_V1] / v[CG_V2] : UNK)
FN(by_v2,  SAME_T ? v[CG_P1] * v[CG_V1] / v[CG_P2] : UNK)
FN(by_p1,  SAME_T ? v[CG_P2] * v[CG_V2] / v[CG_V1] : UNK)
FN(by_v1,  SAME_T ? v[CG_P2] * v[CG_V2] / v[CG_P1] : UNK)
FN(ch_v2,  SAME_P ? v[CG_V1] * v[CG_T2] / v[CG_T1] : UNK)
FN(ch_t2,  SAME_P ? v[CG_T1] * v[CG_V2] / v[CG_V1] : UNK)
FN(ch_v1,  SAME_P ? v[CG_V2] * v[CG_T1] / v[CG_T2] : UNK)
FN(ch_t1,  SAME_P ? v[CG_T2] * v[CG_V1] / v[CG_V2] : UNK)
FN(gl_p2,  SAME_V ? v[CG_P1] * v[CG_T2] / v[CG_T1] : UNK)
FN(gl_t2,  SAME_V ? v[CG_T1] * v[CG_P2] / v[CG_P1] : UNK)
FN(gl_p1,  SAME_V ? v[CG_P2] * v[CG_T1] / v[CG_T2] : UNK)
FN(gl_t1,  SAME_V ? v[CG_T2] * v[CG_P1] / v[CG_P2] : UNK)
static const var_t comb_vars[] = {
    {"P1",  "(any)", "P1: pressure before", UNK},
    {"V1",  "(any)", "V1: volume before", UNK},
    {"T1",  "K",     "T1: kelvin before (same T: blank)", UNK},
    {"P2",  "(any)", "P2: pressure after", UNK},
    {"V2",  "(any)", "V2: volume after", UNK},
    {"T2",  "K",     "T2: kelvin after", UNK},
    {"T1c", "deg C", "T1c: or Celsius before", UNK},
    {"T2c", "deg C", "T2c: or Celsius after", UNK},
};
#define CG_EQ "P1 V1 / T1 = P2 V2 / T2"
#define BY_EQ "same T (Boyle): P1 V1 = P2 V2"
#define CH_EQ "same P (Charles): V1/T1 = V2/T2"
#define GL_EQ "same V: P1 / T1 = P2 / T2"
static const rule_t comb_rules[] = {
    {CG_T1, B(CG_C1), cg_t1c, "T = Tc + 273.15", 0},
    {CG_T2, B(CG_C2), cg_t2c, "T = Tc + 273.15", 0},
    {CG_P2, B(CG_P1)|B(CG_V1)|B(CG_T1)|B(CG_V2)|B(CG_T2), cg_p2, CG_EQ, 0},
    {CG_V2, B(CG_P1)|B(CG_V1)|B(CG_T1)|B(CG_P2)|B(CG_T2), cg_v2, CG_EQ, 0},
    {CG_T2, B(CG_P1)|B(CG_V1)|B(CG_T1)|B(CG_P2)|B(CG_V2), cg_t2, CG_EQ, 0},
    {CG_P1, B(CG_V1)|B(CG_T1)|B(CG_P2)|B(CG_V2)|B(CG_T2), cg_p1, CG_EQ, 0},
    {CG_V1, B(CG_P1)|B(CG_T1)|B(CG_P2)|B(CG_V2)|B(CG_T2), cg_v1, CG_EQ, 0},
    {CG_T1, B(CG_P1)|B(CG_V1)|B(CG_P2)|B(CG_V2)|B(CG_T2), cg_t1, CG_EQ, 0},
    {CG_P2, B(CG_P1)|B(CG_V1)|B(CG_V2), by_p2, BY_EQ, R_NOCHK},
    {CG_V2, B(CG_P1)|B(CG_V1)|B(CG_P2), by_v2, BY_EQ, R_NOCHK},
    {CG_P1, B(CG_P2)|B(CG_V2)|B(CG_V1), by_p1, BY_EQ, R_NOCHK},
    {CG_V1, B(CG_P2)|B(CG_V2)|B(CG_P1), by_v1, BY_EQ, R_NOCHK},
    {CG_V2, B(CG_V1)|B(CG_T1)|B(CG_T2), ch_v2, CH_EQ, R_NOCHK},
    {CG_T2, B(CG_T1)|B(CG_V1)|B(CG_V2), ch_t2, CH_EQ, R_NOCHK},
    {CG_V1, B(CG_V2)|B(CG_T1)|B(CG_T2), ch_v1, CH_EQ, R_NOCHK},
    {CG_T1, B(CG_T2)|B(CG_V1)|B(CG_V2), ch_t1, CH_EQ, R_NOCHK},
    {CG_P2, B(CG_P1)|B(CG_T1)|B(CG_T2), gl_p2, GL_EQ, R_NOCHK},
    {CG_T2, B(CG_T1)|B(CG_P1)|B(CG_P2), gl_t2, GL_EQ, R_NOCHK},
    {CG_P1, B(CG_P2)|B(CG_T1)|B(CG_T2), gl_p1, GL_EQ, R_NOCHK},
    {CG_T1, B(CG_T2)|B(CG_P1)|B(CG_P2), gl_t1, GL_EQ, R_NOCHK},
    {CG_C1, B(CG_T1), cg_c1, "T = Tc + 273.15", 0},
    {CG_C2, B(CG_T2), cg_c2, "T = Tc + 273.15", 0},
};

/* ---------- gas at STP ---------------------------------------------- */
enum { SP_N, SP_V, SP_M, SP_MM, SP_P };
FN(sp_v,  v[SP_N] * VSTP)
FN(sp_nv, v[SP_V] / VSTP)
FN(sp_n,  v[SP_M] / v[SP_MM])
FN(sp_np, v[SP_P] / NA)
FN(sp_m,  v[SP_N] * v[SP_MM])
FN(sp_mm, v[SP_M] / v[SP_N])
FN(sp_p,  v[SP_N] * NA)
static const var_t stp_vars[] = {
    {"n", "mol",   "n: moles of gas", UNK},
    {"V", "L",     "V: liters at STP (0 C, 1 atm)", UNK},
    {"m", "g",     "m: grams of gas", UNK},
    {"M", "g/mol", "M: molar mass (GRAPH = table)", UNK},
    {"N", "",      "N: molecules", UNK},
};
static const rule_t stp_rules[] = {
    {SP_N,  B(SP_V), sp_nv, "V = n x 22.4 L/mol", 0},
    {SP_N,  B(SP_M)|B(SP_MM), sp_n, "n = m / M", 0},
    {SP_N,  B(SP_P), sp_np, "N = n x 6.022E23", 0},
    {SP_V,  B(SP_N), sp_v, "V = n x 22.4 L/mol", 0},
    {SP_M,  B(SP_N)|B(SP_MM), sp_m, "n = m / M", 0},
    {SP_MM, B(SP_M)|B(SP_N), sp_mm, "n = m / M", 0},
    {SP_P,  B(SP_N), sp_p, "N = n x 6.022E23", 0},
};

/* ---------- partial pressures --------------------------------------- */
enum { DP_PT, DP_P1, DP_P2, DP_P3, DP_N1, DP_NT, DP_X1 };
FN(dp_pt, v[DP_P1] + v[DP_P2] + v[DP_P3])
FN(dp_p1, v[DP_PT] - v[DP_P2] - v[DP_P3])
FN(dp_p2, v[DP_PT] - v[DP_P1] - v[DP_P3])
FN(dp_x,  v[DP_N1] / v[DP_NT])
FN(dp_n1, v[DP_X1] * v[DP_NT])
FN(dp_nt, v[DP_N1] / v[DP_X1])
FN(dp_px, v[DP_X1] * v[DP_PT])
FN(dp_xp, v[DP_P1] / v[DP_PT])
FN(dp_tx, v[DP_P1] / v[DP_X1])
static const var_t dal_vars[] = {
    {"Ptot", "(any)", "Ptot: total pressure", UNK},
    {"P1",   "(any)", "P1: gas 1 (the one asked about)", UNK},
    {"P2",   "(any)", "P2: gas 2 (or water vapor)", UNK},
    {"P3",   "(any)", "P3: gas 3 (blank = none)", 0},
    {"n1",   "mol",   "n1: moles of gas 1", UNK},
    {"ntot", "mol",   "ntot: total moles of gas", UNK},
    {"x1",   "",      "x1: mole fraction of gas 1", UNK},
};
#define DP_EQ "Ptot = P1 + P2 + P3"
static const rule_t dal_rules[] = {
    {DP_X1, B(DP_N1)|B(DP_NT), dp_x, "x1 = n1 / ntot", 0},
    {DP_PT, B(DP_P1)|B(DP_P2)|B(DP_P3), dp_pt, DP_EQ, 0},
    {DP_P1, B(DP_PT)|B(DP_P2)|B(DP_P3), dp_p1, DP_EQ, 0},
    {DP_P1, B(DP_X1)|B(DP_PT), dp_px, "P1 = x1 Ptot", 0},
    {DP_P2, B(DP_PT)|B(DP_P1)|B(DP_P3), dp_p2, DP_EQ, 0},
    {DP_X1, B(DP_P1)|B(DP_PT), dp_xp, "P1 = x1 Ptot", 0},
    {DP_PT, B(DP_P1)|B(DP_X1), dp_tx, "P1 = x1 Ptot", 0},
    {DP_N1, B(DP_X1)|B(DP_NT), dp_n1, "x1 = n1 / ntot", 0},
    {DP_NT, B(DP_N1)|B(DP_X1), dp_nt, "x1 = n1 / ntot", 0},
};

/* ---------- stoichiometry: grams of A -> grams of B ------------------ */
enum { SO_MA, SO_MMA, SO_NA, SO_A, SO_B, SO_NB, SO_MMB, SO_MB, SO_VB };
FN(so_na,  v[SO_MA] / v[SO_MMA])
FN(so_ma,  v[SO_NA] * v[SO_MMA])
FN(so_mma, v[SO_MA] / v[SO_NA])
FN(so_nb,  v[SO_NA] * v[SO_B] / v[SO_A])
FN(so_nab, v[SO_NB] * v[SO_A] / v[SO_B])
FN(so_mb,  v[SO_NB] * v[SO_MMB])
FN(so_nbm, v[SO_MB] / v[SO_MMB])
FN(so_mmb, v[SO_MB] / v[SO_NB])
FN(so_vb,  v[SO_NB] * VSTP)
FN(so_nbv, v[SO_VB] / VSTP)
static const var_t sto_vars[] = {
    {"mA", "g",     "mA: grams of what you HAVE", UNK},
    {"MA", "g/mol", "MA: its molar mass (GRAPH)", UNK},
    {"nA", "mol",   "nA: moles of what you have", UNK},
    {"a",  "",      "a: its coefficient (blank = 1)", 1},
    {"b",  "",      "b: WANTED one's coefficient", 1},
    {"nB", "mol",   "nB: moles of what you want", UNK},
    {"MB", "g/mol", "MB: its molar mass (GRAPH)", UNK},
    {"mB", "g",     "mB: grams of what you WANT", UNK},
    {"VB", "L",     "VB: liters, if a gas at STP", UNK},
};
#define SO_EQ "nB = nA x b / a (mole ratio)"
static const rule_t sto_rules[] = {
    {SO_NA,  B(SO_MA)|B(SO_MMA), so_na, "n = m / M", 0},
    {SO_NB,  B(SO_MB)|B(SO_MMB), so_nbm, "n = m / M", 0},
    {SO_NB,  B(SO_VB), so_nbv, "V = n x 22.4 L/mol (STP)", 0},
    {SO_NB,  B(SO_NA)|B(SO_A)|B(SO_B), so_nb, SO_EQ, 0},
    {SO_NA,  B(SO_NB)|B(SO_A)|B(SO_B), so_nab, SO_EQ, 0},
    {SO_MB,  B(SO_NB)|B(SO_MMB), so_mb, "n = m / M", 0},
    {SO_MA,  B(SO_NA)|B(SO_MMA), so_ma, "n = m / M", 0},
    {SO_VB,  B(SO_NB), so_vb, "V = n x 22.4 L/mol (STP)", 0},
    {SO_MMA, B(SO_MA)|B(SO_NA), so_mma, "n = m / M", 0},
    {SO_MMB, B(SO_MB)|B(SO_NB), so_mmb, "n = m / M", 0},
};

/* ---------- limiting reactant --------------------------------------- */
enum { LR_MA, LR_MMA, LR_A, LR_MB, LR_MMB, LR_B, LR_P, LR_MMP,
       LR_NA, LR_NB, LR_NP, LR_MP, LR_LIM, LR_LEFT };
#define LR_RA (v[LR_NA] / v[LR_A])
#define LR_RB (v[LR_NB] / v[LR_B])
FN(lr_na,   v[LR_MA] / v[LR_MMA])
FN(lr_nb,   v[LR_MB] / v[LR_MMB])
FN(lr_np,   min2(LR_RA, LR_RB) * v[LR_P])
FN(lr_mp,   v[LR_NP] * v[LR_MMP])
FN(lr_lim,  LR_RA <= LR_RB ? 1 : 2)
FN(lr_left, LR_RA <= LR_RB ? (v[LR_NB] - LR_RA * v[LR_B]) * v[LR_MMB]
                           : (v[LR_NA] - LR_RB * v[LR_A]) * v[LR_MMA])
static const var_t lim_vars[] = {
    {"mA",   "g",     "mA: grams of reactant A", UNK},
    {"MA",   "g/mol", "MA: molar mass of A (GRAPH)", UNK},
    {"a",    "",      "a: coefficient of A (blank 1)", 1},
    {"mB",   "g",     "mB: grams of reactant B", UNK},
    {"MB",   "g/mol", "MB: molar mass of B (GRAPH)", UNK},
    {"b",    "",      "b: coefficient of B (blank 1)", 1},
    {"p",    "",      "p: coefficient of product", 1},
    {"MP",   "g/mol", "MP: molar mass of product", UNK},
    {"nA",   "mol",   "nA: moles of A", UNK},
    {"nB",   "mol",   "nB: moles of B", UNK},
    {"nP",   "mol",   "nP: most product possible", UNK},
    {"mP",   "g",     "mP: theoretical yield", UNK},
    {"lim",  "",      "lim: limiting one: 1 = A, 2 = B", UNK},
    {"left", "g",     "left: grams of the EXCESS left", UNK},
};
static const rule_t lim_rules[] = {
    {LR_NA,   B(LR_MA)|B(LR_MMA), lr_na, "n = m / M", 0},
    {LR_NB,   B(LR_MB)|B(LR_MMB), lr_nb, "n = m / M", 0},
    {LR_NP,   B(LR_NA)|B(LR_A)|B(LR_NB)|B(LR_B)|B(LR_P), lr_np, "nP = p x smaller of nA/a, nB/b", 0},
    {LR_MP,   B(LR_NP)|B(LR_MMP), lr_mp, "mP = nP x MP", 0},
    {LR_LIM,  B(LR_NA)|B(LR_A)|B(LR_NB)|B(LR_B), lr_lim, "limiting = smaller of nA/a, nB/b", 0},
    {LR_LEFT, B(LR_NA)|B(LR_A)|B(LR_NB)|B(LR_B)|B(LR_MMA)|B(LR_MMB), lr_left,
              "left = (n excess - n used) x M", 0},
};

/* ---------- percent yield ------------------------------------------- */
enum { PY_ACT, PY_THEO, PY_PCT };
FN(py_pct,  100 * v[PY_ACT] / v[PY_THEO])
FN(py_act,  v[PY_PCT] * v[PY_THEO] / 100)
FN(py_theo, 100 * v[PY_ACT] / v[PY_PCT])
static const var_t yld_vars[] = {
    {"actual", "g", "actual: what you really got", UNK},
    {"theo",   "g", "theo: theoretical (stoich)", UNK},
    {"pct",    "%", "pct: percent yield", UNK},
};
#define PY_EQ "pct = actual / theo x 100"
static const rule_t yld_rules[] = {
    {PY_PCT,  B(PY_ACT)|B(PY_THEO), py_pct, PY_EQ, 0},
    {PY_ACT,  B(PY_PCT)|B(PY_THEO), py_act, PY_EQ, 0},
    {PY_THEO, B(PY_ACT)|B(PY_PCT), py_theo, PY_EQ, 0},
};

/* ---------- heat: q = m c dT ---------------------------------------- */
enum { HQ_Q, HQ_M, HQ_C, HQ_TI, HQ_TF, HQ_DT, HQ_KJ };
FN(hq_dt,  v[HQ_TF] - v[HQ_TI])
FN(hq_tf,  v[HQ_TI] + v[HQ_DT])
FN(hq_ti,  v[HQ_TF] - v[HQ_DT])
FN(hq_q,   v[HQ_M] * v[HQ_C] * v[HQ_DT])
FN(hq_m,   v[HQ_Q] / (v[HQ_C] * v[HQ_DT]))
FN(hq_c,   v[HQ_Q] / (v[HQ_M] * v[HQ_DT]))
FN(hq_dtq, v[HQ_Q] / (v[HQ_M] * v[HQ_C]))
FN(hq_kj,  v[HQ_Q] / 1000)
FN(hq_qk,  v[HQ_KJ] * 1000)
static const var_t heat_vars[] = {
    {"q",  "J",     "q: heat (+ in, - given off)", UNK},
    {"m",  "g",     "m: mass", UNK},
    {"c",  "J/g C", "c: specific heat (blank=water)", CWATER},
    {"Ti", "deg C", "Ti: starting temperature", UNK},
    {"Tf", "deg C", "Tf: final temperature", UNK},
    {"dT", "deg C", "dT: change = Tf - Ti", UNK},
    {"kJ", "kJ",    "kJ: same heat in kJ", UNK},
};
#define HQ_EQ "q = m c dT"
#define HT_EQ "dT = Tf - Ti"
static const rule_t heat_rules[] = {
    {HQ_Q,  B(HQ_KJ), hq_qk, "q = kJ x 1000", 0},
    {HQ_DT, B(HQ_TI)|B(HQ_TF), hq_dt, HT_EQ, 0},
    {HQ_Q,  B(HQ_M)|B(HQ_C)|B(HQ_DT), hq_q, HQ_EQ, 0},
    {HQ_C,  B(HQ_Q)|B(HQ_M)|B(HQ_DT), hq_c, HQ_EQ, 0},
    {HQ_M,  B(HQ_Q)|B(HQ_C)|B(HQ_DT), hq_m, HQ_EQ, 0},
    {HQ_DT, B(HQ_Q)|B(HQ_M)|B(HQ_C), hq_dtq, HQ_EQ, 0},
    {HQ_TF, B(HQ_TI)|B(HQ_DT), hq_tf, HT_EQ, 0},
    {HQ_TI, B(HQ_TF)|B(HQ_DT), hq_ti, HT_EQ, 0},
    {HQ_KJ, B(HQ_Q), hq_kj, "q = kJ x 1000", 0},
};

/* ---------- calorimetry: hot object into water ---------------------- */
enum { CA_M1, CA_C1, CA_T1, CA_M2, CA_C2, CA_T2, CA_TF, CA_Q };
#define CA_H1 (v[CA_M1] * v[CA_C1])
#define CA_H2 (v[CA_M2] * v[CA_C2])
FN(ca_tf, (CA_H1 * v[CA_T1] + CA_H2 * v[CA_T2]) / (CA_H1 + CA_H2))
FN(ca_c1, CA_H2 * (v[CA_TF] - v[CA_T2]) / (v[CA_M1] * (v[CA_T1] - v[CA_TF])))
FN(ca_m1, CA_H2 * (v[CA_TF] - v[CA_T2]) / (v[CA_C1] * (v[CA_T1] - v[CA_TF])))
FN(ca_t1, v[CA_TF] + CA_H2 * (v[CA_TF] - v[CA_T2]) / CA_H1)
FN(ca_m2, CA_H1 * (v[CA_T1] - v[CA_TF]) / (v[CA_C2] * (v[CA_TF] - v[CA_T2])))
FN(ca_c2, CA_H1 * (v[CA_T1] - v[CA_TF]) / (v[CA_M2] * (v[CA_TF] - v[CA_T2])))
FN(ca_t2, v[CA_TF] - CA_H1 * (v[CA_T1] - v[CA_TF]) / CA_H2)
FN(ca_q,  CA_H2 * (v[CA_TF] - v[CA_T2]))
static const var_t cal_vars[] = {
    {"m1", "g",     "m1: hot object (metal) mass", UNK},
    {"c1", "J/g C", "c1: its specific heat", UNK},
    {"T1", "deg C", "T1: its starting temperature", UNK},
    {"m2", "g",     "m2: water mass (1 mL = 1 g)", UNK},
    {"c2", "J/g C", "c2: water (blank = 4.184)", CWATER},
    {"T2", "deg C", "T2: water starting temperature", UNK},
    {"Tf", "deg C", "Tf: final temperature of both", UNK},
    {"q",  "J",     "q: heat the water gained", UNK},
};
#define CA_EQ "m1 c1 (T1 - Tf) = m2 c2 (Tf - T2)"
static const rule_t cal_rules[] = {
    {CA_TF, B(CA_M1)|B(CA_C1)|B(CA_T1)|B(CA_M2)|B(CA_C2)|B(CA_T2), ca_tf, CA_EQ, 0},
    {CA_C1, B(CA_M1)|B(CA_T1)|B(CA_M2)|B(CA_C2)|B(CA_T2)|B(CA_TF), ca_c1, CA_EQ, 0},
    {CA_M1, B(CA_C1)|B(CA_T1)|B(CA_M2)|B(CA_C2)|B(CA_T2)|B(CA_TF), ca_m1, CA_EQ, 0},
    {CA_T1, B(CA_M1)|B(CA_C1)|B(CA_M2)|B(CA_C2)|B(CA_T2)|B(CA_TF), ca_t1, CA_EQ, 0},
    {CA_M2, B(CA_M1)|B(CA_C1)|B(CA_T1)|B(CA_C2)|B(CA_T2)|B(CA_TF), ca_m2, CA_EQ, 0},
    {CA_C2, B(CA_M1)|B(CA_C1)|B(CA_T1)|B(CA_M2)|B(CA_T2)|B(CA_TF), ca_c2, CA_EQ, 0},
    {CA_T2, B(CA_M1)|B(CA_C1)|B(CA_T1)|B(CA_M2)|B(CA_C2)|B(CA_TF), ca_t2, CA_EQ, 0},
    {CA_Q,  B(CA_M2)|B(CA_C2)|B(CA_TF)|B(CA_T2), ca_q, "q = m2 c2 (Tf - T2)", 0},
};

/* ---------- heat of reaction / phase change ------------------------- */
enum { HR_M, HR_MM, HR_N, HR_DH, HR_Q, HR_CO };
FN(hr_n,  v[HR_M] / v[HR_MM])
FN(hr_m,  v[HR_N] * v[HR_MM])
FN(hr_mm, v[HR_M] / v[HR_N])
FN(hr_q,  v[HR_N] * v[HR_DH] / v[HR_CO])
FN(hr_nq, v[HR_Q] * v[HR_CO] / v[HR_DH])
FN(hr_dh, v[HR_Q] * v[HR_CO] / v[HR_N])
static const var_t hrx_vars[] = {
    {"m",    "g",      "m: grams of the substance", UNK},
    {"M",    "g/mol",  "M: molar mass (GRAPH = table)", UNK},
    {"n",    "mol",    "n: moles of it", UNK},
    {"dH",   "kJ/mol", "dH: from the equation or table", UNK},
    {"q",    "kJ",     "q: heat (- = given off)", UNK},
    {"coef", "",       "coef: its coefficient (blank 1)", 1},
};
#define HR_EQ "q = n x dH / coef"
static const rule_t hrx_rules[] = {
    {HR_N,  B(HR_M)|B(HR_MM), hr_n, "n = m / M", 0},
    {HR_Q,  B(HR_N)|B(HR_DH)|B(HR_CO), hr_q, HR_EQ, 0},
    {HR_N,  B(HR_Q)|B(HR_DH)|B(HR_CO), hr_nq, HR_EQ, 0},
    {HR_DH, B(HR_Q)|B(HR_N)|B(HR_CO), hr_dh, HR_EQ, 0},
    {HR_M,  B(HR_N)|B(HR_MM), hr_m, "n = m / M", 0},
    {HR_MM, B(HR_M)|B(HR_N), hr_mm, "n = m / M", 0},
};

/* ---------- light and photons --------------------------------------- */
enum { LT_NM, LT_L, LT_NU, LT_E, LT_EM };
FN(lt_l,   v[LT_NM] * 1e-9)
FN(lt_nm,  v[LT_L] * 1e9)
FN(lt_nu,  CLIGHT / v[LT_L])
FN(lt_lnu, CLIGHT / v[LT_NU])
FN(lt_e,   HPL * v[LT_NU])
FN(lt_nue, v[LT_E] / HPL)
FN(lt_em,  v[LT_E] * NA / 1000)
FN(lt_eem, v[LT_EM] * 1000 / NA)
static const var_t light_vars[] = {
    {"nm",   "nm",     "nm: wavelength in nanometers", UNK},
    {"lam",  "m",      "lam: wavelength in meters", UNK},
    {"nu",   "Hz",     "nu: frequency (1/s, Hz)", UNK},
    {"E",    "J",      "E: energy of ONE photon", UNK},
    {"Emol", "kJ/mol", "Emol: energy of 1 mol of photons", UNK},
};
#define LT_EQ "c = lam x nu   (c = 2.998E8)"
#define LE_EQ "E = h nu   (h = 6.626E-34)"
static const rule_t light_rules[] = {
    {LT_L,  B(LT_NM), lt_l, "1 nm = 1E-9 m", 0},
    {LT_E,  B(LT_EM), lt_eem, "Emol = E x 6.022E23 / 1000", 0},
    {LT_NU, B(LT_L), lt_nu, LT_EQ, 0},
    {LT_NU, B(LT_E), lt_nue, LE_EQ, 0},
    {LT_L,  B(LT_NU), lt_lnu, LT_EQ, 0},
    {LT_E,  B(LT_NU), lt_e, LE_EQ, 0},
    {LT_NM, B(LT_L), lt_nm, "1 nm = 1E-9 m", 0},
    {LT_EM, B(LT_E), lt_em, "Emol = E x 6.022E23 / 1000", 0},
};

/* ---------- hydrogen atom (Bohr) ------------------------------------ */
enum { HY_NI, HY_NF, HY_E, HY_NM, HY_NU };
FN(hy_e,  RYD * fabs(1 / SQ(v[HY_NF]) - 1 / SQ(v[HY_NI])))
FN(hy_nm, HPL * CLIGHT / v[HY_E] * 1e9)
FN(hy_nu, v[HY_E] / HPL)
FN(hy_enm, HPL * CLIGHT / (v[HY_NM] * 1e-9))
FN(hy_enu, HPL * v[HY_NU])
/* falling down to nf (emission): 1/ni^2 = 1/nf^2 - E/R */
FN(hy_ni, root(1 / (1 / SQ(v[HY_NF]) - v[HY_E] / RYD)))
static const var_t hyd_vars[] = {
    {"ni", "",   "ni: level it starts at", UNK},
    {"nf", "",   "nf: level it ends at", UNK},
    {"E",  "J",  "E: photon energy", UNK},
    {"nm", "nm", "nm: wavelength of the light", UNK},
    {"nu", "Hz", "nu: frequency of the light", UNK},
};
#define HY_EQ "E = 2.18E-18 |1/nf^2 - 1/ni^2|"
static const rule_t hyd_rules[] = {
    {HY_E,  B(HY_NI)|B(HY_NF), hy_e, HY_EQ, 0},
    {HY_E,  B(HY_NM), hy_enm, "E = h c / lam", 0},
    {HY_E,  B(HY_NU), hy_enu, "E = h nu", 0},
    {HY_NM, B(HY_E), hy_nm, "E = h c / lam", 0},
    {HY_NU, B(HY_E), hy_nu, "E = h nu", 0},
    {HY_NI, B(HY_NF)|B(HY_E), hy_ni, HY_EQ " (falls)", R_NOCHK},
};

/* ---------- pH ------------------------------------------------------ */
enum { PH_PH, PH_POH, PH_H, PH_OH };
FN(ph_ph,  v[PH_H] > 0 ? -log10(v[PH_H]) : UNK)
FN(ph_poh, v[PH_OH] > 0 ? -log10(v[PH_OH]) : UNK)
FN(ph_h,   pow(10, -v[PH_PH]))
FN(ph_oh,  pow(10, -v[PH_POH]))
FN(ph_14a, 14 - v[PH_POH])
FN(ph_14b, 14 - v[PH_PH])
static const var_t ph_vars[] = {
    {"pH",  "",  "pH (below 7 = acid)", UNK},
    {"pOH", "",  "pOH", UNK},
    {"H",   "M", "H: [H+] or [H3O+] in mol/L", UNK},
    {"OH",  "M", "OH: [OH-] in mol/L", UNK},
};
static const rule_t ph_rules[] = {
    {PH_PH,  B(PH_H), ph_ph, "pH = -log[H+]", 0},
    {PH_POH, B(PH_OH), ph_poh, "pOH = -log[OH-]", 0},
    {PH_PH,  B(PH_POH), ph_14a, "pH + pOH = 14", 0},
    {PH_POH, B(PH_PH), ph_14b, "pH + pOH = 14", 0},
    {PH_H,   B(PH_PH), ph_h, "[H+] = 10^-pH", 0},
    {PH_OH,  B(PH_POH), ph_oh, "[OH-] = 10^-pOH", 0},
};

/* ---------- titration ----------------------------------------------- */
enum { TI_CA, TI_VA, TI_CB, TI_VB, TI_K, TI_NA };
FN(ti_ca, v[TI_CB] * v[TI_VB] / (v[TI_K] * v[TI_VA]))
FN(ti_va, v[TI_CB] * v[TI_VB] / (v[TI_K] * v[TI_CA]))
FN(ti_cb, v[TI_K] * v[TI_CA] * v[TI_VA] / v[TI_VB])
FN(ti_vb, v[TI_K] * v[TI_CA] * v[TI_VA] / v[TI_CB])
FN(ti_na, v[TI_CA] * v[TI_VA] / 1000)
static const var_t tit_vars[] = {
    {"Ca", "M",   "Ca: acid molarity", UNK},
    {"Va", "mL",  "Va: acid volume", UNK},
    {"Cb", "M",   "Cb: base molarity", UNK},
    {"Vb", "mL",  "Vb: base volume used", UNK},
    {"k",  "",    "k: mol base per mol acid (1)", 1},
    {"na", "mol", "na: moles of acid", UNK},
};
#define TI_EQ "k Ca Va = Cb Vb"
static const rule_t tit_rules[] = {
    {TI_CA, B(TI_CB)|B(TI_VB)|B(TI_K)|B(TI_VA), ti_ca, TI_EQ, 0},
    {TI_VA, B(TI_CB)|B(TI_VB)|B(TI_K)|B(TI_CA), ti_va, TI_EQ, 0},
    {TI_CB, B(TI_K)|B(TI_CA)|B(TI_VA)|B(TI_VB), ti_cb, TI_EQ, 0},
    {TI_VB, B(TI_K)|B(TI_CA)|B(TI_VA)|B(TI_CB), ti_vb, TI_EQ, 0},
    {TI_NA, B(TI_CA)|B(TI_VA), ti_na, "na = Ca x Va / 1000", 0},
};

/* ---------- topic table -------------------------------------------- */
#define N(a) ((uint8_t)(sizeof(a) / sizeof((a)[0])))
#define TOPIC(name, vars, rules, tips) { name, vars, N(vars), rules, N(rules), tips }
static const topic_t topics[] = {
    /* Measurement */
    TOPIC("Density", dens_vars, dens_rules,
        "water: d = 1.00 g/mL\n"
        "1 mL = 1 cm^3   1 L = 1000 mL\n"
        "kg: x1000 (2.5 kg = 2500 g)\n"
        "dropped in water: Vi and Vf\n"
        "cube 2.0 cm on a side: L = 2.0\n"
        "sinks in water if d > 1.00"),
    TOPIC("Percent error", perr_vars, perr_rules,
        "meas = what you got in lab\n"
        "true = book / accepted value\n"
        "always positive (absolute value)"),
    TOPIC("Temperature units", temp_vars, temp_rules,
        "gas laws need kelvin\n"
        "K = C + 273.15 (some books 273)\n"
        "room 25 C, body 37 C, boils 100 C"),
    /* Atoms + moles */
    TOPIC("Grams, moles, particles", mole_vars, mole_rules,
        "molar mass: press GRAPH, build the\n"
        "  formula, then GRAPH again: Use M\n"
        "1 mol = 6.022E23 particles\n"
        "kg: x1000   mg: x0.001\n"
        "atoms of H in H2O: k = 2\n"
        "'molecules' or 'formula units': N"),
    TOPIC("Percent by mass", pmass_vars, pmass_rules,
        "from a formula: press GRAPH, build\n"
        "  it, TRACE lists every element %\n"
        "% of a sample: part / whole x 100"),
    TOPIC("Average atomic mass", avg_vars, avg_rules,
        "abundance 75.78%: type 75.78\n"
        "given as 0.7578? type 75.78\n"
        "2 isotopes: leave m3, f3 blank\n"
        "f2 blank = 100 - f1\n"
        "find an abundance: type avg"),
    TOPIC("Empirical formula", emp_vars, emp_rules,
        "percents given: pretend 100 g,\n"
        "  so 40.0% C means g1 = 40.0\n"
        "A1 = atomic mass (GRAPH: table)\n"
        "r near x.5: multiply all by 2\n"
        "r near x.33 or x.67: by 3\n"
        "r near x.25 or x.75: by 4\n"
        "molecular formula: Memp + Mmol"),
    /* Solutions */
    TOPIC("Molarity", mlr_vars, mlr_rules,
        "'0.50 M' means c = 0.50\n"
        "volume in mL: use the mL box\n"
        "grams needed to make it: type c,\n"
        "  mL and M, read m"),
    TOPIC("Dilution", dil_vars, dil_rules,
        "any volume unit, same for both\n"
        "'diluted TO 500 mL': V2 = 500\n"
        "'ADD 500 mL water': water = 500\n"
        "stock = the strong one"),
    /* Gases */
    TOPIC("Ideal gas PV = nRT", gas_vars, gas_rules,
        "P in atm here; kPa or mmHg have\n"
        "  their own boxes\n"
        "torr = mmHg   1 atm = 760 mmHg\n"
        "Celsius? type it in Tc\n"
        "STP: P = 1 atm, Tc = 0\n"
        "mL: divide by 1000 for V (L)\n"
        "R = 0.08206 L atm / (mol K)"),
    TOPIC("Combined gas law", comb_vars, comb_rules,
        "before = 1, after = 2\n"
        "T must be kelvin: Celsius goes in\n"
        "  T1c and T2c\n"
        "same T: leave both T blank\n"
        "same P: leave both P blank\n"
        "rigid tank (same V): V blank\n"
        "P, V: any unit, same both times"),
    TOPIC("Gas at STP", stp_vars, stp_rules,
        "STP = 0 C and 1 atm\n"
        "1 mol of any gas = 22.4 L\n"
        "not STP? use Ideal gas PV = nRT"),
    TOPIC("Partial pressures", dal_vars, dal_rules,
        "collected over water: P2 = water\n"
        "  vapor (23.8 mmHg at 25 C)\n"
        "mole fraction x1 = n1 / ntot\n"
        "only 2 gases: P3 blank"),
    /* Reactions */
    TOPIC("Stoichiometry (A to B)", sto_vars, sto_rules,
        "balance the equation first\n"
        "2 H2 + O2 -> 2 H2O, H2 to H2O:\n"
        "  a = 2, b = 2\n"
        "mole ratio = b / a\n"
        "have moles: type nA instead of mA\n"
        "gas at STP: VB in liters"),
    TOPIC("Limiting reactant", lim_vars, lim_rules,
        "TWO starting amounts = this one\n"
        "limiting = the one that runs out\n"
        "left = grams of the OTHER one\n"
        "mP = theoretical yield"),
    TOPIC("Percent yield", yld_vars, yld_rules,
        "theoretical: get it from\n"
        "  Stoichiometry or Limiting first\n"
        "over 100%: wet or impure product"),
    /* Heat */
    TOPIC("Heat q = m c dT", heat_vars, heat_rules,
        "water: c = 4.184 (leave blank)\n"
        "ice 2.09  steam 2.01  Al 0.897\n"
        "Cu 0.385  Fe 0.449 (J/g C)\n"
        "heats up: q > 0   cools: q < 0\n"
        "given in kJ: use the kJ box"),
    TOPIC("Calorimetry (mix)", cal_vars, cal_rules,
        "heat lost by hot = heat gained\n"
        "  by the cold water\n"
        "water: 1 mL = 1 g\n"
        "both end at the SAME Tf\n"
        "c1 is the metal's specific heat"),
    TOPIC("Heat of reaction", hrx_vars, hrx_rules,
        "melting: dH = 6.01 (water)\n"
        "boiling: dH = 40.7 (water)\n"
        "freezing, condensing: negative\n"
        "given off (exothermic): dH < 0\n"
        "dH is for '2 H2': coef = 2"),
    /* Light + electrons */
    TOPIC("Light and photons", light_vars, light_rules,
        "visible light: 400 to 700 nm\n"
        "MHz: x1E6   GHz: x1E9\n"
        "kJ per mole of photons: Emol\n"
        "c = 2.998E8, h = 6.626E-34"),
    TOPIC("Hydrogen atom (Bohr)", hyd_vars, hyd_rules,
        "visible lines (Balmer): nf = 2\n"
        "UV lines (Lyman): nf = 1\n"
        "ionize from level n: nf = 1E9\n"
        "ni > nf: emits light (falls)\n"
        "ni < nf: absorbs light"),
    /* Acids + bases */
    TOPIC("pH and pOH", ph_vars, ph_rules,
        "strong acid HCl 0.010 M:\n"
        "  H = 0.010\n"
        "strong base NaOH 0.010 M:\n"
        "  OH = 0.010\n"
        "Ca(OH)2 0.010 M: OH = 0.020\n"
        "pH + pOH = 14 (at 25 C)"),
    TOPIC("Titration", tit_vars, tit_rules,
        "HCl + NaOH: k = 1 (leave blank)\n"
        "H2SO4 + NaOH: k = 2\n"
        "HCl + Ba(OH)2: k = 0.5\n"
        "endpoint, neutralized: this one\n"
        "any volume unit, same for both"),
};
#define NUM_TOPICS N(topics)

typedef struct { const char *name; uint8_t count; } cat_t;
static const cat_t cats[] = {
    {"Measurement", 3},
    {"Atoms + moles", 4},
    {"Solutions", 2},
    {"Gases", 4},
    {"Reactions", 3},
    {"Heat", 3},
    {"Light + electrons", 2},
    {"Acids + bases", 2},
};
#define NUM_CATS N(cats)
static uint8_t cat_first(uint8_t c) {
    uint8_t f = 0;
    for (uint8_t i = 0; i < c; i++) f = (uint8_t)(f + cats[i].count);
    return f;
}

/* ---------- periodic table (molar masses) -------------------------- */
static const char el_sym[119][3] = { "",
    "H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si","P","S","Cl","Ar",
    "K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",
    "Cs","Ba","La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
    "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn",
    "Fr","Ra","Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No","Lr",
    "Rf","Db","Sg","Bh","Hs","Mt","Ds","Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og" };
/* standard atomic weights; mass number of the longest-lived isotope when none */
static const float el_mass[119] = { 0,
    1.008f, 4.0026f, 6.94f, 9.0122f, 10.81f, 12.011f, 14.007f, 15.999f, 18.998f, 20.180f,
    22.990f, 24.305f, 26.982f, 28.085f, 30.974f, 32.06f, 35.45f, 39.948f, 39.098f, 40.078f,
    44.956f, 47.867f, 50.942f, 51.996f, 54.938f, 55.845f, 58.933f, 58.693f, 63.546f, 65.38f,
    69.723f, 72.630f, 74.922f, 78.971f, 79.904f, 83.798f, 85.468f, 87.62f, 88.906f, 91.224f,
    92.906f, 95.95f, 98, 101.07f, 102.91f, 106.42f, 107.87f, 112.41f, 114.82f, 118.71f,
    121.76f, 127.60f, 126.90f, 131.29f, 132.91f, 137.33f, 138.91f, 140.12f, 140.91f, 144.24f,
    145, 150.36f, 151.96f, 157.25f, 158.93f, 162.50f, 164.93f, 167.26f, 168.93f, 173.05f,
    174.97f, 178.49f, 180.95f, 183.84f, 186.21f, 190.23f, 192.22f, 195.08f, 196.97f, 200.59f,
    204.38f, 207.2f, 208.98f, 209, 210, 222, 223, 226, 227, 232.04f,
    231.04f, 238.03f, 237, 244, 243, 247, 247, 251, 252, 257,
    258, 259, 266, 267, 268, 269, 270, 269, 278, 281,
    282, 285, 286, 289, 290, 293, 294, 294 };
#define NUM_EL 118
/* where element z sits: rows 0-6 = periods 1-7, rows 7-8 = the f-block strips */
static void el_pos(uint8_t z, uint8_t *row, uint8_t *col) {
    static const uint8_t start[8] = { 1, 3, 11, 19, 37, 55, 87, 119 };
    uint8_t p = 0, k;
    while (z >= start[p + 1]) p++;
    k = (uint8_t)(z - start[p]);
    *row = p;
    if (p == 0) *col = z == 1 ? 0 : 17;
    else if (p <= 2) *col = (uint8_t)(k < 2 ? k : k + 10);
    else if (p <= 4) *col = k;
    else if (k < 2) *col = k;
    else if (k < 17) { *row = (uint8_t)(p + 2); *col = (uint8_t)(k + 1); }
    else *col = (uint8_t)(k - 14);
}
/* the formula being built: element numbers and how many of each */
#define FMAX 10
static uint8_t f_z[FMAX];
static uint16_t f_n[FMAX];
static uint8_t f_len;
static int8_t f_find(uint8_t z) {
    for (uint8_t i = 0; i < f_len; i++) if (f_z[i] == z) return (int8_t)i;
    return -1;
}
static void f_set(uint8_t z, uint16_t n) {     /* n = 0 removes it */
    int8_t i = f_find(z);
    if (i < 0) {
        if (!n || f_len >= FMAX) return;
        i = (int8_t)f_len++;
        f_z[i] = z;
    }
    if (n) { f_n[i] = n; return; }
    for (uint8_t j = (uint8_t)i; j + 1 < f_len; j++) { f_z[j] = f_z[j + 1]; f_n[j] = f_n[j + 1]; }
    f_len--;
}
static double f_mass(void) {
    double m = 0;
    for (uint8_t i = 0; i < f_len; i++) m += (double)el_mass[f_z[i]] * f_n[i];
    return m;
}
/* molar masses with 3 decimals, trailing zeros dropped: 180.156, 58.44, 18.015 */
static void fmt_mass(double x, char *o) {
    uint32_t t = (uint32_t)(x * 1000 + 0.5), ip = t / 1000, fp = t % 1000;
    char d[12];
    uint8_t n = 0, k = 0;
    do { d[k++] = (char)('0' + ip % 10); ip /= 10; } while (ip);
    while (k) o[n++] = d[--k];
    if (fp) {
        o[n++] = '.';
        o[n++] = (char)('0' + fp / 100); o[n++] = (char)('0' + fp / 10 % 10); o[n++] = (char)('0' + fp % 10);
        while (o[n - 1] == '0') n--;
    }
    o[n] = 0;
}
/* "C6H12O6" */
static void f_name(char *out, uint8_t max) {
    uint8_t n = 0;
    out[0] = 0;
    for (uint8_t i = 0; i < f_len; i++) {
        const char *s = el_sym[f_z[i]];
        char d[6];
        uint8_t k = 0;
        uint16_t c = f_n[i];
        while (*s && n + 1 < max) out[n++] = *s++;
        if (c > 1) {
            while (c) { d[k++] = (char)('0' + c % 10); c /= 10; }
            while (k && n + 1 < max) out[n++] = d[--k];
        }
    }
    out[n] = 0;
}

/* ---------- form state + solver ------------------------------------ */
enum { S_UNK, S_USER, S_SOLVED, S_DEF };
#ifdef HOST_TEST
#define VAL_LEN 30
#else
#define VAL_LEN 11
#endif

static const topic_t *tp;
static char txt[MAXV][VAL_LEN + 1];
static double val[MAXV];
static uint8_t st[MAXV];        /* declared again: the topics above read it */
static uint8_t by[MAXV];        /* rule that produced a solved value */
static int16_t conflict;        /* rule index that disagrees, or -1 */

static double parse_num(const char *s) {
    double m = 0, scale = 1;
    int e = 0, esign = 1;
    bool neg = false, digits = false, frac = false;
    if (*s == '-') { neg = true; s++; }
    for (; *s && *s != 'E'; s++) {
        if (*s == '.') { frac = true; continue; }
        if (*s < '0' || *s > '9') return UNK;
        m = m * 10 + (*s - '0');
        if (frac) scale *= 10;
        digits = true;
    }
    if (!digits) return UNK;
    m /= scale;
    if (*s == 'E') {
        s++;
        if (*s == '-') { esign = -1; s++; }
        if (!*s) return UNK;
        for (; *s; s++) {
            if (*s < '0' || *s > '9') return UNK;
            e = e * 10 + (*s - '0');
        }
        if (e > 37) return UNK;
        m *= pow(10.0, esign * e);
    }
    return neg ? -m : m;
}

/* 4 significant figures, no printf needed */
static void fmt_num(double x, char *out) {
    char dig[12];
    int e, dec, n = 0, i;
    uint32_t s;
    double ax = fabs(x), m;
    bool sci;

    if (!ok(x)) { strcpy(out, "?"); return; }
    if (ax < 1e-30) { strcpy(out, "0"); return; }
    e = (int)floor(log10(ax));
    m = ax / pow(10.0, e);
    if (m >= 9.9995) { m /= 10; e++; }
    if (m < 0.99995) { m *= 10; e--; }
    sci = (e >= 6 || e <= -4);
    if (sci) { dec = 3; s = (uint32_t)(m * 1000 + 0.5); }
    else {
        dec = 3 - e;
        if (dec < 0) dec = 0;
        s = (uint32_t)(ax * pow(10.0, dec) + 0.5);
    }
    if (x < 0) *out++ = '-';
    do { dig[n++] = (char)('0' + s % 10); s /= 10; } while (s || n <= dec);
    for (i = n - 1; i >= 0; i--) {
        *out++ = dig[i];
        if (i == dec && dec > 0) *out++ = '.';
    }
    if (dec > 0) {                       /* trim trailing zeros */
        while (out[-1] == '0') out--;
        if (out[-1] == '.') out--;
    }
    if (sci) {
        *out++ = 'E';
        if (e < 0) { *out++ = '-'; e = -e; }
        if (e >= 10) *out++ = (char)('0' + e / 10);
        *out++ = (char)('0' + e % 10);
    }
    *out = 0;
}

static uint16_t known_mask(void) {
    uint16_t k = 0;
    for (uint8_t i = 0; i < tp->nvars; i++)
        if (st[i] != S_UNK) k |= B(i);
    return k;
}

/* How much a rule's answer moves when each number that went into it changes by
   the fraction `delta`. Used two ways: with a tiny delta, an answer no bigger
   than this is just rounding noise (so it is really 0); with delta = 0.5%, it
   is how far an answer can drift because the givens were typed rounded. */
#define NOISE (sizeof(double) == 4 ? 6e-7 : 1e-11)
static bool nosnap, snapped;
static double sens(const rule_t *R, double x, double delta, bool strict) {
    double s = 0;
    for (uint8_t j = 0; j < tp->nvars; j++) {
        double keep = val[j], y;
        if (!(R->need & B(j)) || keep == 0) continue;
        val[j] = keep * (1 + delta);
        y = R->fn(val);
        val[j] = keep;
        if (!ok(y)) { if (strict) return 1e37; continue; }
        s += fabs(y - x);
    }
    return s;
}

/* Rules are listed in priority order: after any rule fires, start again from
   the top, so the simplest equation that applies is always the one used.
   An answer that is a small difference of big numbers has few good digits, so
   it waits: it is only used when no better-behaved rule can give something. */
static void propagate(void) {
    for (;;) {
        uint16_t k = known_mask();
        int16_t pick = -1, shaky = -1;
        double px = 0, sx = 0;
        bool shaky_zero = false;
        for (uint8_t r = 0; r < tp->nrules && pick < 0; r++) {
            const rule_t *R = &tp->rules[r];
            double x, sn;
            if (st[R->target] != S_UNK || (R->need & k) != R->need) continue;
            x = R->fn(val);
            if (!ok(x)) continue;
            if (x == 0) { pick = r; px = 0; break; }    /* chemistry has real values far below 1E-9 */
            sn = sens(R, x, NOISE, false);
            if (fabs(x) <= sn) {
                /* no bigger than rounding noise: normally a true zero (a - b with a = b),
                   but it can also be a division that cannot be trusted */
                if (!nosnap && shaky < 0) { shaky = r; sx = 0; shaky_zero = true; }
            } else if (sn > 1e-3 * fabs(x)) {
                if (shaky < 0) { shaky = r; sx = x; shaky_zero = false; }
            } else {
                pick = r; px = x;
            }
        }
        if (pick < 0) {
            if (shaky < 0) return;
            pick = shaky; px = sx;
            if (shaky_zero) snapped = true;
        }
        {
            const rule_t *R = &tp->rules[pick];
            if ((R->flags & R_ABS) && (flip & B(R->target))) px = -px;
            val[R->target] = px;
            st[R->target] = S_SOLVED;
            by[R->target] = (uint8_t)pick;
        }
    }
}

/* Blanks that have a usual value (no friction, no spring, force at 90 deg).
   If those blanks are the only unknowns left, one of them is what the problem
   is asking for: keep the first one that the others let us solve. */
static void apply_defaults(void) {
    uint8_t cand[MAXV], nc = 0, real = 0;
    for (uint8_t i = 0; i < tp->nvars; i++) {
        if (st[i] != S_UNK) continue;
        if (ok(tp->vars[i].def)) cand[nc++] = i; else real++;
    }
    if (!nc) return;
    if (!real) {
        double sv[MAXV];
        uint8_t ss[MAXV], sb[MAXV];
        memcpy(sv, val, sizeof sv); memcpy(ss, st, sizeof ss); memcpy(sb, by, sizeof sb);
        for (uint8_t c = 0; c < nc; c++) {
            for (uint8_t o = 0; o < nc; o++) {
                if (o == c) continue;
                val[cand[o]] = tp->vars[cand[o]].def;
                st[cand[o]] = S_DEF;
            }
            propagate();
            if (st[cand[c]] == S_SOLVED) return;
            memcpy(val, sv, sizeof sv); memcpy(st, ss, sizeof ss); memcpy(by, sb, sizeof sb);
        }
    }
    for (uint8_t c = 0; c < nc; c++) {
        if (st[cand[c]] != S_UNK) continue;
        val[cand[c]] = tp->vars[cand[c]].def;
        st[cand[c]] = S_DEF;
        propagate();
    }
}

static void solve_once(void) {
    uint16_t k;
    for (uint8_t i = 0; i < tp->nvars; i++) {
        double x = txt[i][0] ? parse_num(txt[i]) : UNK;
        st[i] = ok(x) ? S_USER : S_UNK;
        val[i] = ok(x) ? x : 0;
    }
    propagate();
    apply_defaults();
    /* do the givens agree with each other? */
    conflict = -1;
    k = known_mask();
    for (uint8_t r = 0; r < tp->nrules; r++) {
        const rule_t *R = &tp->rules[r];
        double x, y, tol;
        if ((R->flags & R_NOCHK) || !(k & B(R->target)) || (R->need & k) != R->need) continue;
        if (st[R->target] == S_SOLVED && by[R->target] == r) continue;
        x = R->fn(val);
        y = val[R->target];
        if (!ok(x)) continue;
        if (R->flags & R_ABS) { x = fabs(x); y = fabs(y); }
        /* allow for the givens being typed rounded to about 3 figures */
        tol = 5e-3 * (fabs(x) > fabs(y) ? fabs(x) : fabs(y)) + 1e-30;
        if (fabs(x - y) <= tol) continue;
        tol += sens(R, R->fn(val), 5e-3, true);
        if (fabs(x - y) > tol) { conflict = r; break; }
    }
}

/* A square root has two signs. The rules pick the usual one (moving forward,
   launched upward). If that choice contradicts something else that was typed,
   try the other sign before deciding the numbers disagree. */
static void solve(void) {
    uint16_t amb = 0;
    flip = 0;
    nosnap = snapped = false;
    solve_once();
    if (conflict < 0) return;
    if (snapped) {              /* a "zero" that then contradicts the givens was not a zero */
        nosnap = true;
        solve_once();
        if (conflict < 0) return;
        nosnap = false;
        solve_once();
    }
    for (uint8_t i = 0; i < tp->nvars; i++)
        if (st[i] == S_SOLVED && (tp->rules[by[i]].flags & (R_ABS | R_ALT))) amb |= B(i);
    for (uint8_t i = 0; i < tp->nvars; i++) {
        if (!(amb & B(i))) continue;
        flip = B(i);
        solve_once();
        if (conflict < 0) return;
    }
    flip = 0;
    solve_once();
}

/* ---------- text builder (steps, tips, prompts) -------------------- */
#define MAX_TEXT 6144
static char text[MAX_TEXT];
static uint16_t text_len;

static void tput(const char *s) {
    while (*s && text_len < MAX_TEXT - 1) text[text_len++] = *s++;
    text[text_len] = 0;
}
static void tnum(double x) { char b[20]; fmt_num(x, b); tput(b); }
static void tvar(uint8_t i) {
    tput(tp->vars[i].sym); tput(" = "); tnum(val[i]);
    if (tp->vars[i].unit[0]) { tput(" "); tput(tp->vars[i].unit); }
}

static uint16_t step_done;
static uint8_t step_no;
static void step_walk(uint8_t i) {
    const rule_t *R;
    if (st[i] != S_SOLVED || (step_done & B(i))) return;
    step_done |= B(i);
    R = &tp->rules[by[i]];
    for (uint8_t j = 0; j < tp->nvars; j++)
        if (R->need & B(j)) step_walk(j);
    step_no++;
    tput("\n");
    { char n[4] = { (char)('0' + step_no / 10), (char)('0' + step_no % 10), 0, 0 };
      tput(step_no >= 10 ? n : n + 1); }
    tput(") "); tput(R->eq); tput("\n   "); tvar(i); tput("\n");
}

static uint8_t bits(uint16_t m) { uint8_t n = 0; for (; m; m &= (uint16_t)(m - 1)) n++; return n; }
static void put_syms(uint16_t m) {
    bool first = true;
    for (uint8_t j = 0; j < tp->nvars; j++) {
        if (!(m & B(j))) continue;
        if (!first) tput(", ");
        tput(tp->vars[j].sym);
        first = false;
    }
}
/* the smallest set of extra givens that would unlock variable i (0 = no way) */
static uint16_t easiest_missing(uint8_t i) {
    uint16_t k = known_mask(), best = 0;
    uint8_t bestn = 255;
    for (uint8_t r = 0; r < tp->nrules; r++) {
        uint16_t miss = (uint16_t)(tp->rules[r].need & ~k);
        uint8_t n = bits(miss);
        if (tp->rules[r].target == i && n && n < bestn) { bestn = n; best = miss; }
    }
    return best;
}
/* every different way to unlock variable i, fewest extra numbers first */
static void list_ways(uint8_t i) {
    uint16_t k = known_mask(), seen[6];
    uint8_t nseen = 0;
    for (uint8_t n = 1; n <= 4 && nseen < 6; n++) {
        for (uint8_t r = 0; r < tp->nrules && nseen < 6; r++) {
            uint16_t miss = (uint16_t)(tp->rules[r].need & ~k);
            bool dup = false;
            if (tp->rules[r].target != i || bits(miss) != n) continue;
            for (uint8_t q = 0; q < nseen; q++) if (seen[q] == miss) dup = true;
            if (dup) continue;
            seen[nseen++] = miss;
            tput("  enter "); put_syms(miss); tput("\n    ("); tput(tp->rules[r].eq); tput(")\n");
        }
    }
    if (!nseen) tput("  (this one is only an input)\n");
}

/* the equations behind a topic, each listed once */
static void put_formulas(const topic_t *t) {
    for (uint8_t r = 0; r < t->nrules; r++) {
        bool dup = false;
        for (uint8_t q = 0; q < r; q++)
            if (!strcmp(t->rules[q].eq, t->rules[r].eq)) { dup = true; break; }
        if (dup) continue;
        tput("  "); tput(t->rules[r].eq); tput("\n");
    }
}
static void build_formula_sheet(void) {
    uint8_t ti = 0;
    text_len = 0;
    tput("NA = 6.022E23 /mol\nR = 0.08206 L atm / (mol K)\nh = 6.626E-34 J s   c = 2.998E8 m/s\nSTP: 0 C, 1 atm, 22.4 L/mol\nwater: c = 4.184 J/g C, d = 1.00 g/mL\nK = C + 273.15   1 atm = 760 mmHg\n");
    for (uint8_t c = 0; c < N(cats); c++) {
        tput("\n=== "); tput(cats[c].name); tput(" ===\n");
        for (uint8_t i = 0; i < cats[c].count; i++, ti++) {
            tput(topics[ti].name); tput("\n");
            put_formulas(&topics[ti]);
        }
    }
}

/* a remark about what the numbers on screen mean, when there is something to say */
static const char *situation_note(void) {
#define KN(i) (st[i] != S_UNK)
    if (tp->vars == temp_vars) {
        if (KN(TM_K) && val[TM_K] < 0) return "below 0 K is impossible: check it";
    } else if (tp->vars == dens_vars) {
        if (KN(DN_D) && st[DN_D] == S_SOLVED)
            return val[DN_D] > 1.0 ? "d > 1.00 g/mL: sinks in water" : "d < 1.00 g/mL: floats on water";
    } else if (tp->vars == emp_vars) {
        /* smallest whole-number multiplier that makes every ratio whole */
        if (KN(EF_R1) && KN(EF_R2)) {
            static const char *const by[] = { "ratios are whole numbers: done",
                "multiply every r by 2", "multiply every r by 3", "multiply every r by 4",
                "multiply every r by 5", "multiply every r by 6" };
            for (uint8_t k = 1; k <= 6; k++) {
                bool whole = true;
                for (uint8_t i = EF_R1; i <= EF_R3; i++) {
                    double r;
                    if (!KN(i)) continue;
                    r = val[i] * k;
                    if (fabs(r - floor(r + 0.5)) > 0.1) whole = false;
                }
                if (whole) return by[k - 1];
            }
            return "ratios not close to whole: check";
        }
    } else if (tp->vars == lim_vars) {
        if (KN(LR_LIM)) return val[LR_LIM] < 1.5 ? "A runs out first: A is limiting" : "B runs out first: B is limiting";
    } else if (tp->vars == yld_vars) {
        if (KN(PY_PCT) && val[PY_PCT] > 100) return "over 100%: product wet or impure?";
    } else if (tp->vars == heat_vars) {
        if (KN(HQ_Q) && st[HQ_Q] == S_SOLVED)
            return val[HQ_Q] > 0 ? "q > 0: heat absorbed (warms up)" : "q < 0: heat given off (cools)";
    } else if (tp->vars == cal_vars) {
        if (KN(CA_TF) && KN(CA_T1) && KN(CA_T2) &&
            (val[CA_TF] > (val[CA_T1] > val[CA_T2] ? val[CA_T1] : val[CA_T2]) + 1e-6 ||
             val[CA_TF] < (val[CA_T1] < val[CA_T2] ? val[CA_T1] : val[CA_T2]) - 1e-6))
            return "Tf must be between T1 and T2";
    } else if (tp->vars == hrx_vars) {
        if (KN(HR_Q)) return val[HR_Q] < 0 ? "q < 0: heat given off (exothermic)" : "q > 0: heat taken in (endothermic)";
    } else if (tp->vars == hyd_vars) {
        if (KN(HY_NI) && KN(HY_NF))
            return val[HY_NI] > val[HY_NF] ? "ni > nf: falls down, EMITS light" : "ni < nf: jumps up, ABSORBS light";
    } else if (tp->vars == ph_vars) {
        if (KN(PH_PH)) return val[PH_PH] < 6.995 ? "pH < 7: acidic" : val[PH_PH] > 7.005 ? "pH > 7: basic" : "pH = 7: neutral";
    } else if (tp->vars == comb_vars) {
        bool t = !KN(CG_T1) && !KN(CG_T2), p = !KN(CG_P1) && !KN(CG_P2), v = !KN(CG_V1) && !KN(CG_V2);
        if (t && !p && !v) return "T blank: held constant (Boyle)";
        if (p && !t && !v) return "P blank: held constant (Charles)";
        if (v && !t && !p) return "V blank: held constant (Gay-Lussac)";
    }
    return NULL;
}

static void build_steps(uint8_t target) {
    bool any = false;
    text_len = 0;
    tput("GIVEN\n");
    for (uint8_t i = 0; i < tp->nvars; i++) {
        if (st[i] == S_USER) { tput("   "); tvar(i); tput("\n"); any = true; }
        else if (st[i] == S_DEF) { tput("   "); tvar(i); tput(" (assumed)\n"); }
    }
    if (!any) tput("   (nothing entered yet)\n");
    step_done = 0;
    step_no = 0;
    if (st[target] == S_SOLVED) {
        step_walk(target);
        tput("\nANSWER\n   "); tvar(target); tput("\n");
    } else if (st[target] == S_UNK) {
        tput("\nNot enough yet to get ");
        tput(tp->vars[target].sym);
        tput(".\nAny ONE of these would do it:\n");
        list_ways(target);
    } else {
        tput("\n"); tput(tp->vars[target].sym);
        tput(" is a number you typed in.\nMove to a solved (colored) value and press ENTER to see its steps.\n");
    }
    if (situation_note()) { tput("\nNOTE\n   "); tput(situation_note()); tput("\n"); }
    if (conflict >= 0) {
        tput("\nWARNING: your numbers do not agree with\n   ");
        tput(tp->rules[conflict].eq);
        tput("\nCheck signs and what you typed.\n");
    }
}

/* ================================================================== */
/* "Help me pick": a few yes/no questions -> the right topic, then the   */
/* boxes to fill one at a time, each with what to look for in the words  */
/* ================================================================== */
#if 1
#define HAVE_WIZ 1
typedef struct { uint8_t var; const char *ask; const char *look; } wstep_t;
typedef struct {
    const var_t *vars;      /* which topic */
    const char *pre;        /* values the answers already told us, "e=0 v1=0" */
    uint8_t answer;         /* box that holds the answer, or WANY = first solved one */
    const char *done;       /* shown under the answer */
    const wstep_t *steps;
    uint8_t nsteps;
} wplan_t;
typedef struct { const char *label; int8_t go; } wopt_t;   /* go >= 0: next question, < 0: a plan */
typedef struct { const char *q; const char *hint; uint8_t nopt; wopt_t opt[4]; } wq_t;
#define WANY 0xFF
#define WALL 0xFE       /* every box ends up filled: nothing single to point at */
#define WSTEPS(a) a, N(a)
#define GO(p) ((int8_t)(-1 - (p)))

#define LK_G    "grams; kg: x1000   mg: x0.001"
#define LK_MM   "builds the formula; GRAPH = Use M"
#define LK_ASK  "the one it asks for: just ENTER"
static const wstep_t ws_g2mol[] = {
    {MO_M,  "Grams it gives", LK_G},
    {MO_MM, "Molar mass? GRAPH = table", LK_MM},
};
static const wstep_t ws_mol2g[] = {
    {MO_N,  "Moles it gives", "the number in front of 'mol'"},
    {MO_MM, "Molar mass? GRAPH = table", LK_MM},
};
static const wstep_t ws_part[] = {
    {MO_P,  "How many molecules or atoms", "6.02 x 10^23: type 6.02E23"},
    {MO_MM, "Molar mass? GRAPH = table", LK_MM},
};
static const wstep_t ws_molar[] = {
    {ML_M,  "Grams of solute? (asked: ENTER)", LK_G},
    {ML_MM, "Molar mass? GRAPH = table", LK_MM},
    {ML_ML, "Volume of solution in mL", "liters? x1000 (0.250 L = 250)"},
    {ML_C,  "Molarity? (asked: ENTER)", "'0.50 M' means 0.50"},
};
static const wstep_t ws_dil[] = {
    {DL_C1, "Molarity of the STOCK", "the strong one; asked: ENTER"},
    {DL_V1, "Volume of stock used", LK_ASK},
    {DL_C2, "Molarity after diluting", LK_ASK},
    {DL_V2, "Final TOTAL volume", "'diluted to 500 mL' = 500"},
    {DL_W,  "Water added? (none: ENTER)", "'add 250 mL of water' = 250"},
};
static const wstep_t ws_tit[] = {
    {TI_CA, "Acid molarity", LK_ASK},
    {TI_VA, "Acid volume", "mL (any unit, same for both)"},
    {TI_CB, "Base molarity", LK_ASK},
    {TI_VB, "Base volume used", LK_ASK},
    {TI_K,  "Base per acid (H2SO4 = 2)", "HCl + NaOH: just ENTER (1)"},
};
static const wstep_t ws_comb[] = {
    {CG_P1, "Pressure BEFORE", "any unit; same P: ENTER"},
    {CG_V1, "Volume BEFORE", "any unit; same V: ENTER"},
    {CG_T1, "Temp BEFORE in kelvin", "Celsius? ENTER; same T: ENTER"},
    {CG_C1, "Temp BEFORE in Celsius", "same T: ENTER"},
    {CG_P2, "Pressure AFTER", LK_ASK},
    {CG_V2, "Volume AFTER", LK_ASK},
    {CG_T2, "Temp AFTER in kelvin", "Celsius? ENTER"},
    {CG_C2, "Temp AFTER in Celsius", LK_ASK},
};
static const wstep_t ws_ideal[] = {
    {IG_P,    "Pressure in atm", "kPa or mmHg? ENTER, next boxes"},
    {IG_KPA,  "Pressure in kPa", "mmHg or torr? ENTER"},
    {IG_MMHG, "Pressure in mmHg (torr)", LK_ASK},
    {IG_V,    "Volume in liters", "mL: divide by 1000"},
    {IG_T,    "Temperature in kelvin", "Celsius? ENTER"},
    {IG_TC,   "Temperature in Celsius", "STP: 0"},
    {IG_N,    "Moles of gas", "grams given instead? ENTER"},
    {IG_M,    "Grams of gas (none: ENTER)", LK_ASK},
    {IG_MM,   "Molar mass? GRAPH = table", LK_MM},
};
static const wstep_t ws_stp[] = {
    {SP_N,  "Moles of gas", LK_ASK},
    {SP_V,  "Liters at STP", LK_ASK},
    {SP_M,  "Grams of gas", LK_ASK},
    {SP_MM, "Molar mass? GRAPH = table", LK_MM},
};
static const wstep_t ws_sto[] = {
    {SO_MA,  "Grams of what you HAVE", "given moles instead? ENTER"},
    {SO_NA,  "Moles of what you HAVE", "had grams? just ENTER"},
    {SO_MMA, "Its molar mass? GRAPH = table", LK_MM},
    {SO_A,   "Its coefficient (1: ENTER)", "number in front in the equation"},
    {SO_B,   "WANTED one's coefficient", "1: just ENTER"},
    {SO_MMB, "WANTED one's molar mass? GRAPH", LK_MM},
};
static const wstep_t ws_lim[] = {
    {LR_MA,  "Grams of reactant A", LK_G},
    {LR_MMA, "A's molar mass? GRAPH = table", LK_MM},
    {LR_A,   "A's coefficient (1: ENTER)", "number in front in the equation"},
    {LR_MB,  "Grams of reactant B", LK_G},
    {LR_MMB, "B's molar mass? GRAPH = table", LK_MM},
    {LR_B,   "B's coefficient (1: ENTER)", "number in front in the equation"},
    {LR_P,   "Product's coefficient (1: ENTER)", "number in front in the equation"},
    {LR_MMP, "Product's molar mass? GRAPH", LK_MM},
};
static const wstep_t ws_yld[] = {
    {PY_ACT,  "Grams you actually got", "'actual', 'recovered', 'isolated'"},
    {PY_THEO, "Theoretical yield (grams)", "from Stoichiometry first"},
};
static const wstep_t ws_heat[] = {
    {HQ_M,  "Mass in grams", LK_G},
    {HQ_C,  "Specific heat (water: ENTER)", "metals: see Tips (ZOOM)"},
    {HQ_TI, "Starting temperature", LK_ASK},
    {HQ_TF, "Final temperature", LK_ASK},
    {HQ_Q,  "Heat in joules", "kJ: x1000; asked: ENTER"},
};
static const wstep_t ws_cal[] = {
    {CA_M1, "Mass of the hot object", LK_G},
    {CA_C1, "Its specific heat", LK_ASK},
    {CA_T1, "Its starting temperature", "deg C"},
    {CA_M2, "Mass of the water", "1 mL of water = 1 g"},
    {CA_T2, "Water's starting temperature", "deg C"},
    {CA_TF, "Final temperature", LK_ASK},
};
static const wstep_t ws_ph[] = {
    {PH_H,   "[H+] in mol/L", "pH given instead? ENTER"},
    {PH_PH,  "pH", "pOH or [OH-] given? ENTER"},
    {PH_OH,  "[OH-] in mol/L", "strong base: its molarity"},
    {PH_POH, "pOH", LK_ASK},
};
enum { WP_G2MOL, WP_MOL2G, WP_PART, WP_MOLAR, WP_DIL, WP_TIT, WP_COMB, WP_IDEAL,
       WP_STP, WP_STO, WP_LIM, WP_YLD, WP_HEAT, WP_CAL, WP_PH };
#define D_COLOR "the answer is the colored box"
static const wplan_t wplans[] = {
    {mole_vars, NULL, WANY,   "n = moles, N = molecules", WSTEPS(ws_g2mol)},
    {mole_vars, NULL, WANY,   "m = grams, N = molecules", WSTEPS(ws_mol2g)},
    {mole_vars, NULL, WANY,   "m = grams, n = moles", WSTEPS(ws_part)},
    {mlr_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_molar)},
    {dil_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_dil)},
    {tit_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_tit)},
    {comb_vars, NULL, WANY,   D_COLOR, WSTEPS(ws_comb)},
    {gas_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_ideal)},
    {stp_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_stp)},
    {sto_vars,  NULL, SO_MB,  "mB = grams; nB moles, VB liters", WSTEPS(ws_sto)},
    {lim_vars,  NULL, LR_MP,  "theoretical yield; lim, left below", WSTEPS(ws_lim)},
    {yld_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_yld)},
    {heat_vars, NULL, WANY,   D_COLOR, WSTEPS(ws_heat)},
    {cal_vars,  NULL, WANY,   D_COLOR, WSTEPS(ws_cal)},
    {ph_vars,   NULL, WALL,   "DONE: pH, pOH, [H+] and [OH-] filled", WSTEPS(ws_ph)},
};
static const wq_t wqs[] = {
    /* 0 */ {"What is the problem about?", "pick the closest match", 4,
             {{"Grams, moles, or molecules", 1}, {"Solutions: M, dilute, titrate", 2},
              {"A gas (pressure, volume, temp)", 3}, {"Reaction, heat, or pH", 4}}},
    /* 1 */ {"What does it GIVE you?", "look at the unit after the number", 3,
             {{"Grams (g)", GO(WP_G2MOL)}, {"Moles (mol)", GO(WP_MOL2G)},
              {"A number of atoms or molecules", GO(WP_PART)}}},
    /* 2 */ {"Which kind of solution problem?", "molarity: M = mol / L", 3,
             {{"Molarity, or grams to make it", GO(WP_MOLAR)}, {"Diluting (adding water)", GO(WP_DIL)},
              {"Titration (acid + base)", GO(WP_TIT)}}},
    /* 3 */ {"Does the gas CHANGE?", "before AND after = it changes", 3,
             {{"Yes, before and after", GO(WP_COMB)}, {"No, one set of P, V, T", GO(WP_IDEAL)},
              {"It says STP", GO(WP_STP)}}},
    /* 4 */ {"Which one is it?", "pick the closest match", 4,
             {{"Grams of one thing to another", 5}, {"Heat (q, specific heat)", 6},
              {"pH, [H+], or [OH-]", GO(WP_PH)}, {"Percent yield", GO(WP_YLD)}}},
    /* 5 */ {"Amounts of TWO reactants?", "two starting amounts = limiting", 2,
             {{"No, just one", GO(WP_STO)}, {"Yes, two", GO(WP_LIM)}}},
    /* 6 */ {"What kind of heat problem?", "dropped into water = mixing", 2,
             {{"Heating or cooling one thing", GO(WP_HEAT)}, {"Hot object dropped in water", GO(WP_CAL)}}},
};

static const wplan_t *guide;    /* walking through a plan's boxes, or NULL */
static uint8_t gstep;
static uint16_t gblank;         /* boxes the user left blank on the walk */
static int wiz_topic(const wplan_t *p) {
    for (unsigned i = 0; i < NUM_TOPICS; i++) if (topics[i].vars == p->vars) return (int)i;
    return -1;
}
/* "e=0 v1=0": copy each value into its box */
static void wiz_preset(const char *s) {
    while (s && *s) {
        char sym[8];
        uint8_t n = 0, i;
        while (*s == ' ') s++;
        while (*s && *s != '=' && n < sizeof sym - 1) sym[n++] = *s++;
        sym[n] = 0;
        if (*s++ != '=') return;
        for (i = 0; i < tp->nvars && strcmp(tp->vars[i].sym, sym); i++) {}
        n = 0;
        while (*s && *s != ' ') { if (i < tp->nvars && n < VAL_LEN) txt[i][n++] = *s; s++; }
        if (i < tp->nvars) txt[i][n] = 0;
    }
}
/* a box the numbers already filled in needs no typing */
static void wiz_skip(void) {
    while (gstep < guide->nsteps && st[guide->steps[gstep].var] == S_SOLVED) gstep++;
}
static bool wiz_begin(const wplan_t *p) {
    int t = wiz_topic(p);
    if (t < 0) return false;
    tp = &topics[t];
    memset(txt, 0, sizeof txt);
    wiz_preset(p->pre);
    solve();
    guide = p;
    gstep = 0;
    gblank = 0;
    wiz_skip();
    return true;
}
/* the current box is done (typed or left blank): true when every box is done */
static bool wiz_next(void) {
    if (gstep < guide->nsteps && !txt[guide->steps[gstep].var][0]) gblank |= B(guide->steps[gstep].var);
    solve();
    if (gstep < guide->nsteps) gstep++;
    wiz_skip();
    return gstep >= guide->nsteps;
}
/* solved only by putting another typed box in a different unit (Tc -> T) */
static bool wiz_twin(uint8_t i) {
    const rule_t *R = &tp->rules[by[i]];
    if (bits(R->need) != 1) return false;
    for (uint8_t j = 0; j < tp->nvars; j++) if (R->need & B(j)) return st[j] == S_USER;
    return false;
}
/* the box that holds the answer: one the user skipped that got solved,
   else a box of the walk that solved itself, else the first solved one */
static uint8_t wiz_answer(const wplan_t *p) {
    uint8_t k, v;
    if (p->answer < WALL) return p->answer;
    if (p->answer == WANY) {
        for (k = 0; k < p->nsteps; k++) { v = p->steps[k].var; if ((gblank & B(v)) && st[v] == S_SOLVED && !wiz_twin(v)) return v; }
        for (k = 0; k < p->nsteps; k++) { v = p->steps[k].var; if (st[v] == S_SOLVED && !wiz_twin(v)) return v; }
        for (v = 0; v < tp->nvars; v++) if (st[v] == S_SOLVED) return v;
    }
    return p->steps[0].var;
}
#endif

/* ================================================================== */
#ifdef HOST_TEST
/* ================================================================== */
#include <stdio.h>
#include <stdlib.h>

static int fails = 0;
static int find_topic(const char *name) {
    for (unsigned i = 0; i < NUM_TOPICS; i++) if (!strcmp(topics[i].name, name)) return (int)i;
    printf("NO TOPIC %s\n", name); exit(1);
}
static int find_var(const char *sym) {
    for (int i = 0; i < tp->nvars; i++) if (!strcmp(tp->vars[i].sym, sym)) return i;
    printf("NO VAR %s in %s\n", sym, tp->name); exit(1);
}
/* spec: "v0=20 th=30 | H=5.102 R=35.35" */
static void check(const char *topic, const char *spec) {
    char buf[300], *tok;
    bool expect = false;
    tp = &topics[find_topic(topic)];
    memset(txt, 0, sizeof txt);
    strcpy(buf, spec);
    for (tok = strtok(buf, " "); tok; tok = strtok(NULL, " ")) {
        if (!strcmp(tok, "|")) { solve(); expect = true; continue; }
        char *eq = strchr(tok, '='); *eq = 0;
        int i = find_var(tok);
        if (!expect) { strncpy(txt[i], eq + 1, VAL_LEN); continue; }
        double want = atof(eq + 1);
        bool good = st[i] != S_UNK && fabs(val[i] - want) <= 3e-4 * fabs(want) + 1e-30;
        if (!good) { fails++; printf("FAIL %-26s %s: got %g (state %d) want %g   [%s]\n", topic, tok, val[i], st[i], want, spec); }
    }
    if (conflict >= 0) { fails++; printf("CONFLICT %-22s %s   [%s]\n", topic, tp->rules[conflict].eq, spec); }
}
/* every rule of a topic must agree once everything is known */
static void full(const char *topic, const char *spec) {
    char buf[300]; snprintf(buf, sizeof buf, "%s |", spec);
    check(topic, buf);
    for (int i = 0; i < tp->nvars; i++)
        if (st[i] == S_UNK) { fails++; printf("UNSOLVED %-22s %s   [%s]\n", topic, tp->vars[i].sym, spec); }
    double keep[MAXV]; memcpy(keep, val, sizeof keep);
    /* every rule must reproduce the full set of values */
    for (int r = 0; r < tp->nrules; r++) {
        const rule_t *R = &tp->rules[r];
        if (R->flags & R_NOCHK) continue;
        double x = R->fn(keep), y = keep[R->target];
        if (!ok(x)) continue;       /* rule not applicable for these numbers */
        if (R->flags & R_ABS) { x = fabs(x); y = fabs(y); }
        if (fabs(x - y) > 1e-6 * (fabs(x) + fabs(y)) + 1e-9) {
            fails++; printf("RULE MISMATCH %-18s rule %d (%s -> %s): %g vs %g\n", topic, r, R->eq, tp->vars[R->target].sym, x, y);
        }
    }
}

/* "--solve": read "Topic|sym=val sym=val" lines, print everything known */
static int solve_lines(void) {
    char line[1024];
    while (fgets(line, sizeof line, stdin)) {
        char *bar = strchr(line, '|'), *tok;
        if (!bar) continue;
        *bar = 0;
        tp = &topics[find_topic(line)];
        memset(txt, 0, sizeof txt);
        for (tok = strtok(bar + 1, " \n"); tok; tok = strtok(NULL, " \n")) {
            char *eq = strchr(tok, '='); *eq = 0;
            strncpy(txt[find_var(tok)], eq + 1, VAL_LEN);
        }
        solve();
        printf("%d|", conflict);
        for (int i = 0; i < tp->nvars; i++)
            if (st[i] != S_UNK) printf("%s=%.17g:%d ", tp->vars[i].sym, val[i], st[i]);
        printf("\n");
    }
    return 0;
}
static int list_topics(void) {
    for (unsigned c = 0; c < NUM_CATS; c++)
        for (int i = 0; i < cats[c].count; i++) {
            const topic_t *t = &topics[cat_first(c) + i];
            printf("%s|%s|", cats[c].name, t->name);
            for (int j = 0; j < t->nvars; j++) printf("%s%s", t->vars[j].sym, ok(t->vars[j].def) ? "*" : ""), printf(j + 1 < t->nvars ? "," : "\n");
        }
    return 0;
}

/* "Help me pick": the question tree is sound and every path gives the answer */
static void wiz_len(const char *what, const char *t, size_t max) {
    if (t && strlen(t) > max) { fails++; printf("WIZ TEXT too long (%u > %u) %s: %s\n", (unsigned)strlen(t), (unsigned)max, what, t); }
}
static void wiz_walk(int plan, const char *inputs, const char *sym, double want) {
    const wplan_t *p = &wplans[plan];
    char buf[200], *tok;
    uint8_t a;
    strncpy(buf, inputs, sizeof buf - 1); buf[sizeof buf - 1] = 0;
    tok = strtok(buf, ",");
    if (!wiz_begin(p)) { fails++; printf("WIZ plan %d: no topic\n", plan); return; }
    while (gstep < guide->nsteps) {
        if (!tok) { fails++; printf("WIZ plan %d: ran out of inputs at box %s\n", plan, tp->vars[guide->steps[gstep].var].sym); return; }
        if (strcmp(tok, "_")) strncpy(txt[guide->steps[gstep].var], tok, VAL_LEN);
        tok = strtok(NULL, ",");
        if (wiz_next()) break;
    }
    if (tok) { fails++; printf("WIZ plan %d: input '%s' left over (box was skipped?)\n", plan, tok); }
    a = p->answer == WALL ? (uint8_t)find_var(sym) : wiz_answer(p);
    if (p->answer == WALL && st[a] != S_UNK && fabs(val[a] - want) <= 2e-3 * fabs(want) + 1e-30) { guide = NULL; return; }
    if (strcmp(tp->vars[a].sym, sym) || st[a] != S_SOLVED || fabs(val[a] - want) > 2e-3 * fabs(want) + 1e-30) {
        fails++; printf("WIZ plan %d: answer %s = %g (state %d), want %s = %g\n", plan, tp->vars[a].sym, val[a], st[a], sym, want);
    }
    if (conflict >= 0) { fails++; printf("WIZ plan %d: numbers disagree: %s\n", plan, tp->rules[conflict].eq); }
    guide = NULL;
}
static void wiz_tests(void) {
    bool reach[N(wplans)] = { false }, seen[N(wqs)] = { false };
    uint8_t todo[N(wqs)], nt = 0;
    todo[nt++] = 0; seen[0] = true;
    while (nt) {                         /* walk every branch from the first question */
        const wq_t *q = &wqs[todo[--nt]];
        wiz_len("question", q->q, 38); wiz_len("hint", q->hint, 38);
        if (q->nopt < 2 || q->nopt > 4) { fails++; printf("WIZ question with %d options\n", q->nopt); }
        for (uint8_t o = 0; o < q->nopt; o++) {
            int8_t go = q->opt[o].go;
            wiz_len("option", q->opt[o].label, 34);
            if (go >= 0) {
                if (go >= (int8_t)N(wqs)) { fails++; printf("WIZ bad question %d\n", go); continue; }
                if (!seen[go]) { seen[go] = true; todo[nt++] = (uint8_t)go; }
            } else if (-1 - go >= (int)N(wplans)) { fails++; printf("WIZ bad plan %d\n", -1 - go); }
            else reach[-1 - go] = true;
        }
    }
    for (unsigned i = 0; i < N(wqs); i++) if (!seen[i]) { fails++; printf("WIZ question %u cannot be reached\n", i); }
    for (unsigned pl = 0; pl < N(wplans); pl++) {
        const wplan_t *p = &wplans[pl];
        int t = wiz_topic(p);
        if (!reach[pl]) { fails++; printf("WIZ plan %u cannot be reached\n", pl); }
        if (t < 0) { fails++; printf("WIZ plan %u: topic missing\n", pl); continue; }
        if (p->nsteps > 9) { fails++; printf("WIZ plan %u: more than 9 boxes\n", pl); }
        wiz_len("done", p->done, 38);
        if (p->answer < WALL && p->answer >= topics[t].nvars) { fails++; printf("WIZ plan %u: bad answer box\n", pl); }
        for (uint8_t k = 0; k < p->nsteps; k++) {
            wiz_len("ask", p->steps[k].ask, 34); wiz_len("look", p->steps[k].look, 38);
            if (p->steps[k].var >= topics[t].nvars) { fails++; printf("WIZ plan %u: bad box %d\n", pl, p->steps[k].var); }
        }
    }
    /* "_" = box left blank (just ENTER) */
    wiz_walk(WP_G2MOL, "36.0,18.015", "n", 1.99833);
    wiz_walk(WP_MOL2G, "0.250,58.44", "m", 14.61);
    wiz_walk(WP_PART,  "3.011E23,18.015", "m", 9.0075);
    wiz_walk(WP_MOLAR, "5.85,58.44,250", "c", 0.400411);
    wiz_walk(WP_MOLAR, "_,58.44,500,0.400", "m", 11.688);
    wiz_walk(WP_DIL,   "6.0,_,1.5,500", "V1", 125);                       /* water fills itself */
    wiz_walk(WP_TIT,   "_,25.0,0.100,31.4,_", "Ca", 0.1256);
    wiz_walk(WP_TIT,   "_,20.0,0.200,30.0,2", "Ca", 0.15);
    wiz_walk(WP_COMB,  "1.00,2.00,300,2.00,_,400", "V2", 1.33333);      /* Celsius boxes fill themselves */
    wiz_walk(WP_COMB,  "1,4,_,_,2,_,_", "V2", 2);                        /* Boyle: no temperatures */
    wiz_walk(WP_COMB,  "_,2.0,_,25,_,_,_,50", "V2", 2.1677);             /* Charles, in Celsius */
    wiz_walk(WP_IDEAL, "_,_,_,10.0,_,25,2.00,_,_", "P", 4.89324);
    wiz_walk(WP_IDEAL, "1,22.4,_,0,_,_", "n", 0.99932);
    wiz_walk(WP_STP,   "_,_,8.0,32.00", "n", 0.25);
    wiz_walk(WP_STO,   "4.00,_,2.016,2,2,18.015", "mB", 35.744);
    wiz_walk(WP_LIM,   "10.0,2.016,2,64.0,31.998,1,2,18.015", "mP", 72.0645);
    wiz_walk(WP_YLD,   "15.0,20.0", "pct", 75);
    wiz_walk(WP_HEAT,  "100,_,20,80", "q", 25104);                      /* q fills itself */
    wiz_walk(WP_HEAT,  "100,_,20,_,5000", "Tf", 31.9503);
    wiz_walk(WP_CAL,   "50,_,100,100,20,22.0", "c1", 0.214564);
    wiz_walk(WP_PH,    "1E-3", "pH", 3);
    wiz_walk(WP_PH,    "_,9.25", "pH", 9.25);
}

static void table_tests(void) {
    static uint8_t seen[9][18];
    uint8_t r, c;
    char name[24];
    for (uint8_t z = 1; z <= NUM_EL; z++) {
        el_pos(z, &r, &c);
        if (r > 8 || c > 17) { fails++; printf("TABLE %s off the grid (%d,%d)\n", el_sym[z], r, c); continue; }
        if (seen[r][c]) { fails++; printf("TABLE %s and %s share a cell\n", el_sym[z], el_sym[seen[r][c]]); }
        seen[r][c] = z;
        if (!el_mass[z] || !el_sym[z][0]) { fails++; printf("TABLE element %d missing\n", z); }
    }
    /* spot checks of the layout: groups and periods */
    struct { uint8_t z, r, c; } pos[] = { {1,0,0}, {2,0,17}, {5,1,12}, {10,1,17}, {11,2,0}, {17,2,16},
        {21,3,2}, {26,3,7}, {36,3,17}, {55,5,0}, {57,7,3}, {71,7,17}, {72,5,3}, {79,5,10}, {86,5,17},
        {89,8,3}, {103,8,17}, {104,6,3}, {118,6,17} };
    for (unsigned i = 0; i < N(pos); i++) {
        el_pos(pos[i].z, &r, &c);
        if (r != pos[i].r || c != pos[i].c) { fails++; printf("TABLE %s at (%d,%d) want (%d,%d)\n", el_sym[pos[i].z], r, c, pos[i].r, pos[i].c); }
    }
    /* molar masses of formulas built with the picker */
    struct { uint8_t z[3]; uint16_t n[3]; const char *name; double m; } mm[] = {
        { {1, 8, 0},  {2, 1, 0},  "H2O",     18.015 },
        { {11, 17, 0}, {1, 1, 0}, "NaCl",    58.44 },
        { {6, 1, 8},  {6, 12, 6}, "C6H12O6", 180.156 },
        { {20, 6, 8}, {1, 1, 3},  "CaCO3",   100.086 },
        { {1, 16, 8}, {2, 1, 4},  "H2SO4",   98.072 },
    };
    for (unsigned i = 0; i < N(mm); i++) {
        f_len = 0;
        for (int k = 0; k < 3 && mm[i].z[k]; k++) f_set(mm[i].z[k], mm[i].n[k]);
        f_name(name, sizeof name);
        if (strcmp(name, mm[i].name) || fabs(f_mass() - mm[i].m) > 0.01) { fails++; printf("FORMULA %s: got %s %.3f want %.3f\n", mm[i].name, name, f_mass(), mm[i].m); }
    }
    {
        struct { double x; const char *s; } fm[] = { {180.156, "180.156"}, {58.44, "58.44"}, {18.0152, "18.015"}, {4.0026, "4.003"}, {98, "98"}, {12.0105, "12.011"} };
        for (unsigned i = 0; i < N(fm); i++) { fmt_mass(fm[i].x, name); if (strcmp(name, fm[i].s)) { fails++; printf("FMT MASS %g -> %s want %s\n", fm[i].x, name, fm[i].s); } }
    }
    f_set(1, 0);                 /* remove H from H2SO4 */
    f_name(name, sizeof name);
    if (strcmp(name, "SO4")) { fails++; printf("FORMULA remove: %s\n", name); }
    f_len = 0;
}

/* every description, unit, symbol and tip line fits on the calculator screen */
static void fit_tests(void) {
    for (unsigned t = 0; t < NUM_TOPICS; t++) {
        const topic_t *T = &topics[t];
        const char *s = T->tips;
        unsigned len = 0;
        if (strlen(T->name) > 30) { fails++; printf("FIT name too long: %s\n", T->name); }
        for (int i = 0; i < T->nvars; i++) {
            if (strlen(T->vars[i].sym) > 6) { fails++; printf("FIT symbol %s in %s\n", T->vars[i].sym, T->name); }
            if (strlen(T->vars[i].unit) > 9) { fails++; printf("FIT unit %s in %s\n", T->vars[i].unit, T->name); }
            if (strlen(T->vars[i].desc) > 38) { fails++; printf("FIT desc (%u) %s\n", (unsigned)strlen(T->vars[i].desc), T->vars[i].desc); }
        }
        for (int r = 0; r < T->nrules; r++)
            if (strlen(T->rules[r].eq) > 38) { fails++; printf("FIT eq (%u) %s\n", (unsigned)strlen(T->rules[r].eq), T->rules[r].eq); }
        for (; *s; s++) {
            if (*s == '\n') { len = 0; continue; }
            if (++len > 38) { fails++; printf("FIT tip line too long in %s\n", T->name); break; }
        }
    }
    {
        unsigned total = 0;
        for (unsigned c = 0; c < NUM_CATS; c++) total += cats[c].count;
        if (total != NUM_TOPICS) { fails++; printf("CATS add up to %u, topics %u\n", total, (unsigned)NUM_TOPICS); }
    }
}

int main(int argc, char **argv) {
    if (argc > 1 && !strcmp(argv[1], "--solve")) return solve_lines();
    if (argc > 1 && !strcmp(argv[1], "--topics")) return list_topics();
    if (argc > 1 && !strcmp(argv[1], "--sheet")) { tp = &topics[0]; build_formula_sheet(); printf("%s\n[%u chars]\n", text, text_len); return 0; }

    /* every rule consistent (all variables known) */
    full("Density", "m=59.0 Vi=20.0 Vf=27.5");
    full("Density", "m=21.6 L=2.0 Vi=10 Vf=18");
    full("Percent error", "meas=2.65 true=2.70");
    full("Temperature units", "C=25");
    full("Grams, moles, particles", "m=36.0 M=18.015 k=2");
    full("Percent by mass", "part=12.0 whole=44.0");
    full("Average atomic mass", "m1=23.985 f1=78.99 m2=24.986 f2=10.00 m3=25.983");
    full("Empirical formula", "g1=40.0 A1=12.011 g2=6.7 A2=1.008 g3=53.3 A3=15.999 Memp=30.03 Mmol=180.16");
    full("Molarity", "m=5.85 M=58.44 mL=250");
    full("Dilution", "C1=6.0 C2=1.5 V2=500");
    full("Ideal gas PV = nRT", "kPa=101.325 Tc=25 m=4.0 M=4.003");
    full("Combined gas law", "P1=1 V1=2 T1c=27 P2=2 T2c=127");
    full("Gas at STP", "m=8.0 M=32.00");
    full("Partial pressures", "Ptot=2.0 n1=1 ntot=4 P3=0.3");
    full("Stoichiometry (A to B)", "mA=4.00 MA=2.016 a=2 b=2 MB=18.015");
    full("Limiting reactant", "mA=10.0 MA=2.016 a=2 mB=64.0 MB=31.998 b=1 p=2 MP=18.015");
    full("Percent yield", "actual=15.0 theo=20.0");
    full("Heat q = m c dT", "m=100 c=4.184 Ti=20 Tf=80");
    full("Calorimetry (mix)", "m1=50 c1=0.385 T1=100 m2=100 c2=4.184 T2=20");
    full("Heat of reaction", "m=32.0 M=16.04 dH=-890 coef=1");
    full("Light and photons", "nm=500");
    full("Hydrogen atom (Bohr)", "ni=3 nf=2");
    full("pH and pOH", "H=1E-3");
    full("Titration", "Va=25.0 Cb=0.100 Vb=31.4 k=1");

    /* textbook-style problems */
    check("Density", "m=25.0 V=10.0 | d=2.50");
    check("Density", "m=59.0 Vi=20.0 Vf=27.5 | V=7.5 d=7.867");                 /* iron, displacement */
    check("Density", "m=21.6 L=2.00 | V=8.00 d=2.70");                          /* aluminum cube */
    check("Density", "d=13.6 V=5.00 | m=68.0");                                 /* mercury */
    check("Percent error", "meas=2.65 true=2.70 | err=1.852");
    check("Temperature units", "C=25 | K=298.15 F=77");
    check("Temperature units", "F=98.6 | C=37 K=310.15");
    check("Temperature units", "K=0 | C=-273.15 F=-459.67");
    check("Grams, moles, particles", "m=36.0 M=18.015 | n=1.9983 N=1.2034E24");
    check("Grams, moles, particles", "m=36.0 M=18.015 k=2 | Nat=2.4068E24");    /* H atoms */
    check("Grams, moles, particles", "n=0.250 M=58.44 | m=14.61");
    check("Grams, moles, particles", "N=3.011E23 M=18.015 | n=0.5 m=9.0075");
    check("Grams, moles, particles", "m=10.0 n=0.250 | M=40.0");
    check("Percent by mass", "part=12.011 whole=44.009 | pct=27.29");             /* C in CO2 */
    check("Percent by mass", "pct=40.0 whole=250 | part=100");
    check("Average atomic mass", "m1=34.969 f1=75.78 m2=36.966 | avg=35.45 f2=24.22"); /* chlorine */
    check("Average atomic mass", "m1=62.93 m2=64.93 avg=63.55 | f1=69 f2=31");      /* copper, find abundance */
    check("Average atomic mass", "m1=23.985 f1=78.99 m2=24.986 f2=10.00 m3=25.983 f3=11.01 | avg=24.305");
    check("Empirical formula", "g1=40.0 A1=12.011 g2=6.7 A2=1.008 g3=53.3 A3=15.999 | r1=1 r2=1.996 r3=1.0004");
    check("Empirical formula", "g1=40.0 A1=12.011 g2=6.7 A2=1.008 g3=53.3 A3=15.999 Memp=30.026 Mmol=180.16 | x=6.0");
    check("Empirical formula", "g1=69.9 A1=55.845 g2=30.1 A2=15.999 | r1=1 r2=1.503");   /* Fe2O3: x2 */
    check("Molarity", "m=5.85 M=58.44 mL=250 | n=0.1001 c=0.4004");
    check("Molarity", "c=0.400 mL=500 M=58.44 | m=11.688");
    check("Molarity", "n=0.75 c=1.5 | V=0.5 mL=500");
    check("Dilution", "C1=6.0 C2=1.5 V2=500 | V1=125 water=375");
    check("Dilution", "C1=2.0 V1=50 water=150 | V2=200 C2=0.5");
    check("Ideal gas PV = nRT", "n=2.00 Tc=25 V=10.0 | P=4.893 T=298.15");
    check("Ideal gas PV = nRT", "n=1.00 P=1.00 Tc=0 | V=22.41");
    check("Ideal gas PV = nRT", "P=1 Tc=0 M=44.01 | d=1.963");                     /* CO2 density at STP */
    check("Ideal gas PV = nRT", "d=1.25 P=1 Tc=0 | M=28.02");                      /* molar mass from density */
    check("Ideal gas PV = nRT", "kPa=202.65 V=5.0 T=300 | P=2 n=0.4062");
    check("Ideal gas PV = nRT", "mmHg=380 | P=0.5 kPa=50.66");
    check("Combined gas law", "P1=1.00 V1=2.00 T1=300 P2=2.00 T2=400 | V2=1.3333");
    check("Combined gas law", "P1=1 V1=4 P2=2 | V2=2");                             /* Boyle */
    check("Combined gas law", "V1=2.0 T1c=25 T2c=50 | V2=2.1677 T1=298.15");        /* Charles */
    check("Combined gas law", "P1=1.5 T1c=25 T2c=100 | P2=1.8773");                 /* Gay-Lussac */
    check("Combined gas law", "V1=3.0 T1=300 V2=6.0 | T2=600 T2c=326.85");
    check("Gas at STP", "n=0.500 | V=11.2");
    check("Gas at STP", "m=8.0 M=32.00 | n=0.25 V=5.6");
    check("Gas at STP", "V=44.8 M=28.02 | n=2 m=56.04");
    check("Partial pressures", "Ptot=760 P2=23.8 | P1=736.2");                      /* collected over water */
    check("Partial pressures", "Ptot=2.0 n1=1 ntot=4 | x1=0.25 P1=0.5");
    check("Partial pressures", "P1=0.4 P2=0.3 P3=0.2 | Ptot=0.9");
    check("Stoichiometry (A to B)", "mA=4.00 MA=2.016 a=2 b=2 MB=18.015 | nB=1.9841 mB=35.744");
    check("Stoichiometry (A to B)", "nA=3.0 a=1 b=2 | nB=6");
    check("Stoichiometry (A to B)", "mA=10.0 MA=100.09 a=1 b=1 | VB=2.238");        /* CaCO3 -> CO2 at STP */
    check("Stoichiometry (A to B)", "mB=50.0 MB=18.015 a=2 b=2 MA=2.016 | mA=5.595");
    check("Limiting reactant", "mA=10.0 MA=2.016 a=2 mB=64.0 MB=31.998 b=1 p=2 MP=18.015 | nP=4.0003 mP=72.06 lim=2 left=1.9355");
    check("Limiting reactant", "mA=5.0 MA=2.016 a=2 mB=100.0 MB=31.998 b=1 p=2 MP=18.015 | lim=1 mP=44.68");
    check("Percent yield", "actual=15.0 theo=20.0 | pct=75");
    check("Percent yield", "pct=80 theo=50 | actual=40");
    check("Heat q = m c dT", "m=100 Ti=20 Tf=80 | q=25104 kJ=25.104 dT=60");          /* c = water by default */
    check("Heat q = m c dT", "q=1000 m=50 dT=20 | c=1.0");
    check("Heat q = m c dT", "m=100 Ti=20 q=5000 | Tf=31.95");
    check("Heat q = m c dT", "kJ=2.5 m=50 c=0.385 Ti=25 | Tf=154.87");
    check("Calorimetry (mix)", "m1=50 T1=100 m2=100 T2=20 Tf=22.0 | c1=0.2146 q=836.8");
    check("Calorimetry (mix)", "m1=50 c1=0.385 T1=100 m2=100 T2=20 | Tf=23.52");
    check("Heat of reaction", "m=32.0 M=16.04 dH=-890 | n=1.995 q=-1775.6");
    check("Heat of reaction", "m=18.0 M=18.015 dH=6.01 | q=6.005");                   /* melting ice */
    check("Heat of reaction", "m=4.032 M=2.016 dH=-572 coef=2 | q=-572");
    check("Light and photons", "nm=500 | nu=5.996E14 E=3.973E-19 Emol=239.25");
    check("Light and photons", "nu=98.1E6 | lam=3.0561");                            /* FM radio */
    check("Light and photons", "Emol=239.25 | nm=500");
    check("Hydrogen atom (Bohr)", "ni=3 nf=2 | E=3.0278E-19 nm=656.08 nu=4.5695E14");
    check("Hydrogen atom (Bohr)", "ni=1 nf=1E9 | E=2.18E-18");                       /* ionization */
    check("Hydrogen atom (Bohr)", "nf=2 nm=656.08 | ni=3");
    check("pH and pOH", "H=1E-3 | pH=3 pOH=11 OH=1E-11");
    check("pH and pOH", "pH=9.25 | H=5.6234E-10 pOH=4.75 OH=1.7783E-5");
    check("pH and pOH", "OH=0.020 | pOH=1.699 pH=12.301");
    check("Titration", "Va=25.0 Cb=0.100 Vb=31.4 | Ca=0.1256 na=0.00314");
    check("Titration", "Va=20.0 Cb=0.200 Vb=30.0 k=2 | Ca=0.15");
    check("Titration", "Ca=0.10 Va=25 Cb=0.20 | Vb=12.5");

    /* a note that must appear */
    check("Combined gas law", "P1=1 V1=4 P2=2 |");
    if (!situation_note() || strcmp(situation_note(), "T blank: held constant (Boyle)")) { fails++; printf("NOTE Boyle missing\n"); }
    check("Empirical formula", "g1=69.9 A1=55.845 g2=30.1 A2=15.999 |");
    if (!situation_note() || strcmp(situation_note(), "multiply every r by 2")) { fails++; printf("NOTE Fe2O3 x2 missing: %s\n", situation_note() ? situation_note() : "none"); }
    check("Limiting reactant", "mA=10.0 MA=2.016 a=2 mB=64.0 MB=31.998 b=1 p=2 MP=18.015 |");
    if (!situation_note() || strcmp(situation_note(), "B runs out first: B is limiting")) { fails++; printf("NOTE limiting missing\n"); }

    table_tests();
    fit_tests();
    wiz_tests();
    printf(fails ? "\n%d FAILURES\n" : "\nALL TESTS PASSED\n", fails);
    return fails != 0;
}

#else   /* calculator build */
/* ================================================================== */

#include <graphx.h>
#include <keypadc.h>
#if USE_CLAUDE
#include <srldrvce.h>
#include <usbdrvce.h>
#else
#include <time.h>
#endif

/* ---------- time ---------------------------------------------------- */
#if USE_CLAUDE
static uint32_t now_ms(void) { return usb_GetCycleCounter() / 48000UL; }
#else
static uint32_t now_ms(void) { return (uint32_t)clock() / 33UL; }   /* 32768 Hz -> about ms */
#endif

/* ---------- layout + colors ---------------------------------------- */
#define SCREEN_W   320
#define COLS       38
#define TEXT_X     8
#define TITLE_H    18
#define BODY_Y     24
#define LINE_H     10
#define BODY_LINES 18
#define ROW_H      15
#define ROWS       10
#define HELP_Y     194
#define SOFT_Y     222
#define MAX_LINES  640

enum { C_BG, C_TEXT, C_ACCENT, C_DIM, C_ONACC, C_RULE, C_SOLVED, C_WARN, C_BAR };
#define C_NONE 254
#define RGB gfx_RGBTo1555
static const uint16_t pal_paper[] = {
    RGB(250, 249, 246), RGB(28, 28, 30), RGB(46, 92, 184), RGB(120, 120, 126),
    RGB(255, 255, 255), RGB(226, 226, 232), RGB(14, 124, 102), RGB(196, 60, 40), RGB(236, 236, 240) };
static const uint16_t pal_night[] = {
    RGB(22, 23, 28), RGB(236, 236, 240), RGB(110, 160, 255), RGB(140, 142, 152),
    RGB(16, 18, 24), RGB(52, 54, 64), RGB(98, 214, 170), RGB(255, 130, 110), RGB(36, 38, 46) };
static bool night = false;
static void apply_theme(void) { gfx_SetPalette(night ? pal_night : pal_paper, sizeof pal_paper, 0); }

/* ---------- app state ---------------------------------------------- */
enum { ST_HOME, ST_TOPICS, ST_FORM, ST_READ, ST_WAIT, ST_COMPOSE, ST_WIZ, ST_PT };
static uint8_t state = ST_HOME, read_back = ST_FORM;
static bool dirty = true;
static const char *notice = NULL;
static const char *read_title = "";
static bool read_is_ai = false;

static uint8_t home_sel, topic_sel, cur_cat, cur_topic, row, row_top;
static bool pending;            /* typed something; re-solve once the keys go quiet */
static uint32_t solve_at;
#if HAVE_WIZ
static const wplan_t *wiz_done;     /* walk-through finished: show its answer */
#endif

static uint16_t line_start[MAX_LINES];
static uint8_t line_len[MAX_LINES];
static uint16_t line_count, scroll;

/* ---------- keyboard ------------------------------------------------ */
static uint8_t prev_keys[8], new_keys[8];
static kb_lkey_t held_key;
static uint32_t held_since, last_repeat;
#define PRESSED(k) ((new_keys[(k) >> 8] & (uint8_t)(k)) != 0)
#define DOWN(k)    ((prev_keys[(k) >> 8] & (uint8_t)(k)) != 0)

static void scan_keys(void) {
    kb_Scan();
    for (uint8_t g = 1; g < 8; g++) {
        uint8_t cur = kb_Data[g];
        new_keys[g] = cur & (uint8_t)~prev_keys[g];
        prev_keys[g] = cur;
    }
}
static bool pressed_repeat(kb_lkey_t k) {
    uint32_t t = now_ms();
    if (PRESSED(k)) { held_key = k; held_since = last_repeat = t; return true; }
    if (held_key == k) {
        if (!DOWN(k)) { held_key = 0; return false; }
        if (t - held_since > 400 && t - last_repeat > 70) { last_repeat = t; return true; }
    }
    return false;
}
static const kb_lkey_t num_keys[10] = { kb_Key0, kb_Key1, kb_Key2, kb_Key3, kb_Key4,
                                        kb_Key5, kb_Key6, kb_Key7, kb_Key8, kb_Key9 };
static int8_t digit_pressed(void) {
    for (int8_t i = 0; i < 10; i++) if (PRESSED(num_keys[i])) return i;
    return -1;
}

/* ---------- text wrapping + reader ---------------------------------- */
static void wrap_text(void) {
    uint16_t i = 0;
    line_count = 0;
    while (i < text_len && line_count < MAX_LINES) {
        uint16_t start = i, last_space = 0xFFFF;
        uint8_t len = 0;
        while (i < text_len && text[i] != '\n' && len < COLS) {
            if (text[i] == ' ') last_space = i;
            i++; len++;
        }
        if (i < text_len) {
            if (text[i] == '\n') i++;
            else if (len == COLS) {
                if (text[i] == ' ') i++;
                else if (last_space != 0xFFFF && last_space > start) {
                    len = (uint8_t)(last_space - start);
                    i = last_space + 1;
                }
            }
        }
        line_start[line_count] = start;
        line_len[line_count] = len;
        line_count++;
    }
}
static uint16_t max_scroll(void) {
    return line_count > BODY_LINES ? (uint16_t)(line_count - BODY_LINES) : 0;
}
static void open_reader(const char *title, uint8_t back, bool ai) {
    wrap_text();
    scroll = 0;
    read_title = title;
    read_back = back;
    read_is_ai = ai;
    state = ST_READ;
    dirty = true;
}

/* ---------- drawing ------------------------------------------------- */
static void print_n(const char *s, uint8_t n, int x, int y) {
    char tmp[COLS + 3];
    if (n > COLS + 2) n = COLS + 2;
    memcpy(tmp, s, n);
    tmp[n] = 0;
    gfx_PrintStringXY(tmp, x, y);
}
static void print_right(const char *s, int right, int y) {
    gfx_PrintStringXY(s, right - (int)gfx_GetStringWidth(s), y);
}

#if USE_CLAUDE
static bool srl_open = false, bridge_ready = false, usb_ok = true;
#endif

/* soft: five labels above the Y= WINDOW ZOOM TRACE GRAPH keys */
static void draw_frame(const char *title, const char *const soft[5]) {
    gfx_FillScreen(C_BG);
    gfx_SetColor(C_ACCENT);
    gfx_FillRectangle_NoClip(0, 0, SCREEN_W, TITLE_H);
    gfx_SetTextFGColor(C_ONACC);
    gfx_PrintStringXY(title, TEXT_X, 5);
#if USE_CLAUDE
    print_right(!srl_open ? "offline" : bridge_ready ? "AI ready" : "no bridge", SCREEN_W - TEXT_X, 5);
#else
    print_right("CHEM", SCREEN_W - TEXT_X, 5);
#endif
    gfx_SetColor(C_BAR);
    gfx_FillRectangle_NoClip(0, SOFT_Y - 4, SCREEN_W, 240 - SOFT_Y + 4);
    gfx_SetColor(C_RULE);
    gfx_HorizLine_NoClip(0, SOFT_Y - 4, SCREEN_W);
    if (notice) {
        gfx_SetTextFGColor(C_WARN);
        gfx_PrintStringXY(notice, TEXT_X, SOFT_Y + 2);
        return;
    }
    for (uint8_t i = 0; i < 5; i++) {
        if (!soft || !soft[i] || !soft[i][0]) continue;
        if (i) gfx_VertLine_NoClip(i * 64, SOFT_Y - 1, 15);
        gfx_SetTextFGColor(C_DIM);
        gfx_PrintStringXY(soft[i], i * 64 + (64 - (int)gfx_GetStringWidth(soft[i])) / 2, SOFT_Y + 2);
    }
}

static void draw_list(const char *title, const char *heading, const char *const *items,
                      uint8_t count, uint8_t sel, const char *const soft[5]) {
    draw_frame(title, soft);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY(heading, TEXT_X, BODY_Y + 2);
    for (uint8_t i = 0; i < count; i++) {
        int y = BODY_Y + 20 + i * (count > 8 ? 19 : count > 7 ? 20 : 22);
        char num[4] = { (char)('1' + i), 0, 0, 0 };
        if (i == sel) {
            gfx_SetColor(C_ACCENT);
            gfx_FillRectangle_NoClip(TEXT_X - 2, (uint8_t)(y - 6), SCREEN_W - 2 * TEXT_X + 4, 18);
            gfx_SetTextFGColor(C_ONACC);
        } else {
            gfx_SetTextFGColor(C_DIM);
        }
        gfx_PrintStringXY(num, TEXT_X + 6, y);
        if (i != sel) gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY(items[i], TEXT_X + 28, y);
        if (i == sel) print_right(">", SCREEN_W - TEXT_X - 6, y);
    }
}

#define APP_TITLE "CHEM Chemistry 1 solver"

static void draw_home(void) {
#if HAVE_WIZ
#define HELP_SOFT "Help"
#else
#define HELP_SOFT ""
#endif
#if USE_CLAUDE
    static const char *const soft[5] = { "Formula", "Ask AI", HELP_SOFT, "Theme", "Quit" };
#else
    static const char *const soft[5] = { "Formula", "Table", HELP_SOFT, "Theme", "Quit" };
#endif
    const char *items[NUM_CATS];
    for (uint8_t i = 0; i < NUM_CATS; i++) items[i] = cats[i].name;
    draw_list(APP_TITLE, "What is the problem about?", items, NUM_CATS, home_sel, soft);
#if HAVE_WIZ
    gfx_SetTextFGColor(C_ACCENT);
    gfx_PrintStringXY("Not sure? ZOOM = Help me pick", TEXT_X, HELP_Y + 13);
#endif
}
static void draw_topics(void) {
    static const char *const soft[5] = { "Back", "", HELP_SOFT, "", "" };
    const char *items[10];
    const cat_t *c = &cats[cur_cat];
    for (uint8_t i = 0; i < c->count; i++) items[i] = topics[cat_first(cur_cat) + i].name;
    draw_list(c->name, "Pick the kind of problem", items, c->count, topic_sel, soft);
}

static void draw_form(void) {
#if USE_CLAUDE
    static const char *const soft[5] = { "Menu", "Reset", "Tips", "Steps", "Ask AI" };
#else
    static const char *const soft[5] = { "Menu", "Reset", "Tips", "Steps", "Table" };
#endif
    char b[24];
    draw_frame(tp->name, soft);

    for (uint8_t r = 0; r < ROWS; r++) {
        uint8_t i = row_top + r;
        int y = BODY_Y + r * ROW_H;
        const var_t *V;
        if (i >= tp->nvars) break;
        V = &tp->vars[i];
        if (i == row) {
            gfx_SetColor(C_RULE);
            gfx_FillRectangle_NoClip(2, (uint8_t)(y - 1), SCREEN_W - 4, ROW_H - 1);
            gfx_SetColor(C_ACCENT);
            gfx_FillRectangle_NoClip(2, (uint8_t)(y - 1), 3, ROW_H - 1);
        }
        gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY(V->sym, TEXT_X + 2, y + 2);

        if (txt[i][0] || (i == row && st[i] == S_UNK)) {
            uint8_t n = (uint8_t)strlen(txt[i]);
            memcpy(b, txt[i], n);
            if (i == row) b[n++] = '_';
            b[n] = 0;
            gfx_SetTextFGColor(pending || st[i] == S_USER || !txt[i][0] ? C_TEXT : C_WARN);
            gfx_PrintStringXY(b, 72, y + 2);
        } else if (st[i] == S_SOLVED || st[i] == S_DEF) {
            fmt_num(val[i], b);
            gfx_SetTextFGColor(st[i] == S_SOLVED ? C_SOLVED : C_DIM);
            gfx_PrintStringXY("=", 60, y + 2);
            gfx_PrintStringXY(b, 72, y + 2);
        } else {
            gfx_SetTextFGColor(C_DIM);
            gfx_PrintStringXY("-", 72, y + 2);
        }
        gfx_SetTextFGColor(C_DIM);
        gfx_PrintStringXY(V->unit, 172, y + 2);
        print_right(st[i] == S_USER ? "given" : st[i] == S_SOLVED ? "solved" : st[i] == S_DEF ? "assumed" : "",
                    SCREEN_W - TEXT_X, y + 2);
    }
    /* more rows above/below */
    gfx_SetTextFGColor(C_DIM);
    if (row_top > 0) gfx_PrintStringXY("more above", TEXT_X + 2, BODY_Y + ROWS * ROW_H + 3);
    if (row_top + ROWS < tp->nvars) gfx_PrintStringXY("more below", SCREEN_W - TEXT_X - 80, BODY_Y + ROWS * ROW_H + 3);

    /* help line for the selected value, or a warning */
    gfx_SetColor(C_RULE);
    gfx_HorizLine_NoClip(TEXT_X, HELP_Y - 2, SCREEN_W - 2 * TEXT_X);
    if (conflict >= 0 && !pending) {
        gfx_SetTextFGColor(C_WARN);
        gfx_PrintStringXY("Numbers disagree. Check:", TEXT_X, HELP_Y + 2);
        print_n(tp->rules[conflict].eq, (uint8_t)strlen(tp->rules[conflict].eq), TEXT_X, HELP_Y + 13);
    }
#if HAVE_WIZ
    else if (guide && gstep < guide->nsteps) {
        /* walk-through: what to type in this box, and what words to look for */
        const wstep_t *S = &guide->steps[gstep];
        char cnt[4] = { (char)('1' + gstep), '/', (char)('0' + guide->nsteps), 0 };
        gfx_SetTextFGColor(C_ACCENT);
        gfx_PrintStringXY(S->ask, TEXT_X, HELP_Y + 2);
        gfx_SetTextFGColor(C_DIM);
        print_right(cnt, SCREEN_W - TEXT_X, HELP_Y + 2);
        gfx_PrintStringXY(S->look, TEXT_X, HELP_Y + 13);
    } else if (wiz_done && !pending && wiz_done->answer == WALL) {
        gfx_SetTextFGColor(C_SOLVED);
        gfx_PrintStringXY(wiz_done->done, TEXT_X, HELP_Y + 2);
        gfx_SetTextFGColor(C_DIM);
        gfx_PrintStringXY("ENTER on a solved value: its steps", TEXT_X, HELP_Y + 13);
    } else if (wiz_done && !pending) {
        uint8_t a = wiz_answer(wiz_done);
        if (st[a] == S_SOLVED) {
            char ans[COLS + 8] = "ANSWER  ";
            strcat(ans, tp->vars[a].sym);
            strcat(ans, " = ");
            fmt_num(val[a], b);
            strcat(ans, b);
            if (tp->vars[a].unit[0] && tp->vars[a].unit[0] != '(') { strcat(ans, " "); strcat(ans, tp->vars[a].unit); }
            gfx_SetTextFGColor(C_SOLVED);
            print_n(ans, (uint8_t)strlen(ans), TEXT_X, HELP_Y + 2);
            gfx_SetTextFGColor(C_DIM);
            gfx_PrintStringXY(wiz_done->done, TEXT_X, HELP_Y + 13);
        } else {
            gfx_SetTextFGColor(C_WARN);
            gfx_PrintStringXY("Not enough numbers to solve it yet.", TEXT_X, HELP_Y + 2);
            gfx_SetTextFGColor(C_DIM);
            gfx_PrintStringXY("UP: go back a box   ZOOM: tips", TEXT_X, HELP_Y + 13);
        }
    }
#endif
    else {
        gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY(tp->vars[row].desc, TEXT_X, HELP_Y + 2);
        /* second line: what this value still needs, or a remark about the result */
        {
            const char *note = pending ? NULL : situation_note();
            uint16_t miss = (!pending && st[row] == S_UNK && !txt[row][0]) ? easiest_missing(row) : 0;
            if (miss) {
                char hint[COLS + 12] = "to get it, enter: ";
                for (uint8_t j = 0, first = 1; j < tp->nvars; j++) {
                    if (!(miss & B(j)) || strlen(hint) + strlen(tp->vars[j].sym) + 2 > COLS) continue;
                    if (!first) strcat(hint, ", ");
                    strcat(hint, tp->vars[j].sym);
                    first = 0;
                }
                gfx_SetTextFGColor(C_ACCENT);
                gfx_PrintStringXY(hint, TEXT_X, HELP_Y + 13);
            } else if (note) {
                gfx_SetTextFGColor(C_WARN);
                gfx_PrintStringXY(note, TEXT_X, HELP_Y + 13);
            } else {
                gfx_SetTextFGColor(C_DIM);
                gfx_PrintStringXY(st[row] == S_SOLVED ? "ENTER: show the steps for this"
                                  : "type the givens; blanks get solved", TEXT_X, HELP_Y + 13);
            }
        }
    }
}

static void draw_read(void) {
    static const char *const soft_ai[5] = { "Back", "", "", "", "Reply" };
    static const char *const soft[5] = { "Back", "", "", "", "" };
    draw_frame(read_title, read_is_ai ? soft_ai : soft);
    for (uint16_t r = 0; r < BODY_LINES; r++) {
        uint16_t i = scroll + r;
        const char *s;
        if (i >= line_count) break;
        s = text + line_start[i];
        gfx_SetTextFGColor(line_len[i] && s[0] != ' ' && !read_is_ai ? C_ACCENT : C_TEXT);
        if (line_len[i] > 3 && s[0] == ' ' && !read_is_ai) gfx_SetTextFGColor(C_SOLVED);
        print_n(s, line_len[i], TEXT_X, BODY_Y + (int)r * LINE_H);
    }
    if (line_count > BODY_LINES) {
        int track = BODY_LINES * LINE_H;
        int bar = track * BODY_LINES / (int)line_count;
        int top = (int)((uint32_t)track * scroll / line_count);
        if (bar < 6) bar = 6;
        gfx_SetColor(C_RULE);
        gfx_FillRectangle_NoClip(SCREEN_W - 4, BODY_Y, 2, (uint8_t)track);
        gfx_SetColor(C_DIM);
        gfx_FillRectangle_NoClip(SCREEN_W - 4, (uint8_t)(BODY_Y + top), 2, (uint8_t)bar);
    }
}

/* ================================================================== */
/* Claude link (USB serial to bridge.js) - hybrid build only           */
/* ================================================================== */
#if USE_CLAUDE
#define MSG_END      '\0'
#define BRIDGE_HELLO 0x06
#define NEW_CHAT     0x02
#define MAX_INPUT    200

static srl_device_t srl;
static uint8_t srl_buf[1024];
static usb_device_t pending_dev = NULL;
static uint32_t last_hello;
static bool discard_reply = false;
static uint16_t reply_offset;
static uint8_t wait_back = ST_FORM;
static char input[MAX_INPUT + 1];
static uint8_t input_len, mode;      /* mode: 0 ABC, 1 abc, 2 123 */
static char wire[900];
static uint16_t wire_len;

static void close_serial(void) {
    if (srl_open) { srl_Close(&srl); srl_open = false; }
    bridge_ready = false;
    dirty = true;
}
static usb_error_t handle_usb_event(usb_event_t event, void *event_data,
                                    usb_callback_data_t *callback_data) {
    usb_error_t err;
    if ((err = srl_UsbEventCallback(event, event_data, callback_data)) != USB_SUCCESS) return err;
    if (event == USB_HOST_CONFIGURE_EVENT && !srl_open) {
        usb_device_t dev = usb_FindDevice(NULL, NULL, USB_SKIP_HUBS);
        if (dev != NULL) pending_dev = dev;
    }
    if (event == USB_DEVICE_DISCONNECTED_EVENT) {
        if (srl_open && (usb_device_t)event_data == srl.dev) close_serial();
        pending_dev = NULL;
    }
    return USB_SUCCESS;
}
static void open_pending_serial(void) {
    srl_error_t e;
    if (pending_dev == NULL || srl_open) return;
    e = srl_Open(&srl, pending_dev, srl_buf, sizeof srl_buf, SRL_INTERFACE_ANY, 115200);
    pending_dev = NULL;
    srl_open = (e == SRL_SUCCESS);
    dirty = true;
}
static bool send_bytes(const char *data, size_t len) {
    uint32_t start = now_ms();
    while (len > 0) {
        int n;
        if (!srl_open) return false;
        n = srl_Write(&srl, data, len);
        if (n < 0) { close_serial(); return false; }
        data += n;
        len -= (size_t)n;
        usb_HandleEvents();
        if (now_ms() - start > 3000) return false;
    }
    return true;
}

/* send `q` to Claude; the reader shows `echo` then the answer */
static void ask_claude(const char *q, uint16_t qlen, const char *echo, bool fresh, uint8_t back) {
    static const char newchat[2] = { NEW_CHAT, '\n' };
    if (!srl_open)     { notice = "Ask AI: plug into the Chromebook"; dirty = true; return; }
    if (!bridge_ready) { notice = "Ask AI: start bridge.js + link page"; dirty = true; return; }
    if (fresh) send_bytes(newchat, 2);
    if (!send_bytes(q, qlen) || !send_bytes("\n", 1)) { notice = "Send failed - check the cable"; dirty = true; return; }
    text_len = 0;
    tput("> "); tput(echo); tput("\n\n");
    reply_offset = text_len;
    wrap_text();
    wait_back = back;
    state = ST_WAIT;
    dirty = true;
}

static void wput(const char *s) {
    while (*s && wire_len < sizeof wire - 1) wire[wire_len++] = *s++;
    wire[wire_len] = 0;
}
static void ask_ai_form(void) {
    char echo[120], b[20];
    bool any = false;
    wire_len = 0;
    wput("Physics 1 problem (Young and Freedman, SI units, g = 9.80). Type: ");
    wput(tp->name);
    wput(". Given: ");
    for (uint8_t i = 0; i < tp->nvars; i++) {
        if (st[i] != S_USER) continue;
        if (any) wput("; ");
        wput(tp->vars[i].desc); wput(" = "); wput(txt[i]); wput(" "); wput(tp->vars[i].unit);
        any = true;
    }
    if (!any) { notice = "Type at least one number first"; dirty = true; return; }
    wput(". Find ");
    if (st[row] == S_USER) wput("every other quantity you can");
    else wput(tp->vars[row].desc);
    if (st[row] == S_SOLVED) {
        fmt_num(val[row], b);
        wput(" (my calculator got "); wput(b); wput(", confirm or correct it)");
    }
    wput(". Reply in at most 9 short lines: EQN lines, a SUB line with numbers, then ANS with units.");
    strcpy(echo, tp->name);
    strcat(echo, ": find ");
    strcat(echo, st[row] == S_USER ? "the rest" : tp->vars[row].sym);
    ask_claude(wire, wire_len, echo, true, ST_FORM);
}

static void poll_serial(void) {
    char buf[64];
    int n;
    if (!srl_open) return;
    n = srl_Read(&srl, buf, sizeof buf);
    if (n < 0) { close_serial(); return; }
    for (int i = 0; i < n; i++) {
        char c = buf[i];
        if (c == BRIDGE_HELLO) {
            if (!bridge_ready) dirty = true;
            bridge_ready = true;
            last_hello = now_ms();
            continue;
        }
        if (discard_reply) { if (c == MSG_END) discard_reply = false; continue; }
        if (state != ST_WAIT) continue;
        if (c == MSG_END) {
            while (text_len > reply_offset && (text[text_len - 1] == '\n' || text[text_len - 1] == ' ')) text_len--;
            open_reader("Claude", wait_back, true);
        } else if (text_len < MAX_TEXT - 1) {
            text[text_len++] = c;
        }
    }
    if (bridge_ready && state != ST_WAIT && now_ms() - last_hello > 8000) { bridge_ready = false; dirty = true; }
}

typedef struct { kb_lkey_t key; char alpha; char num; } keymap_t;
static const keymap_t keymap[] = {
    { kb_KeyMath,   'A',  0  }, { kb_KeyApps,  'B',  0  }, { kb_KeyPrgm,   'C',  0  },
    { kb_KeyRecip,  'D',  0  }, { kb_KeySin,   'E',  0  }, { kb_KeyCos,    'F',  0  },
    { kb_KeyTan,    'G',  0  }, { kb_KeyPower, 'H', '^' }, { kb_KeySquare, 'I',  0  },
    { kb_KeyComma,  'J', ',' }, { kb_KeyLParen,'K', '(' }, { kb_KeyRParen, 'L', ')' },
    { kb_KeyDiv,    'M', '/' }, { kb_KeyLog,   'N',  0  }, { kb_Key7,      'O', '7' },
    { kb_Key8,      'P', '8' }, { kb_Key9,     'Q', '9' }, { kb_KeyMul,    'R', '*' },
    { kb_KeyLn,     'S',  0  }, { kb_Key4,     'T', '4' }, { kb_Key5,      'U', '5' },
    { kb_Key6,      'V', '6' }, { kb_KeySub,   'W', '-' }, { kb_KeySto,    'X', '=' },
    { kb_Key1,      'Y', '1' }, { kb_Key2,     'Z', '2' }, { kb_Key3,      ',', '3' },
    { kb_Key0,      ' ', '0' }, { kb_KeyDecPnt,'.', '.' }, { kb_KeyChs,    '?', '-' },
    { kb_KeyAdd,    '\'', '+' },
};
static char typed_char(void) {
    for (uint8_t i = 0; i < sizeof keymap / sizeof keymap[0]; i++) {
        char c;
        if (!PRESSED(keymap[i].key)) continue;
        if (mode == 2) return keymap[i].num;
        c = keymap[i].alpha;
        if (mode == 1 && c >= 'A' && c <= 'Z') c = (char)(c + 32);
        return c;
    }
    return 0;
}

static void draw_compose(void) {
    static const char *const soft[5] = { "Back", "", "", "", "Send" };
    static const char *mode_names[] = { "[ABC]", "[abc]", "[123]" };
    char line[COLS + 2];
    uint8_t shown, n;
    int y = BODY_Y + 2;
    draw_frame("Ask Claude", soft);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("Type your question:", TEXT_X, y);
    print_right(mode_names[mode], SCREEN_W - TEXT_X, y);
    y += LINE_H + 6;
    gfx_SetTextFGColor(C_TEXT);
    for (shown = 0; shown + COLS <= input_len; shown += COLS) {
        print_n(input + shown, COLS, TEXT_X, y);
        y += LINE_H;
    }
    n = (uint8_t)(input_len - shown);
    memcpy(line, input + shown, n);
    line[n++] = '_';
    line[n] = 0;
    gfx_PrintStringXY(line, TEXT_X, y);
    y += 2 * LINE_H + 4;
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("ALPHA: ABC / abc / 123", TEXT_X, y);
    gfx_PrintStringXY("0 = space   (-) = ?   ENTER = send", TEXT_X, y + LINE_H + 2);
    gfx_PrintStringXY("Short is fine: WHY, SHOW STEPS", TEXT_X, y + 2 * LINE_H + 4);
}
static void draw_wait(uint8_t dots) {
    static const char *const soft[5] = { "Cancel", "", "", "", "" };
    static const char *spin[] = { "Thinking", "Thinking.", "Thinking..", "Thinking..." };
    draw_frame("Claude", soft);
    gfx_SetTextFGColor(C_ACCENT);
    for (uint16_t i = 0; i < line_count && i < 6; i++)
        print_n(text + line_start[i], line_len[i], TEXT_X, BODY_Y + (int)i * LINE_H);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY(spin[dots & 3], TEXT_X, BODY_Y + 8 * LINE_H);
}
#endif /* USE_CLAUDE */

/* ================================================================== */
/* Screens: key handling                                               */
/* ================================================================== */

#if HAVE_WIZ
static uint8_t wq, wsel, wdepth, wstack[8];
static void show_row(void) {
    if (row < row_top) row_top = row;
    if (row >= row_top + ROWS) row_top = (uint8_t)(row - ROWS + 1);
}
static void wiz_open(void) { wq = 0; wsel = 0; wdepth = 0; state = ST_WIZ; dirty = true; }
static void wiz_finish(void) {
    wiz_done = guide;
    guide = NULL;
    row = wiz_answer(wiz_done);
    show_row();
    dirty = true;
}
static void wiz_go_plan(const wplan_t *p) {
    int t = wiz_topic(p);
    if (t < 0 || !wiz_begin(p)) return;
    cur_topic = (uint8_t)t;
    for (uint8_t c = 0; c < NUM_CATS; c++) {        /* so Back lands in the right list */
        uint8_t f = cat_first(c);
        if (t >= f && t < f + cats[c].count) { cur_cat = home_sel = c; topic_sel = (uint8_t)(t - f); }
    }
    pending = false;
    wiz_done = NULL;
    row_top = 0;
    state = ST_FORM;
    if (gstep >= guide->nsteps) wiz_finish();
    else { row = guide->steps[gstep].var; show_row(); }
    dirty = true;
}
static void draw_wiz(void) {
    static const char *const soft[5] = { "Back", "", "", "", "" };
    const wq_t *q = &wqs[wq];
    draw_frame("Help me pick", soft);
    gfx_SetTextFGColor(C_TEXT);
    gfx_PrintStringXY(q->q, TEXT_X, BODY_Y + 8);
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY(q->hint, TEXT_X, BODY_Y + 22);
    for (uint8_t i = 0; i < q->nopt; i++) {
        int y = BODY_Y + 52 + i * 26;
        char num[2] = { (char)('1' + i), 0 };
        if (i == wsel) {
            gfx_SetColor(C_ACCENT);
            gfx_FillRectangle_NoClip(TEXT_X - 2, (uint8_t)(y - 7), SCREEN_W - 2 * TEXT_X + 4, 20);
            gfx_SetTextFGColor(C_ONACC);
        } else gfx_SetTextFGColor(C_DIM);
        gfx_PrintStringXY(num, TEXT_X + 6, y);
        if (i != wsel) gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY(q->opt[i].label, TEXT_X + 28, y);
    }
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("press the number, or arrows + ENTER", TEXT_X, HELP_Y + 2);
    if (wdepth) gfx_PrintStringXY("CLEAR: back one question", TEXT_X, HELP_Y + 13);
}
static void wiz_keys(void) {
    const wq_t *q = &wqs[wq];
    int8_t d = digit_pressed(), pick = -1, go;
    if (pressed_repeat(kb_KeyDown)) { wsel = (uint8_t)((wsel + 1) % q->nopt); dirty = true; }
    else if (pressed_repeat(kb_KeyUp)) { wsel = (uint8_t)((wsel + q->nopt - 1) % q->nopt); dirty = true; }
    else if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyRight)) pick = (int8_t)wsel;
    else if (d >= 1 && d <= (int8_t)q->nopt) pick = (int8_t)(d - 1);
    else if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyLeft) || PRESSED(kb_KeyYequ)) {
        if (wdepth) { wq = wstack[--wdepth]; wsel = 0; } else state = ST_HOME;
        dirty = true;
    }
    if (pick < 0) return;
    go = q->opt[pick].go;
    if (go >= 0) {
        if (wdepth < sizeof wstack) wstack[wdepth++] = wq;
        wq = (uint8_t)go;
        wsel = 0;
        dirty = true;
    } else wiz_go_plan(&wplans[-1 - go]);
}
#endif

static void open_topic(uint8_t t) {
#if HAVE_WIZ
    guide = NULL;
    wiz_done = NULL;
#endif
    cur_topic = t;
    tp = &topics[t];
    memset(txt, 0, sizeof txt);
    row = row_top = 0;
    pending = false;
    solve();
    state = ST_FORM;
    dirty = true;
}
static void edited(void) {
    pending = true; solve_at = now_ms() + 250; dirty = true;
#if HAVE_WIZ
    wiz_done = NULL;
#endif
}
static void settle(void) { if (pending) { solve(); pending = false; dirty = true; } }

static void show_steps(void) {
    settle();
    build_steps(row);
    open_reader("Steps", ST_FORM, false);
}
/* ================================================================== */
/* Periodic table: build a formula, get its molar mass                 */
/* ================================================================== */
#define PT_X 7
#define PT_W 17
#define PT_H 15
static uint8_t pt_grid[9][18];      /* element number in each cell, 0 = empty */
static uint8_t pt_r, pt_c;          /* cursor */
static uint8_t pt_back;             /* screen to return to */
static int8_t pt_fill = -1;         /* form box that gets the molar mass, -1 = none */
static bool pt_alpha;               /* ALPHA pressed: next key is a letter */
static bool pt_typing;              /* digits are setting the count of this element */
static uint16_t pt_count;

static const struct { kb_lkey_t key; char ch; } pt_letters[] = {
    { kb_KeyMath, 'A' }, { kb_KeyApps, 'B' }, { kb_KeyPrgm, 'C' }, { kb_KeyRecip, 'D' },
    { kb_KeySin, 'E' }, { kb_KeyCos, 'F' }, { kb_KeyTan, 'G' }, { kb_KeyPower, 'H' },
    { kb_KeySquare, 'I' }, { kb_KeyComma, 'J' }, { kb_KeyLParen, 'K' }, { kb_KeyRParen, 'L' },
    { kb_KeyDiv, 'M' }, { kb_KeyLog, 'N' }, { kb_Key7, 'O' }, { kb_Key8, 'P' }, { kb_Key9, 'Q' },
    { kb_KeyMul, 'R' }, { kb_KeyLn, 'S' }, { kb_Key4, 'T' }, { kb_Key5, 'U' }, { kb_Key6, 'V' },
    { kb_KeySub, 'W' }, { kb_KeySto, 'X' }, { kb_Key1, 'Y' }, { kb_Key2, 'Z' },
};

static uint8_t pt_z(void) { return pt_grid[pt_r][pt_c]; }
static void pt_open(uint8_t back, int8_t fill) {
    if (!pt_grid[0][0]) {
        uint8_t r, c;
        for (uint8_t z = 1; z <= NUM_EL; z++) { el_pos(z, &r, &c); pt_grid[r][c] = z; }
    }
    pt_back = back;
    pt_fill = fill;
    pt_alpha = pt_typing = false;
    state = ST_PT;
    dirty = true;
}
static void pt_go(uint8_t z) { el_pos(z, &pt_r, &pt_c); pt_typing = false; dirty = true; }
/* move one step; skip the empty cells */
static void pt_move(int8_t dr, int8_t dc) {
    int8_t r = (int8_t)pt_r, c = (int8_t)pt_c;
    if (dc) {
        for (c = (int8_t)(c + dc); c >= 0 && c < 18; c = (int8_t)(c + dc))
            if (pt_grid[r][c]) { pt_c = (uint8_t)c; break; }
    } else {
        for (r = (int8_t)(r + dr); r >= 0 && r < 9; r = (int8_t)(r + dr)) {
            /* the nearest filled cell in that row */
            for (int8_t k = 0; k < 18; k++) {
                if (c - k >= 0 && pt_grid[r][c - k]) { pt_r = (uint8_t)r; pt_c = (uint8_t)(c - k); goto moved; }
                if (c + k < 18 && pt_grid[r][c + k]) { pt_r = (uint8_t)r; pt_c = (uint8_t)(c + k); goto moved; }
            }
        }
    }
moved:
    pt_typing = false;
    dirty = true;
}
/* ALPHA + letter: next element whose symbol starts with that letter */
static void pt_letter(char ch) {
    uint8_t z0 = pt_z();
    for (uint8_t k = 1; k <= NUM_EL; k++) {
        uint8_t z = (uint8_t)((z0 + k - 1) % NUM_EL + 1);
        if (el_sym[z][0] == ch) { pt_go(z); return; }
    }
}
static void pt_breakdown(void) {
    double M = f_mass();
    char b[24];
    text_len = 0;
    if (!f_len) { tput("Nothing in the formula yet.\n\nMove to an element and press ENTER\n(once per atom), or type how many.\n"); return; }
    f_name(b, sizeof b);
    tput("FORMULA  "); tput(b); tput("\n\n");
    for (uint8_t i = 0; i < f_len; i++) {
        double part = (double)el_mass[f_z[i]] * f_n[i];
        tput(el_sym[f_z[i]]); tput(": ");
        tnum(f_n[i]); tput(" x "); fmt_mass(el_mass[f_z[i]], b); tput(b); tput(" = "); fmt_mass(part, b); tput(b);
        tput("\n   "); tnum(100 * part / M); tput(" % by mass\n");
    }
    fmt_mass(M, b);
    tput("\nMOLAR MASS\n   M = "); tput(b); tput(" g/mol\n");
}
static void draw_pt(void) {
    const char *soft[5] = { "Back", "Clear", "", "Steps", pt_fill >= 0 ? "Use M" : "" };
    char b[40];
    uint8_t zs = pt_z();
    draw_frame("Periodic table: molar mass", soft);
    for (uint8_t r = 0; r < 9; r++)
        for (uint8_t c = 0; c < 18; c++) {
            uint8_t z = pt_grid[r][c];
            int x = PT_X + c * PT_W, y = BODY_Y + 2 + r * PT_H + (r >= 7 ? 5 : 0);
            const char *s;
            if (!z) continue;
            s = el_sym[z];
            if (z == zs) {
                gfx_SetColor(C_ACCENT);
                gfx_FillRectangle_NoClip(x, (uint8_t)y, PT_W - 1, PT_H - 1);
                gfx_SetTextFGColor(C_ONACC);
            } else gfx_SetTextFGColor(f_find(z) >= 0 ? C_SOLVED : z > 103 ? C_DIM : C_TEXT);
            gfx_PrintStringXY(s, x + (s[1] ? 0 : 4), y + 3);
        }
    /* the selected element */
    {
        int8_t i = f_find(zs);
        uint8_t n = 0;
        const char *s = el_sym[zs];
        b[n++] = (char)('0' + zs / 100 % 10); b[n++] = (char)('0' + zs / 10 % 10); b[n++] = (char)('0' + zs % 10);
        while (b[0] == '0' && n > 1) { memmove(b, b + 1, --n); }
        b[n++] = ' '; while (*s) b[n++] = *s++;
        b[n] = 0;
        gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY(b, TEXT_X, BODY_Y + 2 + 9 * PT_H + 8);
        fmt_mass(el_mass[zs], b); strcat(b, " g/mol");
        gfx_PrintStringXY(b, TEXT_X + 64, BODY_Y + 2 + 9 * PT_H + 8);
        if (i >= 0) {
            strcpy(b, "x"); fmt_num(f_n[i], b + 1);
            gfx_SetTextFGColor(C_SOLVED);
            print_right(b, SCREEN_W - TEXT_X, BODY_Y + 2 + 9 * PT_H + 8);
        }
        if (pt_alpha) { gfx_SetTextFGColor(C_ACCENT); print_right("[A] letter?", SCREEN_W - TEXT_X - 40, BODY_Y + 2 + 9 * PT_H + 8); }
    }
    /* the formula so far */
    if (f_len) {
        char m[20];
        f_name(b, 24);
        gfx_SetTextFGColor(C_SOLVED);
        gfx_PrintStringXY(b, TEXT_X, HELP_Y + 2);
        strcpy(m, "M = "); fmt_mass(f_mass(), m + 4);
        print_right(m, SCREEN_W - TEXT_X, HELP_Y + 2);
    } else {
        gfx_SetTextFGColor(C_DIM);
        gfx_PrintStringXY("ENTER adds 1 atom; digits: how many", TEXT_X, HELP_Y + 2);
    }
    gfx_SetTextFGColor(C_DIM);
    gfx_PrintStringXY("ALPHA+letter jumps   DEL removes", TEXT_X, HELP_Y + 13);
}
static void pt_use(void) {
    if (pt_fill >= 0 && f_len) {
        fmt_mass(f_mass(), txt[pt_fill]);
        row = (uint8_t)pt_fill;
        if (row < row_top) row_top = row;
        if (row >= row_top + ROWS) row_top = (uint8_t)(row - ROWS + 1);
        edited();
    }
    state = pt_back;
    dirty = true;
}
static void pt_keys(void) {
    int8_t d = digit_pressed();
    uint8_t z = pt_z();
    if (pt_alpha) {
        for (uint8_t i = 0; i < sizeof pt_letters / sizeof pt_letters[0]; i++)
            if (PRESSED(pt_letters[i].key)) { pt_alpha = false; pt_letter(pt_letters[i].ch); return; }
        if (PRESSED(kb_KeyAlpha) || PRESSED(kb_KeyClear)) { pt_alpha = false; dirty = true; }
        return;
    }
    if (PRESSED(kb_KeyAlpha)) { pt_alpha = true; dirty = true; }
    else if (pressed_repeat(kb_KeyLeft)) pt_move(0, -1);
    else if (pressed_repeat(kb_KeyRight)) pt_move(0, 1);
    else if (pressed_repeat(kb_KeyUp)) pt_move(-1, 0);
    else if (pressed_repeat(kb_KeyDown)) pt_move(1, 0);
    else if (PRESSED(kb_KeyEnter)) {
        int8_t i = f_find(z);
        f_set(z, (uint16_t)(i >= 0 ? f_n[i] + 1 : 1));
        pt_typing = false;
        dirty = true;
    }
    else if (d >= 0) {                          /* type how many of this element */
        pt_count = pt_typing ? (uint16_t)(pt_count * 10 + d) : (uint16_t)d;
        if (pt_count > 999) pt_count = (uint16_t)d;
        f_set(z, pt_count);
        pt_typing = true;
        dirty = true;
    }
    else if (PRESSED(kb_KeyDel)) {
        if (f_find(z) >= 0) f_set(z, 0); else if (f_len) f_len--;
        pt_typing = false;
        dirty = true;
    }
    else if (PRESSED(kb_KeyWindow)) { f_len = 0; pt_typing = false; dirty = true; }
    else if (PRESSED(kb_KeyClear)) {
        if (f_len) f_len = 0; else state = pt_back;
        pt_typing = false;
        dirty = true;
    }
    else if (PRESSED(kb_KeyTrace)) { pt_breakdown(); open_reader("Molar mass", ST_PT, false); }
    else if (PRESSED(kb_KeyGraph) || PRESSED(kb_Key2nd)) pt_use();
    else if (PRESSED(kb_KeyYequ)) { state = pt_back; dirty = true; }
}
/* GRAPH in a topic: the box that wants a molar mass */
static int8_t pt_target(void) {
    if (!strcmp(tp->vars[row].unit, "g/mol")) return (int8_t)row;
    for (uint8_t i = 0; i < tp->nvars; i++) if (!strcmp(tp->vars[i].unit, "g/mol") && st[i] != S_USER) return (int8_t)i;
    for (uint8_t i = 0; i < tp->nvars; i++) if (!strcmp(tp->vars[i].unit, "g/mol")) return (int8_t)i;
    return -1;
}

static void show_tips(void) {
    text_len = 0;
    tput("WORDS IN THE PROBLEM -> WHAT TO TYPE\n\n");
    tput(tp->tips);
    tput("\n\nFORMULAS USED HERE\n");
    put_formulas(tp);
    tput("\nKEYS\n(-) or LEFT: negative  RIGHT: positive\n^ types E (6.02E23)  CLEAR empties a box\nENTER on a solved value = steps\nGRAPH: periodic table, fills molar mass\n");
    open_reader("Tips", ST_FORM, false);
}

static void form_keys(void) {
    int8_t d;
    char c = 0;
    uint8_t n = (uint8_t)strlen(txt[row]);

#if HAVE_WIZ
    if (guide) {
        /* walk-through: ENTER or DOWN = this box is done (blank is fine), UP = back one */
        if (PRESSED(kb_KeyEnter) || pressed_repeat(kb_KeyDown)) {
            settle();
            if (wiz_next()) wiz_finish();
            else { row = guide->steps[gstep].var; show_row(); dirty = true; }
            return;
        }
        if (pressed_repeat(kb_KeyUp)) {
            if (gstep) { gstep--; row = guide->steps[gstep].var; show_row(); dirty = true; }
            return;
        }
        /* any menu key ends the walk-through and then does its usual job */
        if (PRESSED(kb_KeyYequ) || PRESSED(kb_KeyWindow) || (PRESSED(kb_KeyClear) && !n)) { guide = NULL; dirty = true; }
    }
    if (wiz_done && (PRESSED(kb_KeyUp) || PRESSED(kb_KeyDown) || PRESSED(kb_KeyYequ) || PRESSED(kb_KeyWindow) || PRESSED(kb_KeyClear)))
        { wiz_done = NULL; dirty = true; }
#endif
    if (pressed_repeat(kb_KeyDown)) { row = (uint8_t)((row + 1) % tp->nvars); dirty = true; }
    else if (pressed_repeat(kb_KeyUp)) { row = (uint8_t)((row + tp->nvars - 1) % tp->nvars); dirty = true; }
    else if (PRESSED(kb_KeyEnter)) {
        settle();
        if (st[row] == S_SOLVED) show_steps();
        else { row = (uint8_t)((row + 1) % tp->nvars); dirty = true; }
    }
    else if (PRESSED(kb_KeyTrace)) show_steps();
    else if (PRESSED(kb_KeyZoom)) show_tips();
    else if (PRESSED(kb_KeyWindow)) { memset(txt, 0, sizeof txt); pending = false; solve(); dirty = true; }
    else if (PRESSED(kb_KeyYequ)) { state = ST_TOPICS; dirty = true; }
#if USE_CLAUDE
    else if (PRESSED(kb_KeyGraph)) { settle(); ask_ai_form(); }
#else
    else if (PRESSED(kb_KeyGraph)) { settle(); pt_open(ST_FORM, pt_target()); }
#endif
    else if (PRESSED(kb_KeyClear)) {
        if (n) { txt[row][0] = 0; edited(); }
        else state = ST_TOPICS;
        dirty = true;
    }
    else if (pressed_repeat(kb_KeyDel)) {
        if (n) { txt[row][n - 1] = 0; edited(); }
    }
    /* LEFT makes the value negative, RIGHT positive, any time while typing:
       read "0.118 m/s to the left", type .118, press LEFT */
    else if (PRESSED(kb_KeyLeft)) {
        if (txt[row][0] != '-' && n < VAL_LEN) {
            memmove(txt[row] + 1, txt[row], (size_t)n + 1);
            txt[row][0] = '-';
            edited();
        }
    }
    else if (PRESSED(kb_KeyRight)) {
        if (txt[row][0] == '-') { memmove(txt[row], txt[row] + 1, n); edited(); }
    }
    else {
        if ((d = digit_pressed()) >= 0) c = (char)('0' + d);
        else if (PRESSED(kb_KeyDecPnt)) c = '.';
        else if (PRESSED(kb_KeyChs) || PRESSED(kb_KeySub)) c = '-';
        else if (PRESSED(kb_KeyPower) || PRESSED(kb_KeyComma) || PRESSED(kb_Key2nd)) c = 'E';
        if (!c) return;
        /* keep the text a valid number as it is typed */
        if (c == '-' && !(n == 0 || txt[row][n - 1] == 'E')) return;
        if (c == 'E' && (n == 0 || strchr(txt[row], 'E') || txt[row][n - 1] == '-')) return;
        if (c == '.' && (strchr(txt[row], '.') || strchr(txt[row], 'E'))) return;
        if (n >= VAL_LEN) return;
        txt[row][n] = c;
        txt[row][n + 1] = 0;
        edited();
    }
    if (row < row_top) row_top = row;
    if (row >= row_top + ROWS) row_top = (uint8_t)(row - ROWS + 1);
}

static bool home_keys(void) {       /* true = quit */
    int8_t d = digit_pressed(), pick = -1;
    if (pressed_repeat(kb_KeyDown)) { home_sel = (uint8_t)((home_sel + 1) % NUM_CATS); dirty = true; }
    if (pressed_repeat(kb_KeyUp))   { home_sel = (uint8_t)((home_sel + NUM_CATS - 1) % NUM_CATS); dirty = true; }
    if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyRight)) pick = (int8_t)home_sel;
    if (d >= 1 && d <= (int8_t)NUM_CATS) pick = (int8_t)(d - 1);
    if (PRESSED(kb_KeyTrace) || PRESSED(kb_KeyMode)) { night = !night; apply_theme(); dirty = true; }
    if (PRESSED(kb_KeyYequ)) { build_formula_sheet(); open_reader("Formula sheet", ST_HOME, false); return false; }
#if USE_CLAUDE
    if (PRESSED(kb_KeyWindow)) { read_back = ST_HOME; state = ST_COMPOSE; dirty = true; return false; }
#else
    if (PRESSED(kb_KeyWindow)) { pt_open(ST_HOME, -1); return false; }
#endif
    if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyGraph)) return true;
#if HAVE_WIZ
    if (PRESSED(kb_KeyZoom)) { wiz_open(); return false; }
#endif
    if (pick < 0) return false;
    home_sel = cur_cat = (uint8_t)pick;
    topic_sel = 0;
    state = ST_TOPICS;
    dirty = true;
    return false;
}
static void topics_keys(void) {
    const cat_t *c = &cats[cur_cat];
    int8_t d = digit_pressed();
    if (pressed_repeat(kb_KeyDown)) { topic_sel = (uint8_t)((topic_sel + 1) % c->count); dirty = true; }
    else if (pressed_repeat(kb_KeyUp)) { topic_sel = (uint8_t)((topic_sel + c->count - 1) % c->count); dirty = true; }
    else if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyRight)) open_topic((uint8_t)(cat_first(cur_cat) + topic_sel));
    else if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyYequ) || PRESSED(kb_KeyLeft)) { state = ST_HOME; dirty = true; }
    else if (d >= 1 && d <= (int8_t)c->count) { topic_sel = (uint8_t)(d - 1); open_topic((uint8_t)(cat_first(cur_cat) + topic_sel)); }
#if HAVE_WIZ
    else if (PRESSED(kb_KeyZoom)) wiz_open();
#endif
}
static void read_keys(void) {
    uint16_t ms = max_scroll();
    if (pressed_repeat(kb_KeyDown) && scroll < ms) { scroll++; dirty = true; }
    else if (pressed_repeat(kb_KeyUp) && scroll > 0) { scroll--; dirty = true; }
    else if (PRESSED(kb_KeyRight)) { scroll = (uint16_t)(scroll + BODY_LINES - 1); if (scroll > ms) scroll = ms; dirty = true; }
    else if (PRESSED(kb_KeyLeft)) { scroll = scroll > BODY_LINES - 1 ? (uint16_t)(scroll - (BODY_LINES - 1)) : 0; dirty = true; }
    else if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyYequ)) { state = read_back; dirty = true; }
#if USE_CLAUDE
    else if (read_is_ai && (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyGraph))) { state = ST_COMPOSE; dirty = true; }
#endif
    else if (PRESSED(kb_KeyEnter)) { state = read_back; dirty = true; }
}

/* ================================================================== */
int main(void) {
    bool quit = false;
#if USE_CLAUDE
    uint32_t last_anim = 0;
    uint8_t dots = 0;
#endif

    gfx_Begin();
    apply_theme();
    gfx_SetDrawBuffer();
    gfx_SetMonospaceFont(8);
    gfx_SetTextBGColor(C_NONE);
    gfx_SetTextTransparentColor(C_NONE);
    tp = &topics[0];

#if USE_CLAUDE
    /* if USB will not start, the local solver still works */
    if (usb_Init(handle_usb_event, NULL, srl_GetCDCStandardDescriptors(), USB_DEFAULT_INIT_FLAGS) != USB_SUCCESS)
    { usb_ok = false; notice = "USB did not start: solver only"; }
#endif

    while (!quit) {
#if USE_CLAUDE
        if (usb_ok) {
            usb_HandleEvents();
            open_pending_serial();
            poll_serial();
        }
#endif
        scan_keys();
        for (uint8_t g = 1; g < 8; g++)
            if (new_keys[g] && notice) { notice = NULL; dirty = true; break; }
        if (PRESSED(kb_KeyMode) && state != ST_HOME) { night = !night; apply_theme(); dirty = true; }

        switch (state) {
        case ST_HOME:   quit = home_keys(); break;
        case ST_TOPICS: topics_keys(); break;
        case ST_FORM:
            form_keys();
            if (state == ST_FORM && pending && (int32_t)(now_ms() - solve_at) >= 0) settle();
            break;
        case ST_READ:   read_keys(); break;
#if HAVE_WIZ
        case ST_WIZ:    wiz_keys(); break;
#endif
        case ST_PT:     pt_keys(); break;
#if USE_CLAUDE
        case ST_WAIT:
            if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyYequ)) {
                discard_reply = true;
                state = wait_back;
                notice = "Cancelled";
                dirty = true;
            } else if (now_ms() - last_anim > 350) { last_anim = now_ms(); dots++; dirty = true; }
            break;
        case ST_COMPOSE: {
            char c;
            if (PRESSED(kb_KeyAlpha)) { mode = (uint8_t)((mode + 1) % 3); dirty = true; }
            else if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyGraph)) {
                if (input_len) {
                    input[input_len] = 0;
                    ask_claude(input, input_len, input, false, read_back);
                    if (state == ST_WAIT) input_len = 0;
                }
            }
            else if (pressed_repeat(kb_KeyDel)) { if (input_len) { input_len--; dirty = true; } }
            else if (PRESSED(kb_KeyClear)) {
                if (input_len) input_len = 0; else state = read_back;
                dirty = true;
            }
            else if (PRESSED(kb_KeyYequ)) { state = read_back; dirty = true; }
            else if ((c = typed_char()) != 0 && input_len < MAX_INPUT) { input[input_len++] = c; dirty = true; }
            break;
        }
#endif
        default: state = ST_HOME; break;
        }

        if (dirty) {
            dirty = false;
            switch (state) {
            case ST_HOME:   draw_home(); break;
            case ST_TOPICS: draw_topics(); break;
            case ST_FORM:   draw_form(); break;
            case ST_READ:   draw_read(); break;
#if HAVE_WIZ
            case ST_WIZ:    draw_wiz(); break;
#endif
            case ST_PT:     draw_pt(); break;
#if USE_CLAUDE
            case ST_WAIT:   draw_wait(dots); break;
            case ST_COMPOSE: draw_compose(); break;
#endif
            default: break;
            }
            gfx_SwapDraw();
        }
    }

#if USE_CLAUDE
    close_serial();
    usb_Cleanup();
#endif
    gfx_End();
    return 0;
}
#endif /* HOST_TEST */
