/*
 * PHYS1 - Physics 1 problem solver for the TI-84 Plus CE
 *
 * Pick a topic, type the numbers a word problem gives you, and every other
 * quantity fills in by itself. ENTER on a solved value shows the steps.
 *
 * Built from this one file (48 problem types, too big for one program):
 *   PART=1  PHYS1  motion, forces, energy
 *   PART=2  PHYS2  momentum, collisions, changing forces
 *   PART=3  PHYS3  rotation, gravity, springs, fluids
 *   USE_CLAUDE=1 adds the "Ask AI" key over USB (needs bridge.js);
 *   USE_CLAUDE=0 builds the same solver with no USB code (PHL1-3).
 *
 * Keys (form): digits . (-) type a value   ^ or , = E (3E8)
 *   UP/DOWN field   DEL backspace   CLEAR erase field / go back
 *   ENTER steps for this value      MODE light/dark
 *   Y= menu  WINDOW reset  ZOOM tips  TRACE steps  GRAPH Ask AI
 */

#ifndef USE_CLAUDE
#define USE_CLAUDE 1
#endif
#ifndef PART
#define PART 0      /* 0 = every topic (PC tests); 1, 2, 3 = the three calculator programs */
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

/* ---------- 1. straight-line motion (also rotational kinematics) ---- */
/* 0=v0 1=v 2=a 3=t 4=dx */
FN(k_v,     v[0] + v[2] * v[3])
FN(k_v0,    v[1] - v[2] * v[3])
FN(k_a,     (v[1] - v[0]) / v[3])
FN(k_t,     (v[1] - v[0]) / v[2])
FN(k_dxav,  0.5 * (v[0] + v[1]) * v[3])
FN(k_tav,   2 * v[4] / (v[0] + v[1]))
FN(k_v0av,  2 * v[4] / v[3] - v[1])
FN(k_vav,   2 * v[4] / v[3] - v[0])
FN(k_dx,    v[0] * v[3] + 0.5 * v[2] * SQ(v[3]))
FN(k_v0dx,  (v[4] - 0.5 * v[2] * SQ(v[3])) / v[3])
FN(k_adx,   2 * (v[4] - v[0] * v[3]) / SQ(v[3]))
FN(k_dxv2,  (SQ(v[1]) - SQ(v[0])) / (2 * v[2]))
FN(k_av2,   (SQ(v[1]) - SQ(v[0])) / (2 * v[4]))
FN(k_dxvt,  v[1] * v[3] - 0.5 * v[2] * SQ(v[3]))
FN(k_vvt,   (v[4] + 0.5 * v[2] * SQ(v[3])) / v[3])
FN(k_avt,   2 * (v[1] * v[3] - v[4]) / SQ(v[3]))
static double k_tq(const double *v) {
    double a = v[2], b = v[0], dx = v[4], s, t1, t2;
    if (fabs(a) < 1e-9) return dx / b;
    s = rootz(b * b + 2 * a * dx, b * b + fabs(2 * a * dx));
    if (!ok(s)) return UNK;
    t1 = (-b + s) / a; t2 = (-b - s) / a;
    if (t1 > t2) { double x = t1; t1 = t2; t2 = x; }
    if (flip & B(3)) return t1 > 1e-9 ? t2 : t1;   /* the other time */
    return t1 > 1e-9 ? t1 : t2;           /* first time it gets there */
}
static double k_vv2(const double *v) {
    double s = rootz(SQ(v[0]) + 2 * v[2] * v[4], SQ(v[0]) + fabs(2 * v[2] * v[4]));
    return v[4] < 0 ? -s : s;
}
static double k_v0v2(const double *v) {
    double s = rootz(SQ(v[1]) - 2 * v[2] * v[4], SQ(v[1]) + fabs(2 * v[2] * v[4]));
    return v[4] < 0 ? -s : s;
}

#define KIN_RULES(E1, E2, E3, E4, E5) { \
    {1, B(0)|B(2)|B(3), k_v,    E1, 0}, {0, B(1)|B(2)|B(3), k_v0,   E1, 0}, \
    {2, B(0)|B(1)|B(3), k_a,    E1, 0}, {3, B(0)|B(1)|B(2), k_t,    E1, 0}, \
    {4, B(0)|B(1)|B(3), k_dxav, E4, 0}, {3, B(0)|B(1)|B(4), k_tav,  E4, 0}, \
    {0, B(1)|B(3)|B(4), k_v0av, E4, 0}, {1, B(0)|B(3)|B(4), k_vav,  E4, 0}, \
    {4, B(0)|B(2)|B(3), k_dx,   E2, 0}, {0, B(2)|B(3)|B(4), k_v0dx, E2, 0}, \
    {2, B(0)|B(3)|B(4), k_adx,  E2, 0}, {4, B(0)|B(1)|B(2), k_dxv2, E3, 0}, \
    {2, B(0)|B(1)|B(4), k_av2,  E3, 0}, {4, B(1)|B(2)|B(3), k_dxvt, E5, 0}, \
    {1, B(2)|B(3)|B(4), k_vvt,  E5, 0}, {2, B(1)|B(3)|B(4), k_avt,  E5, 0}, \
    {3, B(0)|B(2)|B(4), k_tq,   E2, R_NOCHK|R_ALT}, \
    {1, B(0)|B(2)|B(4), k_vv2,  E3, R_ABS}, {0, B(1)|B(2)|B(4), k_v0v2, E3, R_ABS} }

static const var_t lin_vars[] = {
    {"v0", "m/s",   "v0: velocity at the start", UNK},
    {"v",  "m/s",   "v: velocity at the end", UNK},
    {"a",  "m/s^2", "a: acceleration (free fall: -9.8)", UNK},
    {"t",  "s",     "t: time taken", UNK},
    {"dx", "m",     "dx: displacement (down/back = -)", UNK},
};
static const rule_t lin_rules[] = KIN_RULES("v = v0 + a t", "dx = v0 t + (1/2) a t^2",
    "v^2 = v0^2 + 2 a dx", "dx = (1/2)(v0 + v) t", "dx = v t - (1/2) a t^2");

static const var_t rot_vars[] = {
    {"w0",    "rad/s",   "w0: angular velocity at the start", UNK},
    {"w",     "rad/s",   "w: angular velocity at the end", UNK},
    {"alpha", "rad/s^2", "alpha: angular acceleration", UNK},
    {"t",     "s",       "t: time taken", UNK},
    {"dth",   "rad",     "dth: angle turned (1 rev = 6.283)", UNK},
};
static const rule_t rot_rules[] = KIN_RULES("w = w0 + alpha t", "dth = w0 t + (1/2) alpha t^2",
    "w^2 = w0^2 + 2 alpha dth", "dth = (1/2)(w0 + w) t", "dth = w t - (1/2) alpha t^2");

/* ---------- 2. projectile ------------------------------------------ */
enum { P_V0, P_TH, P_V0X, P_V0Y, P_T, P_X, P_Y, P_VX, P_VY, P_V, P_H, P_R };
FN(p_v0x,   v[P_V0] * cos(v[P_TH] * DEG))
FN(p_v0y,   v[P_V0] * sin(v[P_TH] * DEG))
FN(p_v0,    sqrt(SQ(v[P_V0X]) + SQ(v[P_V0Y])))
FN(p_th,    atan2(v[P_V0Y], v[P_V0X]) / DEG)
FN(p_v0c,   fabs(cos(v[P_TH] * DEG)) < 1e-6 ? UNK : v[P_V0X] / cos(v[P_TH] * DEG))
FN(p_v0s,   fabs(sin(v[P_TH] * DEG)) < 1e-6 ? UNK : v[P_V0Y] / sin(v[P_TH] * DEG))
FN(p_v0yt,  v[P_V0X] * tan(v[P_TH] * DEG))
FN(p_v0xt,  fabs(tan(v[P_TH] * DEG)) < 1e-6 ? UNK : v[P_V0Y] / tan(v[P_TH] * DEG))
FN(p_vx,    v[P_V0X])
FN(p_v0xv,  v[P_VX])
FN(p_x,     v[P_V0X] * v[P_T])
FN(p_tx,    v[P_X] / v[P_V0X])
FN(p_v0xx,  v[P_X] / v[P_T])
FN(p_vy,    v[P_V0Y] - GRAV * v[P_T])
FN(p_tvy,   (v[P_V0Y] - v[P_VY]) / GRAV)
FN(p_v0yv,  v[P_VY] + GRAV * v[P_T])
FN(p_y,     v[P_V0Y] * v[P_T] - 0.5 * GRAV * SQ(v[P_T]))
FN(p_v0yy,  (v[P_Y] + 0.5 * GRAV * SQ(v[P_T])) / v[P_T])
#define PDISC rootz(SQ(v[P_V0Y]) - 2 * GRAV * v[P_Y], SQ(v[P_V0Y]) + fabs(2 * GRAV * v[P_Y]))
FN(p_ty,    (v[P_V0Y] + ((flip & B(P_T)) ? -1 : 1) * PDISC) / GRAV)
FN(p_vyy,   -PDISC)
FN(p_v0yvy, rootz(SQ(v[P_VY]) + 2 * GRAV * v[P_Y], SQ(v[P_VY]) + fabs(2 * GRAV * v[P_Y])))
FN(p_v0yvt, (v[P_TH] < 0 ? -1 : 1) * rootz(SQ(v[P_VY]) + 2 * GRAV * v[P_Y], SQ(v[P_VY]) + fabs(2 * GRAV * v[P_Y])))
FN(p_v0yht, (v[P_TH] < 0 ? -1 : 1) * root(2 * GRAV * v[P_H]))
FN(p_v,     sqrt(SQ(v[P_VX]) + SQ(v[P_VY])))
FN(p_h,     SQ(v[P_V0Y]) / (2 * GRAV))
FN(p_v0yh,  root(2 * GRAV * v[P_H]))
FN(p_r,     SQ(v[P_V0]) * sin(2 * v[P_TH] * DEG) / GRAV)
FN(p_v0r,   fabs(sin(2 * v[P_TH] * DEG)) < 1e-6 ? UNK : root(v[P_R] * GRAV / sin(2 * v[P_TH] * DEG)))

static const var_t proj_vars[] = {
    {"v0",  "m/s", "v0: launch speed", UNK},
    {"th",  "deg", "th: launch angle above horizontal", UNK},
    {"v0x", "m/s", "v0x: horizontal part of v0", UNK},
    {"v0y", "m/s", "v0y: vertical part of v0 (up = +)", UNK},
    {"t",   "s",   "t: time in the air", UNK},
    {"x",   "m",   "x: horizontal distance at time t", UNK},
    {"y",   "m",   "y: height vs launch (below = -)", UNK},
    {"vx",  "m/s", "vx: horizontal velocity (constant)", UNK},
    {"vy",  "m/s", "vy: vertical velocity at time t", UNK},
    {"v",   "m/s", "v: speed at time t", UNK},
    {"H",   "m",   "H: max height above launch", UNK},
    {"R",   "m",   "R: range, ONLY if lands at launch ht", UNK},
};
static const rule_t proj_rules[] = {
    {P_V0X, B(P_V0)|B(P_TH), p_v0x, "v0x = v0 cos(th)", 0},
    {P_V0Y, B(P_V0)|B(P_TH), p_v0y, "v0y = v0 sin(th)", 0},
    {P_V0,  B(P_V0X)|B(P_V0Y), p_v0, "v0 = sqrt(v0x^2 + v0y^2)", 0},
    {P_TH,  B(P_V0X)|B(P_V0Y), p_th, "th = atan(v0y / v0x)", 0},
    {P_V0,  B(P_V0X)|B(P_TH), p_v0c, "v0x = v0 cos(th)", 0},
    {P_V0,  B(P_V0Y)|B(P_TH), p_v0s, "v0y = v0 sin(th)", 0},
    {P_V0Y, B(P_V0X)|B(P_TH), p_v0yt, "v0y = v0x tan(th)", 0},
    {P_V0X, B(P_V0Y)|B(P_TH), p_v0xt, "v0y = v0x tan(th)", 0},
    {P_VX,  B(P_V0X), p_vx, "vx = v0x (no horizontal accel)", 0},
    {P_V0X, B(P_VX), p_v0xv, "vx = v0x (no horizontal accel)", 0},
    {P_X,   B(P_V0X)|B(P_T), p_x, "x = v0x t", 0},
    {P_T,   B(P_X)|B(P_V0X), p_tx, "x = v0x t", 0},
    {P_V0X, B(P_X)|B(P_T), p_v0xx, "x = v0x t", 0},
    {P_VY,  B(P_V0Y)|B(P_T), p_vy, "vy = v0y - g t", 0},
    {P_T,   B(P_V0Y)|B(P_VY), p_tvy, "vy = v0y - g t", 0},
    {P_V0Y, B(P_VY)|B(P_T), p_v0yv, "vy = v0y - g t", 0},
    {P_Y,   B(P_V0Y)|B(P_T), p_y, "y = v0y t - (1/2) g t^2", 0},
    {P_V0Y, B(P_Y)|B(P_T), p_v0yy, "y = v0y t - (1/2) g t^2", 0},
    {P_T,   B(P_V0Y)|B(P_Y), p_ty, "y = v0y t - (1/2) g t^2", R_NOCHK|R_ALT},
    {P_VY,  B(P_V0Y)|B(P_Y), p_vyy, "vy^2 = v0y^2 - 2 g y (going down)", R_ABS},
    {P_V0Y, B(P_VY)|B(P_Y)|B(P_TH), p_v0yvt, "vy^2 = v0y^2 - 2 g y", R_ABS},
    {P_V0Y, B(P_VY)|B(P_Y), p_v0yvy, "vy^2 = v0y^2 - 2 g y", R_ABS},
    {P_V,   B(P_VX)|B(P_VY), p_v, "v = sqrt(vx^2 + vy^2)", 0},
    {P_H,   B(P_V0Y), p_h, "H = v0y^2 / (2 g)", 0},
    {P_V0Y, B(P_H)|B(P_TH), p_v0yht, "H = v0y^2 / (2 g)", R_ABS},
    {P_V0Y, B(P_H), p_v0yh, "H = v0y^2 / (2 g)", R_ABS},
    {P_R,   B(P_V0)|B(P_TH), p_r, "R = v0^2 sin(2 th) / g", 0},
    {P_V0,  B(P_R)|B(P_TH), p_v0r, "R = v0^2 sin(2 th) / g", 0},
};

/* ---------- 3. circular motion ------------------------------------- */
enum { C_R, C_V, C_T, C_W, C_RPM, C_AC, C_M, C_FC, C_MU, C_VMAX };
FN(c_v,    2 * PI * v[C_R] / v[C_T])
FN(c_t,    2 * PI * v[C_R] / v[C_V])
FN(c_rt,   v[C_V] * v[C_T] / (2 * PI))
FN(c_w,    v[C_V] / v[C_R])
FN(c_vw,   v[C_W] * v[C_R])
FN(c_rw,   v[C_V] / v[C_W])
FN(c_wt,   2 * PI / v[C_T])
FN(c_tw,   2 * PI / v[C_W])
FN(c_wrpm, v[C_RPM] * 2 * PI / 60)
FN(c_rpm,  v[C_W] * 60 / (2 * PI))
FN(c_ac,   SQ(v[C_V]) / v[C_R])
FN(c_vac,  root(v[C_AC] * v[C_R]))
FN(c_rac,  SQ(v[C_V]) / v[C_AC])
FN(c_fc,   v[C_M] * v[C_AC])
FN(c_acf,  v[C_FC] / v[C_M])
FN(c_mf,   v[C_FC] / v[C_AC])
FN(c_vmax, root(v[C_MU] * GRAV * v[C_R]))
FN(c_mu,   SQ(v[C_VMAX]) / (GRAV * v[C_R]))

static const var_t circ_vars[] = {
    {"r",    "m",     "r: radius of the circle", UNK},
    {"v",    "m/s",   "v: speed around the circle", UNK},
    {"T",    "s",     "T: time for one full circle", UNK},
    {"w",    "rad/s", "w: angular speed", UNK},
    {"rpm",  "rpm",   "rpm: revolutions per minute", UNK},
    {"ac",   "m/s^2", "ac: centripetal acceleration", UNK},
    {"m",    "kg",    "m: mass", UNK},
    {"Fc",   "N",     "Fc: net force toward the center", UNK},
    {"mu_s", "",      "mu_s: static friction (flat curve)", UNK},
    {"vmax", "m/s",   "vmax: top speed, flat curve, no skid", UNK},
};
static const rule_t circ_rules[] = {
    {C_W, B(C_RPM), c_wrpm, "w = rpm * 2 pi / 60", 0},
    {C_W, B(C_T), c_wt, "w = 2 pi / T", 0},
    {C_T, B(C_W), c_tw, "w = 2 pi / T", 0},
    {C_RPM, B(C_W), c_rpm, "w = rpm * 2 pi / 60", 0},
    {C_V, B(C_R)|B(C_T), c_v, "v = 2 pi r / T", 0},
    {C_T, B(C_R)|B(C_V), c_t, "v = 2 pi r / T", 0},
    {C_R, B(C_V)|B(C_T), c_rt, "v = 2 pi r / T", 0},
    {C_W, B(C_V)|B(C_R), c_w, "v = w r", 0},
    {C_V, B(C_W)|B(C_R), c_vw, "v = w r", 0},
    {C_R, B(C_V)|B(C_W), c_rw, "v = w r", 0},
    {C_AC, B(C_V)|B(C_R), c_ac, "ac = v^2 / r", 0},
    {C_V, B(C_AC)|B(C_R), c_vac, "ac = v^2 / r", R_ABS},
    {C_R, B(C_V)|B(C_AC), c_rac, "ac = v^2 / r", 0},
    {C_FC, B(C_M)|B(C_AC), c_fc, "Fc = m ac = m v^2 / r", 0},
    {C_AC, B(C_FC)|B(C_M), c_acf, "Fc = m ac", 0},
    {C_M, B(C_FC)|B(C_AC), c_mf, "Fc = m ac", 0},
    {C_VMAX, B(C_MU)|B(C_R), c_vmax, "vmax = sqrt(mu_s g r)", 0},
    {C_MU, B(C_VMAX)|B(C_R), c_mu, "vmax = sqrt(mu_s g r)", 0},
};

/* ---------- 4. vector components ----------------------------------- */
FN(vc_x,  v[0] * cos(v[1] * DEG))
FN(vc_y,  v[0] * sin(v[1] * DEG))
FN(vc_a,  sqrt(SQ(v[2]) + SQ(v[3])))
FN(vc_th, atan2(v[3], v[2]) / DEG)
static const var_t vec_vars[] = {
    {"A",  "",    "A: size (magnitude) of the vector", UNK},
    {"th", "deg", "th: angle from +x axis", UNK},
    {"Ax", "",    "Ax: x-component", UNK},
    {"Ay", "",    "Ay: y-component", UNK},
};
static const rule_t vec_rules[] = {
    {2, B(0)|B(1), vc_x, "Ax = A cos(th)", 0},
    {3, B(0)|B(1), vc_y, "Ay = A sin(th)", 0},
    {0, B(2)|B(3), vc_a, "A = sqrt(Ax^2 + Ay^2)", 0},
    {1, B(2)|B(3), vc_th, "th = atan(Ay / Ax)", 0},
};

/* ---------- 5. F = ma on a flat surface ---------------------------- */
enum { F_M, F_A, F_NET, F_W, F_N, F_MU, F_FK, F_APP };
FN(f_net,  v[F_M] * v[F_A])
FN(f_a,    v[F_NET] / v[F_M])
FN(f_m,    v[F_NET] / v[F_A])
FN(f_w,    v[F_M] * GRAV)
FN(f_mw,   v[F_W] / GRAV)
FN(f_n,    v[F_W])
FN(f_wn,   v[F_N])
FN(f_fk,   v[F_MU] * v[F_N])
FN(f_mu,   v[F_FK] / v[F_N])
FN(f_net2, v[F_APP] - v[F_FK])
FN(f_app,  v[F_NET] + v[F_FK])
FN(f_fk2,  v[F_APP] - v[F_NET])
static const var_t force_vars[] = {
    {"m",    "kg",    "m: mass", UNK},
    {"a",    "m/s^2", "a: acceleration", UNK},
    {"Fnet", "N",     "Fnet: net force (total)", UNK},
    {"w",    "N",     "w: weight = m g", UNK},
    {"N",    "N",     "N: normal force (flat floor)", UNK},
    {"mu_k", "",      "mu_k: kinetic friction (none: 0)", UNK},
    {"fk",   "N",     "fk: friction force", UNK},
    {"F",    "N",     "F: horizontal push or pull", UNK},
};
static const rule_t force_rules[] = {
    {F_NET, B(F_M)|B(F_A), f_net, "Fnet = m a", 0},
    {F_A, B(F_NET)|B(F_M), f_a, "Fnet = m a", 0},
    {F_M, B(F_NET)|B(F_A), f_m, "Fnet = m a", 0},
    {F_W, B(F_M), f_w, "w = m g", 0},
    {F_M, B(F_W), f_mw, "w = m g", 0},
    {F_N, B(F_W), f_n, "N = w (flat floor)", 0},
    {F_W, B(F_N), f_wn, "N = w (flat floor)", 0},
    {F_FK, B(F_MU)|B(F_N), f_fk, "fk = mu_k N", 0},
    {F_MU, B(F_FK)|B(F_N), f_mu, "fk = mu_k N", 0},
    {F_NET, B(F_APP)|B(F_FK), f_net2, "Fnet = F - fk", 0},
    {F_APP, B(F_NET)|B(F_FK), f_app, "Fnet = F - fk", 0},
    {F_FK, B(F_APP)|B(F_NET), f_fk2, "Fnet = F - fk", 0},
};

/* ---------- 6. incline --------------------------------------------- */
enum { I_M, I_TH, I_MUK, I_MUS, I_A, I_AUP, I_N, I_WP, I_FK, I_FS, I_THS };
#define ISIN sin(v[I_TH] * DEG)
#define ICOS cos(v[I_TH] * DEG)
FN(i_a,    GRAV * (ISIN - v[I_MUK] * ICOS))
FN(i_aup,  GRAV * (ISIN + v[I_MUK] * ICOS))
FN(i_muk,  (GRAV * ISIN - v[I_A]) / (GRAV * ICOS))
FN(i_muk2, (v[I_AUP] - GRAV * ISIN) / (GRAV * ICOS))
FN(i_n,    v[I_M] * GRAV * ICOS)
FN(i_wp,   v[I_M] * GRAV * ISIN)
FN(i_fk,   v[I_MUK] * v[I_N])
FN(i_fs,   v[I_MUS] * v[I_N])
FN(i_mukf, v[I_FK] / v[I_N])
FN(i_musf, v[I_FS] / v[I_N])
FN(i_mn,   v[I_N] / (GRAV * ICOS))
FN(i_mw,   v[I_WP] / (GRAV * ISIN))
FN(i_tha,  (asin(v[I_A] / GRAV / sqrt(1 + SQ(v[I_MUK]))) + atan(v[I_MUK])) / DEG)
FN(i_thu,  (asin(v[I_AUP] / GRAV / sqrt(1 + SQ(v[I_MUK]))) - atan(v[I_MUK])) / DEG)
FN(i_thn,  acosz(v[I_N] / (v[I_M] * GRAV)) / DEG)
FN(i_thw,  asinz(v[I_WP] / (v[I_M] * GRAV)) / DEG)
FN(i_ths,  atan(v[I_MUS]) / DEG)
FN(i_mus,  tan(v[I_THS] * DEG))
static const var_t incl_vars[] = {
    {"m",    "kg",    "m: mass of the block", UNK},
    {"th",   "deg",   "th: angle of the ramp", UNK},
    {"mu_k", "",      "mu_k: kinetic friction (none: 0)", UNK},
    {"mu_s", "",      "mu_s: static friction", UNK},
    {"a",    "m/s^2", "a: accel while sliding DOWN ramp", UNK},
    {"a_up", "m/s^2", "a_up: slowing rate sliding UP ramp", UNK},
    {"N",    "N",     "N: normal force = m g cos(th)", UNK},
    {"w_par","N",     "w_par: weight along ramp", UNK},
    {"fk",   "N",     "fk: kinetic friction force", UNK},
    {"fsmax","N",     "fsmax: most static friction can hold", UNK},
    {"th_s", "deg",   "th_s: angle where it starts to slip", UNK},
};
static const rule_t incl_rules[] = {
    {I_A, B(I_TH)|B(I_MUK), i_a, "a = g (sin th - mu_k cos th)", 0},
    {I_AUP, B(I_TH)|B(I_MUK), i_aup, "a_up = g (sin th + mu_k cos th)", 0},
    {I_MUK, B(I_TH)|B(I_A), i_muk, "a = g (sin th - mu_k cos th)", 0},
    {I_MUK, B(I_TH)|B(I_AUP), i_muk2, "a_up = g (sin th + mu_k cos th)", 0},
    {I_N, B(I_M)|B(I_TH), i_n, "N = m g cos(th)", 0},
    {I_WP, B(I_M)|B(I_TH), i_wp, "w_par = m g sin(th)", 0},
    {I_M, B(I_N)|B(I_TH), i_mn, "N = m g cos(th)", 0},
    {I_FK, B(I_MUK)|B(I_N), i_fk, "fk = mu_k N", 0},
    {I_FS, B(I_MUS)|B(I_N), i_fs, "fsmax = mu_s N", 0},
    {I_MUK, B(I_FK)|B(I_N), i_mukf, "fk = mu_k N", 0},
    {I_MUS, B(I_FS)|B(I_N), i_musf, "fsmax = mu_s N", 0},
    {I_THS, B(I_MUS), i_ths, "tan(th_s) = mu_s", 0},
    {I_MUS, B(I_THS), i_mus, "tan(th_s) = mu_s", 0},
    {I_M, B(I_WP)|B(I_TH), i_mw, "w_par = m g sin(th)", 0},
    {I_TH, B(I_A)|B(I_MUK), i_tha, "a = g (sin th - mu_k cos th)", 0},
    {I_TH, B(I_AUP)|B(I_MUK), i_thu, "a_up = g (sin th + mu_k cos th)", 0},
    {I_TH, B(I_N)|B(I_M), i_thn, "N = m g cos(th)", 0},
    {I_TH, B(I_WP)|B(I_M), i_thw, "w_par = m g sin(th)", 0},
};

/* ---------- 7. Atwood machine -------------------------------------- */
FN(at_a,  (v[1] - v[0]) * GRAV / (v[0] + v[1]))
FN(at_t,  2 * v[0] * v[1] * GRAV / (v[0] + v[1]))
FN(at_t1, v[0] * (GRAV + v[2]))
FN(at_t2, v[1] * (GRAV - v[2]))
FN(at_m1, v[3] / (GRAV + v[2]))
FN(at_m2, v[3] / (GRAV - v[2]))
FN(at_m2a, v[0] * (GRAV + v[2]) / (GRAV - v[2]))
FN(at_m1a, v[1] * (GRAV - v[2]) / (GRAV + v[2]))
static const var_t atw_vars[] = {
    {"m1", "kg",    "m1: lighter hanging mass", UNK},
    {"m2", "kg",    "m2: heavier hanging mass", UNK},
    {"a",  "m/s^2", "a: acceleration (m2 moves down)", UNK},
    {"T",  "N",     "T: tension in the string", UNK},
};
static const rule_t atw_rules[] = {
    {2, B(0)|B(1), at_a, "a = (m2 - m1) g / (m1 + m2)", 0},
    {3, B(0)|B(1), at_t, "T = 2 m1 m2 g / (m1 + m2)", 0},
    {3, B(0)|B(2), at_t1, "T - m1 g = m1 a", 0},
    {3, B(1)|B(2), at_t2, "m2 g - T = m2 a", 0},
    {0, B(3)|B(2), at_m1, "T - m1 g = m1 a", 0},
    {1, B(3)|B(2), at_m2, "m2 g - T = m2 a", 0},
    {1, B(0)|B(2), at_m2a, "a = (m2 - m1) g / (m1 + m2)", 0},
    {0, B(1)|B(2), at_m1a, "a = (m2 - m1) g / (m1 + m2)", 0},
};

/* ---------- 8. block on table + hanging mass ----------------------- */
/* 0=m1 (table) 1=m2 (hanging) 2=mu_k 3=a 4=T */
FN(tb_a,  (v[1] - v[2] * v[0]) * GRAV / (v[0] + v[1]))
FN(tb_t,  v[1] * (GRAV - v[3]))
FN(tb_m2, v[4] / (GRAV - v[3]))
FN(tb_t1, v[0] * (v[3] + v[2] * GRAV))
FN(tb_mu, (v[1] * GRAV - (v[0] + v[1]) * v[3]) / (v[0] * GRAV))
FN(tb_m1, v[4] / (v[3] + v[2] * GRAV))
static const var_t tab_vars[] = {
    {"m1",   "kg",    "m1: block on the table", UNK},
    {"m2",   "kg",    "m2: mass hanging off the edge", UNK},
    {"mu_k", "",      "mu_k: table friction (none: 0)", UNK},
    {"a",    "m/s^2", "a: acceleration of both", UNK},
    {"T",    "N",     "T: tension in the string", UNK},
};
static const rule_t tab_rules[] = {
    {3, B(0)|B(1)|B(2), tb_a, "a = (m2 - mu_k m1) g / (m1 + m2)", 0},
    {4, B(1)|B(3), tb_t, "m2 g - T = m2 a", 0},
    {1, B(4)|B(3), tb_m2, "m2 g - T = m2 a", 0},
    {4, B(0)|B(2)|B(3), tb_t1, "T - mu_k m1 g = m1 a", 0},
    {2, B(0)|B(1)|B(3), tb_mu, "a = (m2 - mu_k m1) g / (m1 + m2)", 0},
    {0, B(4)|B(2)|B(3), tb_m1, "T - mu_k m1 g = m1 a", 0},
};

/* ---------- 9. work and power -------------------------------------- */
enum { W_F, W_D, W_PHI, W_W, W_T, W_P, W_V };
#define WCOS cos(v[W_PHI] * DEG)
FN(w_w,   v[W_F] * v[W_D] * WCOS)
FN(w_f,   v[W_W] / (v[W_D] * WCOS))
FN(w_d,   v[W_W] / (v[W_F] * WCOS))
FN(w_p,   v[W_W] / v[W_T])
FN(w_wp,  v[W_P] * v[W_T])
FN(w_t,   v[W_W] / v[W_P])
FN(w_pv,  v[W_F] * v[W_V] * WCOS)
FN(w_vp,  v[W_P] / (v[W_F] * WCOS))
FN(w_fp,  v[W_P] / (v[W_V] * WCOS))
FN(w_phi, acosz(v[W_W] / (v[W_F] * v[W_D])) / DEG)
FN(w_phi2, acosz(v[W_P] / (v[W_F] * v[W_V])) / DEG)
static const var_t work_vars[] = {
    {"F",   "N",   "F: force doing the work", UNK},
    {"d",   "m",   "d: distance moved", UNK},
    {"phi", "deg", "phi: angle between F and motion", 0},
    {"W",   "J",   "W: work done", UNK},
    {"t",   "s",   "t: time taken", UNK},
    {"P",   "W",   "P: power (1 hp = 746 W)", UNK},
    {"v",   "m/s", "v: steady speed (for P = F v)", UNK},
};
static const rule_t work_rules[] = {
    {W_W, B(W_F)|B(W_D)|B(W_PHI), w_w, "W = F d cos(phi)", 0},
    {W_P, B(W_W)|B(W_T), w_p, "P = W / t", 0},
    {W_W, B(W_P)|B(W_T), w_wp, "P = W / t", 0},
    {W_T, B(W_W)|B(W_P), w_t, "P = W / t", 0},
    {W_P, B(W_F)|B(W_V)|B(W_PHI), w_pv, "P = F v cos(phi)", 0},
    {W_F, B(W_W)|B(W_D)|B(W_PHI), w_f, "W = F d cos(phi)", 0},
    {W_D, B(W_W)|B(W_F)|B(W_PHI), w_d, "W = F d cos(phi)", 0},
    {W_V, B(W_P)|B(W_F)|B(W_PHI), w_vp, "P = F v cos(phi)", 0},
    {W_F, B(W_P)|B(W_V)|B(W_PHI), w_fp, "P = F v cos(phi)", 0},
    {W_PHI, B(W_W)|B(W_F)|B(W_D), w_phi, "W = F d cos(phi)", 0},
    {W_PHI, B(W_P)|B(W_F)|B(W_V), w_phi2, "P = F v cos(phi)", 0},
};

/* ---------- 10. conservation of energy ----------------------------- */
enum { E_M, E_V1, E_H1, E_V2, E_H2, E_X1, E_X2, E_K, E_W };
#define E_EQ "K1 + Ug1 + Us1 + W = K2 + Ug2 + Us2"
/* extra energy per kg from spring + other work; UNK if m is needed but absent */
#define ESPR (v[E_K] * (SQ(v[E_X1]) - SQ(v[E_X2])) + 2 * v[E_W])
#define EALL (v[E_M] * (SQ(v[E_V1]) - SQ(v[E_V2])) + 2 * v[E_M] * GRAV * (v[E_H1] - v[E_H2]) + 2 * v[E_W])
#define EGH  (2 * GRAV * (fabs(v[E_H1]) + fabs(v[E_H2])))
#define ESCL (SQ(v[E_V1]) + SQ(v[E_V2]) + EGH + fabs(ESPR / v[E_M]))
#define EABS (v[E_M] * (SQ(v[E_V1]) + SQ(v[E_V2]) + EGH) + 2 * fabs(v[E_W]))
FN(e_v2,  rootz(SQ(v[E_V1]) + 2 * GRAV * (v[E_H1] - v[E_H2]) + ESPR / v[E_M], ESCL))
FN(e_v1,  rootz(SQ(v[E_V2]) - 2 * GRAV * (v[E_H1] - v[E_H2]) - ESPR / v[E_M], ESCL))
FN(e_h2,  v[E_H1] + (SQ(v[E_V1]) - SQ(v[E_V2]) + ESPR / v[E_M]) / (2 * GRAV))
FN(e_h1,  v[E_H2] - (SQ(v[E_V1]) - SQ(v[E_V2]) + ESPR / v[E_M]) / (2 * GRAV))
FN(e_v2n, ESPR != 0 ? UNK : rootz(SQ(v[E_V1]) + 2 * GRAV * (v[E_H1] - v[E_H2]), SQ(v[E_V1]) + EGH))
FN(e_v1n, ESPR != 0 ? UNK : rootz(SQ(v[E_V2]) - 2 * GRAV * (v[E_H1] - v[E_H2]), SQ(v[E_V2]) + EGH))
FN(e_h2n, ESPR != 0 ? UNK : v[E_H1] + (SQ(v[E_V1]) - SQ(v[E_V2])) / (2 * GRAV))
FN(e_h1n, ESPR != 0 ? UNK : v[E_H2] - (SQ(v[E_V1]) - SQ(v[E_V2])) / (2 * GRAV))
FN(e_w,   0.5 * v[E_M] * (SQ(v[E_V2]) - SQ(v[E_V1])) + v[E_M] * GRAV * (v[E_H2] - v[E_H1])
          + 0.5 * v[E_K] * (SQ(v[E_X2]) - SQ(v[E_X1])))
static double e_m(const double *v) {
    double num = ESPR, den = SQ(v[E_V1]) - SQ(v[E_V2]) + 2 * GRAV * (v[E_H1] - v[E_H2]);
    double nscale = fabs(v[E_K]) * (SQ(v[E_X1]) + SQ(v[E_X2])) + 2 * fabs(v[E_W]);
    if (fabs(num) <= 1e-6 * nscale || nscale == 0) return UNK;      /* mass cancels out */
    if (fabs(den) <= 1e-6 * (SQ(v[E_V1]) + SQ(v[E_V2]) + EGH)) return UNK;
    return -num / den;
}
FN(e_k,   fabs(SQ(v[E_X2]) - SQ(v[E_X1])) <= 1e-9 * (SQ(v[E_X2]) + SQ(v[E_X1])) ? UNK : EALL / (SQ(v[E_X2]) - SQ(v[E_X1])))
FN(e_x2,  rootz(SQ(v[E_X1]) + EALL / v[E_K], SQ(v[E_X1]) + EABS / fabs(v[E_K])))
FN(e_x1,  rootz(SQ(v[E_X2]) - EALL / v[E_K], SQ(v[E_X2]) + EABS / fabs(v[E_K])))
#define E_ALLBUT(x) ((uint16_t)(0x1FF & ~B(x)))
static const var_t en_vars[] = {
    {"m",  "kg",  "m: mass (often not needed)", UNK},
    {"v1", "m/s", "v1: speed at start (rest: 0)", UNK},
    {"h1", "m",   "h1: height at start", UNK},
    {"v2", "m/s", "v2: speed at end (stops: 0)", UNK},
    {"h2", "m",   "h2: height at end (ground: 0)", UNK},
    {"x1", "m",   "x1: spring squeeze at start", 0},
    {"x2", "m",   "x2: spring squeeze at end", 0},
    {"k",  "N/m", "k: spring constant", 0},
    {"W",  "J",   "W: work by friction (-) or a push", 0},
};
static const rule_t en_rules[] = {
    {E_W,  E_ALLBUT(E_W), e_w, E_EQ, 0},
    {E_V2, E_ALLBUT(E_V2), e_v2, E_EQ, R_ABS},
    {E_V1, E_ALLBUT(E_V1), e_v1, E_EQ, R_ABS},
    {E_H2, E_ALLBUT(E_H2), e_h2, E_EQ, 0},
    {E_H1, E_ALLBUT(E_H1), e_h1, E_EQ, 0},
    {E_K,  E_ALLBUT(E_K), e_k, E_EQ, 0},
    {E_X2, E_ALLBUT(E_X2), e_x2, E_EQ, R_ABS},
    {E_X1, E_ALLBUT(E_X1), e_x1, E_EQ, R_ABS},
    {E_M,  E_ALLBUT(E_M), e_m, E_EQ, 0},
    {E_V2, E_ALLBUT(E_V2) & ~B(E_M), e_v2n, "v2^2 = v1^2 + 2 g (h1 - h2)", R_ABS},
    {E_V1, E_ALLBUT(E_V1) & ~B(E_M), e_v1n, "v2^2 = v1^2 + 2 g (h1 - h2)", R_ABS},
    {E_H2, E_ALLBUT(E_H2) & ~B(E_M), e_h2n, "v2^2 = v1^2 + 2 g (h1 - h2)", 0},
    {E_H1, E_ALLBUT(E_H1) & ~B(E_M), e_h1n, "v2^2 = v1^2 + 2 g (h1 - h2)", 0},
};

/* ---------- 11-13. momentum: 0=m1 1=v1 2=m2 3=v2 4,5 = after -------- */
/* stick together: 4=vf 5=dK */
FN(s_vf,  (v[0] * v[1] + v[2] * v[3]) / (v[0] + v[2]))
FN(s_v1,  ((v[0] + v[2]) * v[4] - v[2] * v[3]) / v[0])
FN(s_v2,  ((v[0] + v[2]) * v[4] - v[0] * v[1]) / v[2])
FN(s_m1,  v[2] * (v[3] - v[4]) / (v[4] - v[1]))
FN(s_m2,  v[0] * (v[1] - v[4]) / (v[4] - v[3]))
FN(s_dk,  0.5 * (v[0] + v[2]) * SQ(v[4]) - 0.5 * v[0] * SQ(v[1]) - 0.5 * v[2] * SQ(v[3]))
static const var_t stick_vars[] = {
    {"m1", "kg",  "m1: first object", UNK},
    {"v1", "m/s", "v1: its velocity before (left = -)", UNK},
    {"m2", "kg",  "m2: second object", UNK},
    {"v2", "m/s", "v2: its velocity before (rest: 0)", UNK},
    {"vf", "m/s", "vf: shared velocity after sticking", UNK},
    {"dK", "J",   "dK: kinetic energy change (- = lost)", UNK},
};
#define S_EQ "m1 v1 + m2 v2 = (m1 + m2) vf"
static const rule_t stick_rules[] = {
    {4, B(0)|B(1)|B(2)|B(3), s_vf, S_EQ, 0},
    {1, B(0)|B(2)|B(3)|B(4), s_v1, S_EQ, 0},
    {3, B(0)|B(1)|B(2)|B(4), s_v2, S_EQ, 0},
    {0, B(1)|B(2)|B(3)|B(4), s_m1, S_EQ, 0},
    {2, B(0)|B(1)|B(3)|B(4), s_m2, S_EQ, 0},
    {5, B(0)|B(1)|B(2)|B(3)|B(4), s_dk, "dK = K after - K before", 0},
};
/* elastic: 4=v1f 5=v2f */
FN(el_v1f, ((v[0] - v[2]) * v[1] + 2 * v[2] * v[3]) / (v[0] + v[2]))
FN(el_v2f, (2 * v[0] * v[1] + (v[2] - v[0]) * v[3]) / (v[0] + v[2]))
FN(el_r2f, v[1] - v[3] + v[4])
FN(el_r1f, v[5] - v[1] + v[3])
FN(el_r1,  v[3] + v[5] - v[4])
FN(el_r2,  v[1] - v[5] + v[4])
FN(el_dp1, v[0] * (v[4] - v[1]))
FN(el_dk1, 0.5 * v[0] * (SQ(v[4]) - SQ(v[1])))
FN(g_m1,   v[2] * (v[5] - v[3]) / (v[1] - v[4]))
FN(g_m2,   v[0] * (v[1] - v[4]) / (v[5] - v[3]))
static const var_t elas_vars[] = {
    {"m1",  "kg",  "m1: first object", UNK},
    {"v1",  "m/s", "v1: its velocity before (left = -)", UNK},
    {"m2",  "kg",  "m2: second object", UNK},
    {"v2",  "m/s", "v2: its velocity before (rest: 0)", UNK},
    {"v1f", "m/s", "v1f: velocity of m1 after", UNK},
    {"v2f", "m/s", "v2f: velocity of m2 after", UNK},
    {"dp1", "kg m/s", "dp1: momentum change of m1 (m2: -dp1)", UNK},
    {"dK1", "J",   "dK1: KE change of m1 (m2 gets -dK1)", UNK},
};
#define G_EQ "m1 v1 + m2 v2 = m1 v1f + m2 v2f"
#define EL_REL "v1 - v2 = v2f - v1f (elastic)"
static const rule_t elas_rules[] = {
    {4, B(0)|B(1)|B(2)|B(3), el_v1f, "v1f = ((m1-m2) v1 + 2 m2 v2)/(m1+m2)", 0},
    {5, B(0)|B(1)|B(2)|B(3), el_v2f, "v2f = (2 m1 v1 + (m2-m1) v2)/(m1+m2)", 0},
    {5, B(1)|B(3)|B(4), el_r2f, EL_REL, 0},
    {4, B(1)|B(3)|B(5), el_r1f, EL_REL, 0},
    {1, B(3)|B(4)|B(5), el_r1, EL_REL, 0},
    {3, B(1)|B(4)|B(5), el_r2, EL_REL, 0},
    {0, B(1)|B(2)|B(3)|B(4)|B(5), g_m1, G_EQ, 0},
    {2, B(0)|B(1)|B(3)|B(4)|B(5), g_m2, G_EQ, 0},
    {6, B(0)|B(1)|B(4), el_dp1, "dp1 = m1 (v1f - v1)", 0},
    {7, B(0)|B(1)|B(4), el_dk1, "dK1 = (1/2) m1 (v1f^2 - v1^2)", 0},
};
/* any collision / recoil / explosion: one equation, six slots */
FN(g_v1,  (v[0] * v[4] + v[2] * v[5] - v[2] * v[3]) / v[0])
FN(g_v2,  (v[0] * v[4] + v[2] * v[5] - v[0] * v[1]) / v[2])
FN(g_v1f, (v[0] * v[1] + v[2] * v[3] - v[2] * v[5]) / v[0])
FN(g_v2f, (v[0] * v[1] + v[2] * v[3] - v[0] * v[4]) / v[2])
static const var_t gen_vars[] = {
    {"m1",  "kg",  "m1: first object", UNK},
    {"v1",  "m/s", "v1: before (at rest: 0, left = -)", UNK},
    {"m2",  "kg",  "m2: second object", UNK},
    {"v2",  "m/s", "v2: before (at rest: 0)", UNK},
    {"v1f", "m/s", "v1f: velocity of m1 after", UNK},
    {"v2f", "m/s", "v2f: velocity of m2 after", UNK},
    {"dK",  "J",   "dK: kinetic energy change (- = lost)", UNK},
};
FN(g_dk,  0.5 * v[0] * (SQ(v[4]) - SQ(v[1])) + 0.5 * v[2] * (SQ(v[5]) - SQ(v[3])))
static const rule_t gen_rules[] = {
    {5, B(0)|B(1)|B(2)|B(3)|B(4), g_v2f, G_EQ, 0},
    {4, B(0)|B(1)|B(2)|B(3)|B(5), g_v1f, G_EQ, 0},
    {1, B(0)|B(2)|B(3)|B(4)|B(5), g_v1, G_EQ, 0},
    {3, B(0)|B(1)|B(2)|B(4)|B(5), g_v2, G_EQ, 0},
    {0, B(1)|B(2)|B(3)|B(4)|B(5), g_m1, G_EQ, 0},
    {2, B(0)|B(1)|B(3)|B(4)|B(5), g_m2, G_EQ, 0},
    {6, B(0)|B(1)|B(2)|B(3)|B(4)|B(5), g_dk, "dK = K after - K before", 0},
};

/* ---------- 14. impulse -------------------------------------------- */
enum { J_M, J_V1, J_V2, J_J, J_F, J_T, J_HD, J_HU };
FN(j_j,   v[J_M] * (v[J_V2] - v[J_V1]))
FN(j_v2,  v[J_V1] + v[J_J] / v[J_M])
FN(j_v1,  v[J_V2] - v[J_J] / v[J_M])
FN(j_m,   v[J_J] / (v[J_V2] - v[J_V1]))
FN(j_jf,  v[J_F] * v[J_T])
FN(j_f,   v[J_J] / v[J_T])
FN(j_t,   v[J_J] / v[J_F])
FN(j_vd,  -root(2 * GRAV * v[J_HD]))
FN(j_vu,  root(2 * GRAV * v[J_HU]))
FN(j_hd,  v[J_V1] < 0 ? SQ(v[J_V1]) / (2 * GRAV) : UNK)
FN(j_hu,  v[J_V2] > 0 ? SQ(v[J_V2]) / (2 * GRAV) : UNK)
static const var_t imp_vars[] = {
    {"m",  "kg",   "m: mass", UNK},
    {"v1", "m/s",  "v1: velocity before (left = -)", UNK},
    {"v2", "m/s",  "v2: velocity after (bounce back = -)", UNK},
    {"J",  "kg m/s", "J: impulse = change in momentum", UNK},
    {"F",  "N",    "F: average force", UNK},
    {"t",  "s",    "t: contact time (ms: use E-3)", UNK},
    {"h_drop", "m", "h_drop: fell from this height", UNK},
    {"h_up", "m",  "h_up: bounces back up this high", UNK},
};
static const rule_t imp_rules[] = {
    {J_V1, B(J_HD), j_vd, "v1 = -sqrt(2 g h_drop) (downward)", 0},
    {J_V2, B(J_HU), j_vu, "v2 = +sqrt(2 g h_up) (upward)", 0},
    {J_J, B(J_M)|B(J_V1)|B(J_V2), j_j, "J = m (v2 - v1)", 0},
    {J_J, B(J_F)|B(J_T), j_jf, "J = F t", 0},
    {J_F, B(J_J)|B(J_T), j_f, "J = F t", 0},
    {J_T, B(J_J)|B(J_F), j_t, "J = F t", 0},
    {J_V2, B(J_M)|B(J_V1)|B(J_J), j_v2, "J = m (v2 - v1)", 0},
    {J_V1, B(J_M)|B(J_V2)|B(J_J), j_v1, "J = m (v2 - v1)", 0},
    {J_M, B(J_V1)|B(J_V2)|B(J_J), j_m, "J = m (v2 - v1)", 0},
    {J_HD, B(J_V1), j_hd, "v1 = -sqrt(2 g h_drop) (downward)", 0},
    {J_HU, B(J_V2), j_hu, "v2 = +sqrt(2 g h_up) (upward)", 0},
};

/* ---------- momentum and kinetic energy of one object -------------- */
FN(mk_p,   v[0] * v[1])
FN(mk_v,   v[2] / v[0])
FN(mk_m,   v[2] / v[1])
FN(mk_k,   0.5 * v[0] * SQ(v[1]))
FN(mk_kp,  SQ(v[2]) / (2 * v[0]))
FN(mk_vk,  root(2 * v[3] / v[0]))
FN(mk_mk,  2 * v[3] / SQ(v[1]))
FN(mk_pk,  root(2 * v[0] * v[3]))
FN(mk_mp,  SQ(v[2]) / (2 * v[3]))
static const var_t mk_vars[] = {
    {"m", "kg",     "m: mass (weighs 625 N: 625/9.8)", UNK},
    {"v", "m/s",    "v: velocity", UNK},
    {"p", "kg m/s", "p: momentum", UNK},
    {"K", "J",      "K: kinetic energy", UNK},
};
static const rule_t mk_rules[] = {
    {2, B(0)|B(1), mk_p, "p = m v", 0},
    {1, B(2)|B(0), mk_v, "p = m v", 0},
    {0, B(2)|B(1), mk_m, "p = m v", 0},
    {3, B(0)|B(1), mk_k, "K = (1/2) m v^2", 0},
    {3, B(2)|B(0), mk_kp, "K = p^2 / (2 m)", 0},
    {1, B(3)|B(0), mk_vk, "K = (1/2) m v^2", R_ABS},
    {0, B(3)|B(1), mk_mk, "K = (1/2) m v^2", 0},
    {2, B(0)|B(3), mk_pk, "K = p^2 / (2 m)", R_ABS},
    {0, B(2)|B(3), mk_mp, "K = p^2 / (2 m)", 0},
};

/* ---------- impulse in two dimensions ------------------------------ */
enum { D_M, D_V1, D_A1, D_V2, D_A2, D_T, D_JX, D_JY, D_J, D_FX, D_FY, D_F };
FN(d_jx,  v[D_M] * (v[D_V2] * cos(v[D_A2] * DEG) - v[D_V1] * cos(v[D_A1] * DEG)))
FN(d_jy,  v[D_M] * (v[D_V2] * sin(v[D_A2] * DEG) - v[D_V1] * sin(v[D_A1] * DEG)))
FN(d_j,   sqrt(SQ(v[D_JX]) + SQ(v[D_JY])))
FN(d_fx,  v[D_JX] / v[D_T])
FN(d_fy,  v[D_JY] / v[D_T])
FN(d_f,   v[D_J] / v[D_T])
#define D2X (v[D_JX] / v[D_M] + v[D_V1] * cos(v[D_A1] * DEG))
#define D2Y (v[D_JY] / v[D_M] + v[D_V1] * sin(v[D_A1] * DEG))
#define D1X (v[D_V2] * cos(v[D_A2] * DEG) - v[D_JX] / v[D_M])
#define D1Y (v[D_V2] * sin(v[D_A2] * DEG) - v[D_JY] / v[D_M])
#define DDX (v[D_V2] * cos(v[D_A2] * DEG) - v[D_V1] * cos(v[D_A1] * DEG))
#define DDY (v[D_V2] * sin(v[D_A2] * DEG) - v[D_V1] * sin(v[D_A1] * DEG))
FN(d_jxf, v[D_FX] * v[D_T])
FN(d_jyf, v[D_FY] * v[D_T])
FN(d_jf,  v[D_F] * v[D_T])
FN(d_t,   v[D_J] / v[D_F])
FN(d_v2,  sqrt(SQ(D2X) + SQ(D2Y)))
FN(d_a2,  atan2(D2Y, D2X) / DEG)
FN(d_v1,  sqrt(SQ(D1X) + SQ(D1Y)))
FN(d_a1,  atan2(D1Y, D1X) / DEG)
FN(d_m,   v[D_J] / sqrt(SQ(DDX) + SQ(DDY)))
static const var_t imp2_vars[] = {
    {"m",  "kg",  "m: mass", UNK},
    {"v1", "m/s", "v1: speed before", UNK},
    {"a1", "deg", "a1: direction before (right = 0)", 0},
    {"v2", "m/s", "v2: speed after", UNK},
    {"a2", "deg", "a2: direction after (left = 180)", UNK},
    {"t",  "s",   "t: contact time (ms: use E-3)", UNK},
    {"Jx", "kg m/s", "Jx: horizontal impulse", UNK},
    {"Jy", "kg m/s", "Jy: vertical impulse", UNK},
    {"J",  "kg m/s", "J: size of the impulse", UNK},
    {"Fx", "N",   "Fx: horizontal average force", UNK},
    {"Fy", "N",   "Fy: vertical average force", UNK},
    {"F",  "N",   "F: size of the average force", UNK},
};
#define D_BOTH (B(D_M)|B(D_V1)|B(D_A1)|B(D_V2)|B(D_A2))
static const rule_t imp2_rules[] = {
    {D_JX, D_BOTH, d_jx, "Jx = m (v2x - v1x)", 0},
    {D_JY, D_BOTH, d_jy, "Jy = m (v2y - v1y)", 0},
    {D_J,  B(D_JX)|B(D_JY), d_j, "J = sqrt(Jx^2 + Jy^2)", 0},
    {D_FX, B(D_JX)|B(D_T), d_fx, "Fx = Jx / t", 0},
    {D_FY, B(D_JY)|B(D_T), d_fy, "Fy = Jy / t", 0},
    {D_F,  B(D_J)|B(D_T), d_f, "F = J / t", 0},
    {D_JX, B(D_FX)|B(D_T), d_jxf, "Jx = Fx t", 0},
    {D_JY, B(D_FY)|B(D_T), d_jyf, "Jy = Fy t", 0},
    {D_J,  B(D_F)|B(D_T), d_jf, "J = F t", 0},
    {D_T,  B(D_J)|B(D_F), d_t, "J = F t", 0},
    {D_V2, B(D_M)|B(D_V1)|B(D_A1)|B(D_JX)|B(D_JY), d_v2, "v2 = v1 + J / m (x and y parts)", 0},
    {D_A2, B(D_M)|B(D_V1)|B(D_A1)|B(D_JX)|B(D_JY), d_a2, "tan(a2) = v2y / v2x", R_NOCHK},
    {D_V1, B(D_M)|B(D_V2)|B(D_A2)|B(D_JX)|B(D_JY), d_v1, "v1 = v2 - J / m (x and y parts)", 0},
    {D_A1, B(D_M)|B(D_V2)|B(D_A2)|B(D_JX)|B(D_JY), d_a1, "tan(a1) = v1y / v1x", R_NOCHK},
    {D_M,  B(D_J)|B(D_V1)|B(D_A1)|B(D_V2)|B(D_A2), d_m, "J = m |v2 - v1|", 0},
};

/* ---------- collision in two dimensions ---------------------------- */
enum { T_M1, T_V1, T_A1, T_M2, T_V2, T_A2, T_V1F, T_B1, T_V2F, T_B2, T_DK };
#define TPX (v[T_M1] * v[T_V1] * cos(v[T_A1] * DEG) + v[T_M2] * v[T_V2] * cos(v[T_A2] * DEG))
#define TPY (v[T_M1] * v[T_V1] * sin(v[T_A1] * DEG) + v[T_M2] * v[T_V2] * sin(v[T_A2] * DEG))
#define T2X (TPX - v[T_M1] * v[T_V1F] * cos(v[T_B1] * DEG))
#define T2Y (TPY - v[T_M1] * v[T_V1F] * sin(v[T_B1] * DEG))
#define T1X (TPX - v[T_M2] * v[T_V2F] * cos(v[T_B2] * DEG))
#define T1Y (TPY - v[T_M2] * v[T_V2F] * sin(v[T_B2] * DEG))
#define TSIN sin((v[T_B2] - v[T_B1]) * DEG)
FN(t_v2f,  sqrt(SQ(T2X) + SQ(T2Y)) / v[T_M2])
FN(t_b2,   atan2(T2Y, T2X) / DEG)
FN(t_v1f,  sqrt(SQ(T1X) + SQ(T1Y)) / v[T_M1])
FN(t_b1,   atan2(T1Y, T1X) / DEG)
FN(t_v1fa, fabs(TSIN) < 1e-6 ? UNK : (TPX * sin(v[T_B2] * DEG) - TPY * cos(v[T_B2] * DEG)) / (v[T_M1] * TSIN))
FN(t_v2fa, fabs(TSIN) < 1e-6 ? UNK : (TPY * cos(v[T_B1] * DEG) - TPX * sin(v[T_B1] * DEG)) / (v[T_M2] * TSIN))
FN(t_dk,   0.5 * v[T_M1] * (SQ(v[T_V1F]) - SQ(v[T_V1])) + 0.5 * v[T_M2] * (SQ(v[T_V2F]) - SQ(v[T_V2])))
/* change in velocity of each object, x and y */
#define TD1X (v[T_V1] * cos(v[T_A1] * DEG) - v[T_V1F] * cos(v[T_B1] * DEG))
#define TD1Y (v[T_V1] * sin(v[T_A1] * DEG) - v[T_V1F] * sin(v[T_B1] * DEG))
#define TD2X (v[T_V2F] * cos(v[T_B2] * DEG) - v[T_V2] * cos(v[T_A2] * DEG))
#define TD2Y (v[T_V2F] * sin(v[T_B2] * DEG) - v[T_V2] * sin(v[T_A2] * DEG))
FN(t_m2,   fabs(TD2X) > fabs(TD2Y) ? v[T_M1] * TD1X / TD2X : v[T_M1] * TD1Y / TD2Y)
FN(t_m1,   fabs(TD1X) > fabs(TD1Y) ? v[T_M2] * TD2X / TD1X : v[T_M2] * TD2Y / TD1Y)
static double t_v1(const double *v) {
    double c = cos(v[T_A1] * DEG), s = sin(v[T_A1] * DEG);
    double X = v[T_M1] * v[T_V1F] * cos(v[T_B1] * DEG) + v[T_M2] * TD2X;
    double Y = v[T_M1] * v[T_V1F] * sin(v[T_B1] * DEG) + v[T_M2] * TD2Y;
    return fabs(c) > fabs(s) ? X / (v[T_M1] * c) : Y / (v[T_M1] * s);
}
static double t_v2(const double *v) {
    double c = cos(v[T_A2] * DEG), s = sin(v[T_A2] * DEG);
    double X = v[T_M2] * v[T_V2F] * cos(v[T_B2] * DEG) - v[T_M1] * TD1X;
    double Y = v[T_M2] * v[T_V2F] * sin(v[T_B2] * DEG) - v[T_M1] * TD1Y;
    return fabs(c) > fabs(s) ? X / (v[T_M2] * c) : Y / (v[T_M2] * s);
}
static const var_t col2_vars[] = {
    {"m1",  "kg",  "m1: first object (equal masses: 1)", UNK},
    {"v1",  "m/s", "v1: its speed before", UNK},
    {"a1",  "deg", "a1: its direction before (blank = 0)", 0},
    {"m2",  "kg",  "m2: second object", UNK},
    {"v2",  "m/s", "v2: its speed before (blank = rest)", 0},
    {"a2",  "deg", "a2: its direction before", 0},
    {"v1f", "m/s", "v1f: speed of m1 after", UNK},
    {"b1",  "deg", "b1: direction of m1 after", UNK},
    {"v2f", "m/s", "v2f: speed of m2 after", UNK},
    {"b2",  "deg", "b2: direction of m2 after (below: -)", UNK},
    {"dK",  "J",   "dK: kinetic energy change (- = lost)", UNK},
};
#define T_BEF (B(T_M1)|B(T_V1)|B(T_A1)|B(T_M2)|B(T_V2)|B(T_A2))
#define T_EQ "px and py are each conserved"
static const rule_t col2_rules[] = {
    {T_V2F, T_BEF|B(T_V1F)|B(T_B1), t_v2f, T_EQ, 0},
    {T_B2,  T_BEF|B(T_V1F)|B(T_B1), t_b2, "tan(b2) = p2y / p2x", R_NOCHK},
    {T_V1F, T_BEF|B(T_V2F)|B(T_B2), t_v1f, T_EQ, 0},
    {T_B1,  T_BEF|B(T_V2F)|B(T_B2), t_b1, "tan(b1) = p1y / p1x", R_NOCHK},
    {T_V1F, T_BEF|B(T_B1)|B(T_B2), t_v1fa, T_EQ " (2 eqns)", 0},
    {T_V2F, T_BEF|B(T_B1)|B(T_B2), t_v2fa, T_EQ " (2 eqns)", 0},
    {T_DK,  T_BEF|B(T_V1F)|B(T_V2F), t_dk, "dK = K after - K before", 0},
#define T_AFT (B(T_V1F)|B(T_B1)|B(T_V2F)|B(T_B2))
    {T_V1,  (T_BEF & ~B(T_V1))|T_AFT, t_v1, T_EQ, 0},
    {T_V2,  (T_BEF & ~B(T_V2))|T_AFT, t_v2, T_EQ, 0},
    {T_M2,  (T_BEF & ~B(T_M2))|T_AFT, t_m2, T_EQ, 0},
    {T_M1,  (T_BEF & ~B(T_M1))|T_AFT, t_m1, T_EQ, 0},
};

/* ---------- stick together in two dimensions ----------------------- */
enum { X_M1, X_V1, X_A1, X_M2, X_V2, X_A2, X_VF, X_AF, X_DK };
#define XPX (v[X_M1] * v[X_V1] * cos(v[X_A1] * DEG) + v[X_M2] * v[X_V2] * cos(v[X_A2] * DEG))
#define XPY (v[X_M1] * v[X_V1] * sin(v[X_A1] * DEG) + v[X_M2] * v[X_V2] * sin(v[X_A2] * DEG))
#define XS1 sin((v[X_AF] - v[X_A1]) * DEG)
#define XS2 sin((v[X_AF] - v[X_A2]) * DEG)
FN(x_vf,  sqrt(SQ(XPX) + SQ(XPY)) / (v[X_M1] + v[X_M2]))
FN(x_af,  atan2(XPY, XPX) / DEG)
FN(x_v1,  fabs(XS1) < 2e-3 ? UNK : -v[X_M2] * v[X_V2] * XS2 / (v[X_M1] * XS1))
FN(x_v2,  fabs(XS2) < 2e-3 ? UNK : -v[X_M1] * v[X_V1] * XS1 / (v[X_M2] * XS2))
#define XFX ((v[X_M1] + v[X_M2]) * v[X_VF] * cos(v[X_AF] * DEG))
#define XFY ((v[X_M1] + v[X_M2]) * v[X_VF] * sin(v[X_AF] * DEG))
#define X2X (XFX - v[X_M1] * v[X_V1] * cos(v[X_A1] * DEG))
#define X2Y (XFY - v[X_M1] * v[X_V1] * sin(v[X_A1] * DEG))
#define X1X (XFX - v[X_M2] * v[X_V2] * cos(v[X_A2] * DEG))
#define X1Y (XFY - v[X_M2] * v[X_V2] * sin(v[X_A2] * DEG))
FN(x_v2b, sqrt(SQ(X2X) + SQ(X2Y)) / v[X_M2])
FN(x_a2b, atan2(X2Y, X2X) / DEG)
FN(x_v1b, sqrt(SQ(X1X) + SQ(X1Y)) / v[X_M1])
FN(x_a1b, atan2(X1Y, X1X) / DEG)
FN(x_m1,  -v[X_M2] * v[X_V2] * XS2 / (v[X_V1] * XS1))
FN(x_m2,  -v[X_M1] * v[X_V1] * XS1 / (v[X_V2] * XS2))
FN(x_dk,  0.5 * (v[X_M1] + v[X_M2]) * SQ(v[X_VF]) - 0.5 * v[X_M1] * SQ(v[X_V1]) - 0.5 * v[X_M2] * SQ(v[X_V2]))
static const var_t stk2_vars[] = {
    {"m1", "kg",  "m1: first object", UNK},
    {"v1", "m/s", "v1: its speed before", UNK},
    {"a1", "deg", "a1: its direction (east = 0)", 0},
    {"m2", "kg",  "m2: second object", UNK},
    {"v2", "m/s", "v2: its speed before", UNK},
    {"a2", "deg", "a2: its direction (north = 90)", UNK},
    {"vf", "m/s", "vf: speed together after", UNK},
    {"af", "deg", "af: direction together after", UNK},
    {"dK", "J",   "dK: kinetic energy change (- = lost)", UNK},
};
#define X_BEF (B(X_M1)|B(X_V1)|B(X_A1)|B(X_M2)|B(X_V2)|B(X_A2))
static const rule_t stk2_rules[] = {
    {X_VF, X_BEF, x_vf, "(m1+m2) vf = sqrt(px^2 + py^2)", 0},
    {X_AF, X_BEF, x_af, "tan(af) = py / px", R_NOCHK},
    {X_DK, X_BEF|B(X_VF), x_dk, "dK = K after - K before", 0},
    {X_M1, (X_BEF & ~B(X_M1))|B(X_AF), x_m1, "no momentum sideways to af", 0},
    {X_M2, (X_BEF & ~B(X_M2))|B(X_AF), x_m2, "no momentum sideways to af", 0},
#define X_AFT (B(X_VF)|B(X_AF)|B(X_M1)|B(X_M2))
    {X_V2, X_AFT|B(X_V1)|B(X_A1), x_v2b, "p2 = p after - p1 (x and y parts)", 0},
    {X_A2, X_AFT|B(X_V1)|B(X_A1), x_a2b, "tan(a2) = p2y / p2x", R_NOCHK},
    {X_V1, X_AFT|B(X_V2)|B(X_A2), x_v1b, "p1 = p after - p2 (x and y parts)", 0},
    {X_A1, X_AFT|B(X_V2)|B(X_A2), x_a1b, "tan(a1) = p1y / p1x", R_NOCHK},
    {X_V1, (X_BEF & ~B(X_V1))|B(X_AF), x_v1, "no momentum sideways to af", 0},
    {X_V2, (X_BEF & ~B(X_V2))|B(X_AF), x_v2, "no momentum sideways to af", 0},
};

/* ---------- spring pushing two blocks apart ------------------------ */
/* 0=m1 1=m2 2=k 3=x 4=U 5=v1 6=v2 (speeds, opposite directions) */
FN(sp_u,   0.5 * v[2] * SQ(v[3]))
FN(sp_k,   2 * v[4] / SQ(v[3]))
FN(sp_x,   root(2 * v[4] / v[2]))
FN(sp_v1,  root(2 * v[4] * v[1] / (v[0] * (v[0] + v[1]))))
FN(sp_v2,  root(2 * v[4] * v[0] / (v[1] * (v[0] + v[1]))))
FN(sp_v1m, v[1] * v[6] / v[0])
FN(sp_v2m, v[0] * v[5] / v[1])
FN(sp_uk,  0.5 * v[0] * SQ(v[5]) + 0.5 * v[1] * SQ(v[6]))
FN(sp_m1,  v[1] * v[6] / v[5])
FN(sp_m2,  v[0] * v[5] / v[6])
static const var_t spr2_vars[] = {
    {"m1", "kg",  "m1: first block", UNK},
    {"m2", "kg",  "m2: second block", UNK},
    {"k",  "N/m", "k: spring (1 N/cm = 100 N/m)", UNK},
    {"x",  "m",   "x: how far it was squeezed", UNK},
    {"U",  "J",   "U: energy stored in the spring", UNK},
    {"v1", "m/s", "v1: speed of m1 after release", UNK},
    {"v2", "m/s", "v2: speed of m2 (other way)", UNK},
};
static const rule_t spr2_rules[] = {
    {4, B(2)|B(3), sp_u, "U = (1/2) k x^2", 0},
    {5, B(0)|B(1)|B(6), sp_v1m, "m1 v1 = m2 v2", 0},
    {6, B(0)|B(1)|B(5), sp_v2m, "m1 v1 = m2 v2", 0},
    {4, B(0)|B(1)|B(5)|B(6), sp_uk, "U = (1/2) m1 v1^2 + (1/2) m2 v2^2", 0},
    {5, B(0)|B(1)|B(4), sp_v1, "U = K1 + K2 with m1 v1 = m2 v2", R_ABS},
    {6, B(0)|B(1)|B(4), sp_v2, "U = K1 + K2 with m1 v1 = m2 v2", R_ABS},
    {2, B(4)|B(3), sp_k, "U = (1/2) k x^2", 0},
    {3, B(4)|B(2), sp_x, "U = (1/2) k x^2", R_ABS},
    {0, B(1)|B(5)|B(6), sp_m1, "m1 v1 = m2 v2", 0},
    {1, B(0)|B(5)|B(6), sp_m2, "m1 v1 = m2 v2", 0},
};

/* ---------- stick together, then swing / slide / spring / fall ------ */
enum { Z_M1, Z_V1, Z_M2, Z_VF, Z_H, Z_MU, Z_D, Z_K, Z_X, Z_HT, Z_XR, Z_L, Z_T, Z_K1, Z_KF };
#define ZM (v[Z_M1] + v[Z_M2])
FN(z_vf,   v[Z_M1] * v[Z_V1] / ZM)
FN(z_v1,   ZM * v[Z_VF] / v[Z_M1])
FN(z_m2,   v[Z_M1] * v[Z_V1] / v[Z_VF] - v[Z_M1])
FN(z_h,    SQ(v[Z_VF]) / (2 * GRAV))
FN(z_vfh,  root(2 * GRAV * v[Z_H]))
FN(z_d,    SQ(v[Z_VF]) / (2 * v[Z_MU] * GRAV))
FN(z_vfd,  root(2 * v[Z_MU] * GRAV * v[Z_D]))
FN(z_mu,   SQ(v[Z_VF]) / (2 * GRAV * v[Z_D]))
FN(z_x,    v[Z_VF] * root(ZM / v[Z_K]))
FN(z_vfx,  v[Z_X] * root(v[Z_K] / ZM))
FN(z_k,    ZM * SQ(v[Z_VF]) / SQ(v[Z_X]))
FN(z_xr,   v[Z_VF] * root(2 * v[Z_HT] / GRAV))
FN(z_vfr,  v[Z_XR] / root(2 * v[Z_HT] / GRAV))
FN(z_t,    ZM * (GRAV + SQ(v[Z_VF]) / v[Z_L]))
FN(z_ht,   GRAV * SQ(v[Z_XR]) / (2 * SQ(v[Z_VF])))
FN(z_l,    ZM * SQ(v[Z_VF]) / (v[Z_T] - ZM * GRAV))
FN(z_m1,   v[Z_M2] * v[Z_VF] / (v[Z_V1] - v[Z_VF]))
FN(z_k1,   0.5 * v[Z_M1] * SQ(v[Z_V1]))
FN(z_kf,   0.5 * ZM * SQ(v[Z_VF]))
static const var_t then_vars[] = {
    {"m1",  "kg",  "m1: bullet / clay / moving object", UNK},
    {"v1",  "m/s", "v1: its speed before it hits", UNK},
    {"m2",  "kg",  "m2: block it sticks in (at rest)", UNK},
    {"vf",  "m/s", "vf: speed together right after", UNK},
    {"h",   "m",   "h: THEN swings/climbs up this high", UNK},
    {"mu_k","",    "mu_k: THEN slides with this friction", UNK},
    {"d",   "m",   "d: ...and stops after this distance", UNK},
    {"k",   "N/m", "k: THEN hits a spring this stiff", UNK},
    {"x",   "m",   "x: ...and squeezes it this far", UNK},
    {"Ht",  "m",   "Ht: THEN flies off a table this high", UNK},
    {"xr",  "m",   "xr: ...and lands this far out", UNK},
    {"L",   "m",   "L: hangs from a cord this long", UNK},
    {"T",   "N",   "T: cord tension right after impact", UNK},
    {"K1",  "J",   "K1: kinetic energy of m1 before", UNK},
    {"Kf",  "J",   "Kf: kinetic energy together after", UNK},
};
static const rule_t then_rules[] = {
    {Z_VF, B(Z_M1)|B(Z_V1)|B(Z_M2), z_vf, "m1 v1 = (m1 + m2) vf", 0},
    {Z_VF, B(Z_H), z_vfh, "(1/2) vf^2 = g h", R_ABS},
    {Z_VF, B(Z_MU)|B(Z_D), z_vfd, "(1/2) vf^2 = mu_k g d", R_ABS},
    {Z_VF, B(Z_X)|B(Z_K)|B(Z_M1)|B(Z_M2), z_vfx, "(1/2)(m1+m2) vf^2 = (1/2) k x^2", R_ABS},
    {Z_VF, B(Z_XR)|B(Z_HT), z_vfr, "xr = vf sqrt(2 Ht / g)", 0},
    {Z_V1, B(Z_M1)|B(Z_M2)|B(Z_VF), z_v1, "m1 v1 = (m1 + m2) vf", 0},
    {Z_M2, B(Z_M1)|B(Z_V1)|B(Z_VF), z_m2, "m1 v1 = (m1 + m2) vf", 0},
    {Z_H,  B(Z_VF), z_h, "(1/2) vf^2 = g h", 0},
    {Z_D,  B(Z_VF)|B(Z_MU), z_d, "(1/2) vf^2 = mu_k g d", 0},
    {Z_MU, B(Z_VF)|B(Z_D), z_mu, "(1/2) vf^2 = mu_k g d", 0},
    {Z_X,  B(Z_VF)|B(Z_K)|B(Z_M1)|B(Z_M2), z_x, "(1/2)(m1+m2) vf^2 = (1/2) k x^2", 0},
    {Z_K,  B(Z_VF)|B(Z_X)|B(Z_M1)|B(Z_M2), z_k, "(1/2)(m1+m2) vf^2 = (1/2) k x^2", 0},
    {Z_XR, B(Z_VF)|B(Z_HT), z_xr, "xr = vf sqrt(2 Ht / g)", 0},
    {Z_T,  B(Z_VF)|B(Z_L)|B(Z_M1)|B(Z_M2), z_t, "T - M g = M vf^2 / L", 0},
    {Z_K1, B(Z_M1)|B(Z_V1), z_k1, "K1 = (1/2) m1 v1^2", 0},
    {Z_KF, B(Z_M1)|B(Z_M2)|B(Z_VF), z_kf, "Kf = (1/2)(m1 + m2) vf^2", 0},
    {Z_HT, B(Z_XR)|B(Z_VF), z_ht, "xr = vf sqrt(2 Ht / g)", 0},
    {Z_L,  B(Z_T)|B(Z_VF)|B(Z_M1)|B(Z_M2), z_l, "T - M g = M vf^2 / L", 0},
    {Z_M1, B(Z_V1)|B(Z_M2)|B(Z_VF), z_m1, "m1 v1 = (m1 + m2) vf", 0},
};

/* ---------- center of mass of two objects -------------------------- */
enum { O_M1, O_X1, O_M2, O_X2, O_XCM, O_V1, O_V2, O_VCM, O_P };
#define OM (v[O_M1] + v[O_M2])
FN(o_xcm,  v[O_X1] + (v[O_X2] - v[O_X1]) * (v[O_M2] / OM))
FN(o_x2,   v[O_X1] + (v[O_XCM] - v[O_X1]) * (OM / v[O_M2]))
FN(o_x1,   v[O_X2] - (v[O_X2] - v[O_XCM]) * (OM / v[O_M1]))
FN(o_m1,   v[O_M2] * (v[O_X2] - v[O_XCM]) / (v[O_XCM] - v[O_X1]))
FN(o_m2,   v[O_M1] * (v[O_XCM] - v[O_X1]) / (v[O_X2] - v[O_XCM]))
FN(o_vcm,  v[O_V1] + (v[O_V2] - v[O_V1]) * (v[O_M2] / OM))
FN(o_v1,   v[O_V2] - (v[O_V2] - v[O_VCM]) * (OM / v[O_M1]))
FN(o_v2,   v[O_V1] + (v[O_VCM] - v[O_V1]) * (OM / v[O_M2]))
FN(o_p,    OM * v[O_VCM])
FN(o_vp,   v[O_P] / OM)
static const var_t com_vars[] = {
    {"m1",  "kg",  "m1: first object", UNK},
    {"x1",  "m",   "x1: where it is (blank = 0)", 0},
    {"m2",  "kg",  "m2: second object", UNK},
    {"x2",  "m",   "x2: where it is (same axis)", UNK},
    {"xcm", "m",   "xcm: center of mass position", UNK},
    {"v1",  "m/s", "v1: velocity of m1", UNK},
    {"v2",  "m/s", "v2: velocity of m2", UNK},
    {"vcm", "m/s", "vcm: velocity of the center of mass", UNK},
    {"P",   "kg m/s", "P: total momentum = M vcm", UNK},
};
static const rule_t com_rules[] = {
    {O_XCM, B(O_M1)|B(O_X1)|B(O_M2)|B(O_X2), o_xcm, "xcm = (m1 x1 + m2 x2)/(m1 + m2)", 0},
    {O_VCM, B(O_M1)|B(O_V1)|B(O_M2)|B(O_V2), o_vcm, "vcm = (m1 v1 + m2 v2)/(m1 + m2)", 0},
    {O_P,   B(O_M1)|B(O_M2)|B(O_VCM), o_p, "P = (m1 + m2) vcm", 0},
    {O_VCM, B(O_M1)|B(O_M2)|B(O_P), o_vp, "P = (m1 + m2) vcm", 0},
    {O_M1,  B(O_X1)|B(O_M2)|B(O_X2)|B(O_XCM), o_m1, "xcm = (m1 x1 + m2 x2)/(m1 + m2)", 0},
    {O_M2,  B(O_M1)|B(O_X1)|B(O_X2)|B(O_XCM), o_m2, "xcm = (m1 x1 + m2 x2)/(m1 + m2)", 0},
    {O_X2,  B(O_M1)|B(O_X1)|B(O_M2)|B(O_XCM), o_x2, "xcm = (m1 x1 + m2 x2)/(m1 + m2)", 0},
    {O_X1,  B(O_M1)|B(O_M2)|B(O_X2)|B(O_XCM), o_x1, "xcm = (m1 x1 + m2 x2)/(m1 + m2)", 0},
    {O_V1,  B(O_M1)|B(O_M2)|B(O_V2)|B(O_VCM), o_v1, "vcm = (m1 v1 + m2 v2)/(m1 + m2)", 0},
    {O_V2,  B(O_M1)|B(O_M2)|B(O_V1)|B(O_VCM), o_v2, "vcm = (m1 v1 + m2 v2)/(m1 + m2)", 0},
};

/* ---------- rocket thrust ------------------------------------------ */
enum { RK_VEX, RK_R, RK_F, RK_M, RK_A, RK_T, RK_DM };
FN(rk_f,   v[RK_VEX] * v[RK_R])
FN(rk_r,   v[RK_F] / v[RK_VEX])
FN(rk_vex, v[RK_F] / v[RK_R])
FN(rk_a,   v[RK_F] / v[RK_M])
FN(rk_fa,  v[RK_M] * v[RK_A])
FN(rk_m,   v[RK_F] / v[RK_A])
FN(rk_dm,  v[RK_R] * v[RK_T])
FN(rk_rt,  v[RK_DM] / v[RK_T])
FN(rk_t,   v[RK_DM] / v[RK_R])
static const var_t rock_vars[] = {
    {"vex", "m/s",  "vex: exhaust speed relative to rocket", UNK},
    {"R",   "kg/s", "R: fuel burned per second", UNK},
    {"F",   "N",    "F: thrust", UNK},
    {"m",   "kg",   "m: total mass right now", UNK},
    {"a",   "m/s^2","a: acceleration (no gravity)", UNK},
    {"t",   "s",    "t: how long it fires", UNK},
    {"dm",  "kg",   "dm: fuel or gas used in time t", UNK},
};
static const rule_t rock_rules[] = {
    {RK_F, B(RK_VEX)|B(RK_R), rk_f, "F = vex R", 0},
    {RK_F, B(RK_M)|B(RK_A), rk_fa, "F = m a", 0},
    {RK_R, B(RK_DM)|B(RK_T), rk_rt, "dm = R t", 0},
    {RK_R, B(RK_F)|B(RK_VEX), rk_r, "F = vex R", 0},
    {RK_VEX, B(RK_F)|B(RK_R), rk_vex, "F = vex R", 0},
    {RK_A, B(RK_F)|B(RK_M), rk_a, "F = m a", 0},
    {RK_M, B(RK_F)|B(RK_A), rk_m, "F = m a", 0},
    {RK_DM, B(RK_R)|B(RK_T), rk_dm, "dm = R t", 0},
    {RK_T, B(RK_DM)|B(RK_R), rk_t, "dm = R t", 0},
};

/* ---------- force that changes with time: F = A t^n ---------------- */
enum { PW_A, PW_N, PW_T1, PW_T2, PW_F2, PW_J, PW_M, PW_V1, PW_V2, PW_P1, PW_P2 };
#define PWN1 (v[PW_N] + 1)
FN(pw_j,   v[PW_A] * (pow(v[PW_T2], PWN1) - pow(v[PW_T1], PWN1)) / PWN1)
FN(pw_f2,  v[PW_A] * pow(v[PW_T2], v[PW_N]))
FN(pw_af,  v[PW_F2] / pow(v[PW_T2], v[PW_N]))
FN(pw_jv,  v[PW_M] * (v[PW_V2] - v[PW_V1]))
FN(pw_jp,  v[PW_P2] - v[PW_P1])
FN(pw_v2,  v[PW_V1] + v[PW_J] / v[PW_M])
FN(pw_v1,  v[PW_V2] - v[PW_J] / v[PW_M])
FN(pw_m,   v[PW_J] / (v[PW_V2] - v[PW_V1]))
FN(pw_p2,  v[PW_P1] + v[PW_J])
FN(pw_p1,  v[PW_P2] - v[PW_J])
FN(pw_t2,  pow(v[PW_J] * PWN1 / v[PW_A] + pow(v[PW_T1], PWN1), 1.0 / PWN1))
FN(pw_aj,  v[PW_J] * PWN1 / (pow(v[PW_T2], PWN1) - pow(v[PW_T1], PWN1)))
FN(pw_t2f, pow(v[PW_F2] / v[PW_A], 1.0 / v[PW_N]))
FN(pw_n,   fabs(log(v[PW_T2])) < 1e-6 ? UNK : log(v[PW_F2] / v[PW_A]) / log(v[PW_T2]))
static double pw_t1(const double *v) {
    double hi = pow(v[PW_T2], PWN1), x = hi - v[PW_J] * PWN1 / v[PW_A];
    if (x < 0 && -x <= 1e-6 * hi) x = 0;
    return x < 0 ? UNK : pow(x, 1.0 / PWN1);
}
static const var_t pow_vars[] = {
    {"A",  "",    "A: the number in F = A t^n", UNK},
    {"n",  "",    "n: the power of t (F = 3t: n = 1)", UNK},
    {"t1", "s",   "t1: start time (blank = 0)", 0},
    {"t2", "s",   "t2: end time", UNK},
    {"F2", "N",   "F2: size of the force at time t2", UNK},
    {"J",  "N s", "J: impulse = integral of F dt", UNK},
    {"m",  "kg",  "m: mass", UNK},
    {"v1", "m/s", "v1: velocity at t1 (blank = 0)", 0},
    {"v2", "m/s", "v2: velocity at t2", UNK},
    {"p1", "kg m/s", "p1: momentum at t1 (blank = 0)", 0},
    {"p2", "kg m/s", "p2: momentum at t2", UNK},
};
#define PW_INT "J = integral of A t^n dt, t1 to t2"
static const rule_t pow_rules[] = {
    {PW_A,  B(PW_F2)|B(PW_N)|B(PW_T2), pw_af, "F = A t^n", 0},
    {PW_J,  B(PW_A)|B(PW_N)|B(PW_T1)|B(PW_T2), pw_j, PW_INT, 0},
    {PW_J,  B(PW_M)|B(PW_V1)|B(PW_V2), pw_jv, "J = m (v2 - v1)", 0},
    {PW_J,  B(PW_P1)|B(PW_P2), pw_jp, "J = p2 - p1", 0},
    {PW_T2, B(PW_J)|B(PW_A)|B(PW_N)|B(PW_T1), pw_t2, PW_INT, R_ABS},
    {PW_A,  B(PW_J)|B(PW_N)|B(PW_T1)|B(PW_T2), pw_aj, PW_INT, 0},
    {PW_F2, B(PW_A)|B(PW_N)|B(PW_T2), pw_f2, "F = A t^n", 0},
    {PW_V2, B(PW_M)|B(PW_V1)|B(PW_J), pw_v2, "J = m (v2 - v1)", 0},
    {PW_P2, B(PW_P1)|B(PW_J), pw_p2, "J = p2 - p1", 0},
    {PW_M,  B(PW_J)|B(PW_V1)|B(PW_V2), pw_m, "J = m (v2 - v1)", 0},
    {PW_V1, B(PW_M)|B(PW_V2)|B(PW_J), pw_v1, "J = m (v2 - v1)", 0},
    {PW_P1, B(PW_P2)|B(PW_J), pw_p1, "J = p2 - p1", 0},
    {PW_T2, B(PW_F2)|B(PW_A)|B(PW_N), pw_t2f, "F = A t^n", R_ABS},
    {PW_T1, B(PW_J)|B(PW_A)|B(PW_N)|B(PW_T2), pw_t1, PW_INT, R_ABS},
    {PW_N,  B(PW_F2)|B(PW_A)|B(PW_T2), pw_n, "F = A t^n", 0},
};

/* ---------- rate of change: q = B t^n, dq/dt = n B t^(n-1) ---------- */
FN(dv_q,   v[0] * pow(v[2], v[1]))
FN(dv_r,   v[1] * v[0] * pow(v[2], v[1] - 1))
FN(dv_bq,  v[3] / pow(v[2], v[1]))
FN(dv_br,  v[4] / (v[1] * pow(v[2], v[1] - 1)))
FN(dv_t,   pow(v[3] / v[0], 1.0 / v[1]))
FN(dv_tr,  fabs(v[1] - 1) < 1e-6 ? UNK : pow(v[4] / (v[1] * v[0]), 1.0 / (v[1] - 1)))
FN(dv_n,   fabs(log(v[2])) < 1e-6 ? UNK : log(v[3] / v[0]) / log(v[2]))
static const var_t der_vars[] = {
    {"B",  "",  "B: the number in q = B t^n", UNK},
    {"n",  "",  "n: the power of t", UNK},
    {"t",  "s", "t: the time you care about", UNK},
    {"q",  "",  "q: value at t (p, x or v)", UNK},
    {"dq/dt", "", "rate at t: F from p, v from x, a from v", UNK},
};
static const rule_t der_rules[] = {
    {3, B(0)|B(1)|B(2), dv_q, "q = B t^n", 0},
    {4, B(0)|B(1)|B(2), dv_r, "dq/dt = n B t^(n-1)", 0},
    {0, B(3)|B(1)|B(2), dv_bq, "q = B t^n", 0},
    {0, B(4)|B(1)|B(2), dv_br, "dq/dt = n B t^(n-1)", 0},
    {2, B(3)|B(0)|B(1), dv_t, "q = B t^n", R_ABS},
    {2, B(4)|B(0)|B(1), dv_tr, "dq/dt = n B t^(n-1)", R_ABS},
    {1, B(3)|B(0)|B(2), dv_n, "q = B t^n", 0},
};

/* ---------- 16. torque and rotational dynamics --------------------- */
enum { Q_F, Q_R, Q_PHI, Q_TAU, Q_I, Q_AL, Q_W, Q_K, Q_L };
#define QSIN sin(v[Q_PHI] * DEG)
FN(q_tau,  v[Q_F] * v[Q_R] * QSIN)
FN(q_f,    v[Q_TAU] / (v[Q_R] * QSIN))
FN(q_r,    v[Q_TAU] / (v[Q_F] * QSIN))
FN(q_tau2, v[Q_I] * v[Q_AL])
FN(q_al,   v[Q_TAU] / v[Q_I])
FN(q_i,    v[Q_TAU] / v[Q_AL])
FN(q_k,    0.5 * v[Q_I] * SQ(v[Q_W]))
FN(q_wk,   root(2 * v[Q_K] / v[Q_I]))
FN(q_ik,   2 * v[Q_K] / SQ(v[Q_W]))
FN(q_l,    v[Q_I] * v[Q_W])
FN(q_wl,   v[Q_L] / v[Q_I])
FN(q_il,   v[Q_L] / v[Q_W])
FN(q_phi,  asinz(v[Q_TAU] / (v[Q_F] * v[Q_R])) / DEG)
static const var_t torq_vars[] = {
    {"F",   "N",       "F: force applied", UNK},
    {"r",   "m",       "r: distance from the pivot", UNK},
    {"phi", "deg",     "phi: angle between F and r", 90},
    {"tau", "N m",     "tau: torque", UNK},
    {"I",   "kg m^2",  "I: moment of inertia", UNK},
    {"alpha","rad/s^2","alpha: angular acceleration", UNK},
    {"w",   "rad/s",   "w: angular speed", UNK},
    {"Krot","J",       "Krot: rotational kinetic energy", UNK},
    {"L",   "kg m^2/s","L: angular momentum", UNK},
};
static const rule_t torq_rules[] = {
    {Q_TAU, B(Q_F)|B(Q_R)|B(Q_PHI), q_tau, "tau = F r sin(phi)", 0},
    {Q_TAU, B(Q_I)|B(Q_AL), q_tau2, "tau = I alpha", 0},
    {Q_AL, B(Q_TAU)|B(Q_I), q_al, "tau = I alpha", 0},
    {Q_I, B(Q_TAU)|B(Q_AL), q_i, "tau = I alpha", 0},
    {Q_K, B(Q_I)|B(Q_W), q_k, "Krot = (1/2) I w^2", 0},
    {Q_L, B(Q_I)|B(Q_W), q_l, "L = I w", 0},
    {Q_W, B(Q_L)|B(Q_I), q_wl, "L = I w", 0},
    {Q_I, B(Q_L)|B(Q_W), q_il, "L = I w", 0},
    {Q_W, B(Q_K)|B(Q_I), q_wk, "Krot = (1/2) I w^2", R_ABS},
    {Q_I, B(Q_K)|B(Q_W), q_ik, "Krot = (1/2) I w^2", 0},
    {Q_F, B(Q_TAU)|B(Q_R)|B(Q_PHI), q_f, "tau = F r sin(phi)", 0},
    {Q_R, B(Q_TAU)|B(Q_F)|B(Q_PHI), q_r, "tau = F r sin(phi)", 0},
    {Q_PHI, B(Q_TAU)|B(Q_F)|B(Q_R), q_phi, "tau = F r sin(phi)", R_NOCHK},
};

/* ---------- 17. moment of inertia of common shapes ----------------- */
FN(mi_hoop,  v[0] * SQ(v[1]))
FN(mi_disk,  0.5 * v[0] * SQ(v[1]))
FN(mi_sph,   0.4 * v[0] * SQ(v[1]))
FN(mi_shell, 2.0 / 3.0 * v[0] * SQ(v[1]))
FN(mi_rodc,  v[0] * SQ(v[2]) / 12.0)
FN(mi_rode,  v[0] * SQ(v[2]) / 3.0)
#define MI_INV(nm, idx, c, dim) \
    FN(mi_##nm##_d, root(v[idx] / ((c) * v[0]))) \
    FN(mi_##nm##_m, v[idx] / ((c) * SQ(v[dim])))
MI_INV(hoop, 3, 1.0, 1) MI_INV(disk, 4, 0.5, 1) MI_INV(sph, 5, 0.4, 1) MI_INV(shell, 6, 2.0 / 3.0, 1)
MI_INV(rodc, 7, 1.0 / 12.0, 2) MI_INV(rode, 8, 1.0 / 3.0, 2)
static const var_t moi_vars[] = {
    {"M",     "kg",     "M: mass of the object", UNK},
    {"R",     "m",      "R: radius (round shapes)", UNK},
    {"L",     "m",      "L: length (rods)", UNK},
    {"hoop",  "kg m^2", "hoop or thin ring: M R^2", UNK},
    {"disk",  "kg m^2", "solid disk/cylinder: (1/2) M R^2", UNK},
    {"sphere","kg m^2", "solid sphere: (2/5) M R^2", UNK},
    {"shell", "kg m^2", "hollow sphere: (2/3) M R^2", UNK},
    {"rod_c", "kg m^2", "rod about center: (1/12) M L^2", UNK},
    {"rod_e", "kg m^2", "rod about one end: (1/3) M L^2", UNK},
};
static const rule_t moi_rules[] = {
    {3, B(0)|B(1), mi_hoop, "I = M R^2", 0},
    {4, B(0)|B(1), mi_disk, "I = (1/2) M R^2", 0},
    {5, B(0)|B(1), mi_sph, "I = (2/5) M R^2", 0},
    {6, B(0)|B(1), mi_shell, "I = (2/3) M R^2", 0},
    {7, B(0)|B(2), mi_rodc, "I = (1/12) M L^2", 0},
    {8, B(0)|B(2), mi_rode, "I = (1/3) M L^2", 0},
#define MI_BACK(nm, idx, dim, eq) {dim, B(idx)|B(0), mi_##nm##_d, eq, R_ABS}, {0, B(idx)|B(dim), mi_##nm##_m, eq, 0}
    MI_BACK(hoop, 3, 1, "I = M R^2"), MI_BACK(disk, 4, 1, "I = (1/2) M R^2"),
    MI_BACK(sph, 5, 1, "I = (2/5) M R^2"), MI_BACK(shell, 6, 1, "I = (2/3) M R^2"),
    MI_BACK(rodc, 7, 2, "I = (1/12) M L^2"), MI_BACK(rode, 8, 2, "I = (1/3) M L^2"),
};

/* ---------- 18. gravitation ---------------------------------------- */
enum { G_M1, G_M2, G_R, G_F, G_G, G_VO, G_T, G_VE, G_U };
#define GM (BIGG * v[G_M1])
FN(gr_f,   (GM / v[G_R]) * (v[G_M2] / v[G_R]))
FN(gr_m2,  v[G_F] * v[G_R] / GM * v[G_R])
FN(gr_g,   GM / v[G_R] / v[G_R])
FN(gr_m1g, v[G_G] * v[G_R] / BIGG * v[G_R])
FN(gr_rg,  root(GM / v[G_G]))
FN(gr_vo,  root(GM / v[G_R]))
FN(gr_rvo, GM / SQ(v[G_VO]))
FN(gr_m1v, SQ(v[G_VO]) * v[G_R] / BIGG)
FN(gr_t,   2 * PI * v[G_R] / v[G_VO])
FN(gr_rt,  pow(GM / (4 * PI * PI) * SQ(v[G_T]), 1.0 / 3.0))
FN(gr_ve,  root(2 * GM / v[G_R]))
FN(gr_u,   -(GM / v[G_R]) * v[G_M2])
static const var_t grav_vars[] = {
    {"M",    "kg",    "M: planet (Earth: 5.97E24)", UNK},
    {"m",    "kg",    "m: small object or satellite", UNK},
    {"r",    "m",     "r: center to center (Earth R 6.37E6)", UNK},
    {"F",    "N",     "F: gravitational force", UNK},
    {"g",    "m/s^2", "g: gravity strength at r", UNK},
    {"v_orb","m/s",   "v_orb: speed of a circular orbit", UNK},
    {"T",    "s",     "T: time for one orbit", UNK},
    {"v_esc","m/s",   "v_esc: escape speed from r", UNK},
    {"U",    "J",     "U: gravitational PE = -G M m / r", UNK},
};
static const rule_t grav_rules[] = {
    {G_F, B(G_M1)|B(G_M2)|B(G_R), gr_f, "F = G M m / r^2", 0},
    {G_G, B(G_M1)|B(G_R), gr_g, "g = G M / r^2", 0},
    {G_VO, B(G_M1)|B(G_R), gr_vo, "v_orb = sqrt(G M / r)", 0},
    {G_VE, B(G_M1)|B(G_R), gr_ve, "v_esc = sqrt(2 G M / r)", 0},
    {G_T, B(G_R)|B(G_VO), gr_t, "T = 2 pi r / v_orb", 0},
    {G_R, B(G_M1)|B(G_T), gr_rt, "T^2 = 4 pi^2 r^3 / (G M)", 0},
    {G_R, B(G_M1)|B(G_G), gr_rg, "g = G M / r^2", 0},
    {G_R, B(G_M1)|B(G_VO), gr_rvo, "v_orb = sqrt(G M / r)", 0},
    {G_M1, B(G_G)|B(G_R), gr_m1g, "g = G M / r^2", 0},
    {G_M1, B(G_VO)|B(G_R), gr_m1v, "v_orb = sqrt(G M / r)", 0},
    {G_M2, B(G_F)|B(G_M1)|B(G_R), gr_m2, "F = G M m / r^2", 0},
    {G_U, B(G_M1)|B(G_M2)|B(G_R), gr_u, "U = -G M m / r", 0},
};

/* ---------- 19. springs and pendulums ------------------------------ */
enum { H_M, H_K, H_L, H_T, H_F, H_W, H_A, H_VM, H_AM, H_E, H_X, H_V };
FN(h_w,    root(v[H_K] / v[H_M]))
FN(h_k,    v[H_M] * SQ(v[H_W]))
FN(h_m,    v[H_K] / SQ(v[H_W]))
FN(h_t,    2 * PI / v[H_W])
FN(h_wt,   2 * PI / v[H_T])
FN(h_f,    1 / v[H_T])
FN(h_tf,   1 / v[H_F])
FN(h_vm,   v[H_W] * v[H_A])
FN(h_a,    v[H_VM] / v[H_W])
FN(h_am,   SQ(v[H_W]) * v[H_A])
FN(h_tl,   2 * PI * root(v[H_L] / GRAV))
FN(h_l,    GRAV * SQ(v[H_T] / (2 * PI)))
FN(h_e,    0.5 * v[H_K] * SQ(v[H_A]))
FN(h_ae,   root(2 * v[H_E] / v[H_K]))
FN(h_ke,   2 * v[H_E] / SQ(v[H_A]))
FN(h_vx,   v[H_W] * rootz(SQ(v[H_A]) - SQ(v[H_X]), SQ(v[H_A]) + SQ(v[H_X])))
FN(h_xv,   rootz(SQ(v[H_A]) - SQ(v[H_V] / v[H_W]), SQ(v[H_A]) + SQ(v[H_V] / v[H_W])))
static const var_t shm_vars[] = {
    {"m",    "kg",    "m: mass on the spring", UNK},
    {"k",    "N/m",   "k: spring constant", UNK},
    {"L",    "m",     "L: pendulum length with this period", UNK},
    {"T",    "s",     "T: period, time for one cycle", UNK},
    {"f",    "Hz",    "f: frequency, cycles per second", UNK},
    {"w",    "rad/s", "w: angular frequency", UNK},
    {"A",    "m",     "A: amplitude (max stretch)", UNK},
    {"vmax", "m/s",   "vmax: top speed (at the middle)", UNK},
    {"amax", "m/s^2", "amax: top acceleration (at the ends)", UNK},
    {"E",    "J",     "E: total energy = (1/2) k A^2", UNK},
    {"x",    "m",     "x: a position away from the middle", UNK},
    {"v",    "m/s",   "v: speed when it is at x", UNK},
};
static const rule_t shm_rules[] = {
    {H_W, B(H_K)|B(H_M), h_w, "w = sqrt(k / m)", 0},
    {H_T, B(H_L), h_tl, "T = 2 pi sqrt(L / g)", 0},
    {H_T, B(H_F), h_tf, "f = 1 / T", 0},
    {H_T, B(H_W), h_t, "T = 2 pi / w", 0},
    {H_W, B(H_T), h_wt, "T = 2 pi / w", 0},
    {H_F, B(H_T), h_f, "f = 1 / T", 0},
    {H_K, B(H_M)|B(H_W), h_k, "w = sqrt(k / m)", 0},
    {H_M, B(H_K)|B(H_W), h_m, "w = sqrt(k / m)", 0},
    {H_VM, B(H_W)|B(H_A), h_vm, "vmax = w A", 0},
    {H_AM, B(H_W)|B(H_A), h_am, "amax = w^2 A", 0},
    {H_A, B(H_VM)|B(H_W), h_a, "vmax = w A", 0},
    {H_L, B(H_T), h_l, "T = 2 pi sqrt(L / g)  (pendulum)", 0},
    {H_E, B(H_K)|B(H_A), h_e, "E = (1/2) k A^2", 0},
    {H_V, B(H_W)|B(H_A)|B(H_X), h_vx, "v = w sqrt(A^2 - x^2)", R_ABS},
    {H_A, B(H_E)|B(H_K), h_ae, "E = (1/2) k A^2", R_ABS},
    {H_K, B(H_E)|B(H_A), h_ke, "E = (1/2) k A^2", 0},
    {H_X, B(H_W)|B(H_A)|B(H_V), h_xv, "v = w sqrt(A^2 - x^2)", R_ABS},
};


/* =================== more problem types (chapters 1-14) ============ */
#define RULE(t, need, f, eq) {t, need, f, eq, 0}

/* ---------- add two vectors ---------------------------------------- */
enum { AV_A, AV_TA, AV_B, AV_TB, AV_RX, AV_RY, AV_R, AV_TR };
#define AVC(m, a) (v[m] * cos(v[a] * DEG))
#define AVS(m, a) (v[m] * sin(v[a] * DEG))
FN(av_rx,  AVC(AV_A, AV_TA) + AVC(AV_B, AV_TB))
FN(av_ry,  AVS(AV_A, AV_TA) + AVS(AV_B, AV_TB))
FN(av_r,   sqrt(SQ(v[AV_RX]) + SQ(v[AV_RY])))
FN(av_tr,  atan2(v[AV_RY], v[AV_RX]) / DEG)
FN(av_rxr, AVC(AV_R, AV_TR))
FN(av_ryr, AVS(AV_R, AV_TR))
FN(av_b,   sqrt(SQ(v[AV_RX] - AVC(AV_A, AV_TA)) + SQ(v[AV_RY] - AVS(AV_A, AV_TA))))
FN(av_tb,  atan2(v[AV_RY] - AVS(AV_A, AV_TA), v[AV_RX] - AVC(AV_A, AV_TA)) / DEG)
FN(av_a,   sqrt(SQ(v[AV_RX] - AVC(AV_B, AV_TB)) + SQ(v[AV_RY] - AVS(AV_B, AV_TB))))
FN(av_ta,  atan2(v[AV_RY] - AVS(AV_B, AV_TB), v[AV_RX] - AVC(AV_B, AV_TB)) / DEG)
static const var_t addv_vars[] = {
    {"A",  "",    "A: size of the first vector", UNK},
    {"thA","deg", "thA: its direction (east 0, north 90)", UNK},
    {"B",  "",    "B: size of the second vector", UNK},
    {"thB","deg", "thB: its direction", UNK},
    {"Rx", "",    "Rx: x-part of the sum", UNK},
    {"Ry", "",    "Ry: y-part of the sum", UNK},
    {"R",  "",    "R: size of the sum (resultant)", UNK},
    {"thR","deg", "thR: direction of the sum", UNK},
};
#define AV_AB (B(AV_A)|B(AV_TA)|B(AV_B)|B(AV_TB))
static const rule_t addv_rules[] = {
    RULE(AV_RX, AV_AB, av_rx, "Rx = A cos(thA) + B cos(thB)"),
    RULE(AV_RY, AV_AB, av_ry, "Ry = A sin(thA) + B sin(thB)"),
    RULE(AV_RX, B(AV_R)|B(AV_TR), av_rxr, "Rx = R cos(thR)"),
    RULE(AV_RY, B(AV_R)|B(AV_TR), av_ryr, "Ry = R sin(thR)"),
    RULE(AV_R,  B(AV_RX)|B(AV_RY), av_r, "R = sqrt(Rx^2 + Ry^2)"),
    {AV_TR, B(AV_RX)|B(AV_RY), av_tr, "tan(thR) = Ry / Rx", R_NOCHK},
    RULE(AV_B,  B(AV_A)|B(AV_TA)|B(AV_RX)|B(AV_RY), av_b, "B = R - A (x and y parts)"),
    {AV_TB, B(AV_A)|B(AV_TA)|B(AV_RX)|B(AV_RY), av_tb, "tan(thB) = By / Bx", R_NOCHK},
    RULE(AV_A,  B(AV_B)|B(AV_TB)|B(AV_RX)|B(AV_RY), av_a, "A = R - B (x and y parts)"),
    {AV_TA, B(AV_B)|B(AV_TB)|B(AV_RX)|B(AV_RY), av_ta, "tan(thA) = Ay / Ax", R_NOCHK},
};

/* ---------- unit conversions --------------------------------------- */
enum { CU_MS, CU_KMH, CU_MPH, CU_RAD, CU_DEG, CU_REV, CU_M, CU_FT, CU_MI, CU_N, CU_LB, CU_W, CU_HP };
#define CONV(nm, a, b, k) FN(cu_##nm##_f, v[a] * (k)) FN(cu_##nm##_b, v[b] / (k))
CONV(kmh, CU_MS, CU_KMH, 3.6) CONV(mph, CU_MS, CU_MPH, 2.236936)
CONV(deg, CU_RAD, CU_DEG, 180.0 / PI) CONV(rev, CU_RAD, CU_REV, 1.0 / (2 * PI))
CONV(ft, CU_M, CU_FT, 3.280840) CONV(mi, CU_M, CU_MI, 1.0 / 1609.344)
CONV(lb, CU_N, CU_LB, 1.0 / 4.448222) CONV(hp, CU_W, CU_HP, 1.0 / 746.0)
static const var_t conv_vars[] = {
    {"m/s",  "", "speed in meters per second", UNK},
    {"km/h", "", "speed in kilometers per hour", UNK},
    {"mi/h", "", "speed in miles per hour", UNK},
    {"rad",  "", "angle in radians", UNK},
    {"deg",  "", "angle in degrees", UNK},
    {"rev",  "", "angle in revolutions (turns)", UNK},
    {"m",    "", "length in meters", UNK},
    {"ft",   "", "length in feet", UNK},
    {"mi",   "", "length in miles", UNK},
    {"N",    "", "force in newtons", UNK},
    {"lb",   "", "force in pounds", UNK},
    {"W",    "", "power in watts", UNK},
    {"hp",   "", "power in horsepower", UNK},
};
#define CPAIR(nm, a, b, eq) RULE(b, B(a), cu_##nm##_f, eq), RULE(a, B(b), cu_##nm##_b, eq)
static const rule_t conv_rules[] = {
    CPAIR(kmh, CU_MS, CU_KMH, "1 m/s = 3.6 km/h"), CPAIR(mph, CU_MS, CU_MPH, "1 m/s = 2.237 mi/h"),
    CPAIR(deg, CU_RAD, CU_DEG, "pi rad = 180 deg"), CPAIR(rev, CU_RAD, CU_REV, "1 rev = 2 pi rad"),
    CPAIR(ft, CU_M, CU_FT, "1 m = 3.281 ft"), CPAIR(mi, CU_M, CU_MI, "1 mi = 1609 m"),
    CPAIR(lb, CU_N, CU_LB, "1 lb = 4.448 N"), CPAIR(hp, CU_W, CU_HP, "1 hp = 746 W"),
};

/* ---------- vertical circles and banked curves --------------------- */
enum { VC_M, VC_R, VC_V, VC_NT, VC_NB, VC_NH, VC_VMIN, VC_TH, VC_VB };
#define VCA (SQ(v[VC_V]) / v[VC_R])
FN(vc_nt,   v[VC_M] * (VCA - GRAV))
FN(vc_nb,   v[VC_M] * (VCA + GRAV))
FN(vc_nh,   v[VC_M] * (GRAV - VCA))
FN(vc_vmin, root(GRAV * v[VC_R]))
FN(vc_rmin, SQ(v[VC_VMIN]) / GRAV)
FN(vc_vb,   root(GRAV * v[VC_R] * tan(v[VC_TH] * DEG)))
FN(vc_thb,  atan(SQ(v[VC_VB]) / (GRAV * v[VC_R])) / DEG)
FN(vc_rb,   SQ(v[VC_VB]) / (GRAV * tan(v[VC_TH] * DEG)))
FN(vc_vnt,  root(v[VC_R] * (v[VC_NT] / v[VC_M] + GRAV)))
FN(vc_vnb,  root(v[VC_R] * (v[VC_NB] / v[VC_M] - GRAV)))
FN(vc_mnb,  v[VC_NB] / (VCA + GRAV))
FN(vc_mnt,  v[VC_NT] / (VCA - GRAV))
static const var_t vcirc_vars[] = {
    {"m",    "kg",  "m: mass", UNK},
    {"r",    "m",   "r: radius of the circle or curve", UNK},
    {"v",    "m/s", "v: speed at that point", UNK},
    {"N_top","N",   "N_top: push/tension at top (inside)", UNK},
    {"N_bot","N",   "N_bot: push/tension at the bottom", UNK},
    {"N_hill","N",  "N_hill: normal force on top of a hill", UNK},
    {"v_min","m/s", "v_min: slowest at top to stay on", UNK},
    {"th",   "deg", "th: banking angle of the road", UNK},
    {"v_bank","m/s","v_bank: speed needing no friction", UNK},
};
static const rule_t vcirc_rules[] = {
    RULE(VC_NT, B(VC_M)|B(VC_R)|B(VC_V), vc_nt, "N_top + m g = m v^2 / r"),
    RULE(VC_NB, B(VC_M)|B(VC_R)|B(VC_V), vc_nb, "N_bot - m g = m v^2 / r"),
    RULE(VC_NH, B(VC_M)|B(VC_R)|B(VC_V), vc_nh, "m g - N_hill = m v^2 / r"),
    RULE(VC_VMIN, B(VC_R), vc_vmin, "v_min = sqrt(g r)"),
    RULE(VC_VB, B(VC_R)|B(VC_TH), vc_vb, "tan(th) = v^2 / (g r)"),
    RULE(VC_R, B(VC_VMIN), vc_rmin, "v_min = sqrt(g r)"),
    RULE(VC_TH, B(VC_VB)|B(VC_R), vc_thb, "tan(th) = v^2 / (g r)"),
    RULE(VC_R, B(VC_VB)|B(VC_TH), vc_rb, "tan(th) = v^2 / (g r)"),
    {VC_V, B(VC_NB)|B(VC_M)|B(VC_R), vc_vnb, "N_bot - m g = m v^2 / r", R_ABS},
    {VC_V, B(VC_NT)|B(VC_M)|B(VC_R), vc_vnt, "N_top + m g = m v^2 / r", R_ABS},
    RULE(VC_M, B(VC_NB)|B(VC_V)|B(VC_R), vc_mnb, "N_bot - m g = m v^2 / r"),
    RULE(VC_M, B(VC_NT)|B(VC_V)|B(VC_R), vc_mnt, "N_top + m g = m v^2 / r"),
};

/* ---------- pulling or pushing at an angle, with friction ---------- */
enum { PA_M, PA_F, PA_TH, PA_MU, PA_N, PA_FK, PA_A };
#define PAC cos(v[PA_TH] * DEG)
#define PAS sin(v[PA_TH] * DEG)
FN(pa_n,   v[PA_M] * GRAV - v[PA_F] * PAS)
FN(pa_fk,  v[PA_MU] * v[PA_N])
FN(pa_a,   (v[PA_F] * PAC - v[PA_FK]) / v[PA_M])
FN(pa_f,   v[PA_M] * (v[PA_A] + v[PA_MU] * GRAV) / (PAC + v[PA_MU] * PAS))
FN(pa_mu,  v[PA_FK] / v[PA_N])
FN(pa_fka, v[PA_F] * PAC - v[PA_M] * v[PA_A])
FN(pa_m,   v[PA_F] * (PAC + v[PA_MU] * PAS) / (v[PA_A] + v[PA_MU] * GRAV))
FN(pa_thn, asinz((v[PA_M] * GRAV - v[PA_N]) / v[PA_F]) / DEG)
/* two angles work; the usual one first, the other if it fits the rest better */
FN(pa_th,  (atan(v[PA_MU]) + ((flip & B(PA_TH)) ? -1 : 1)
            * acosz(v[PA_M] * (v[PA_A] + v[PA_MU] * GRAV) / (v[PA_F] * sqrt(1 + SQ(v[PA_MU]))))) / DEG)
static const var_t pull_vars[] = {
    {"m",    "kg",    "m: mass", UNK},
    {"F",    "N",     "F: the pull (or push)", UNK},
    {"th",   "deg",   "th: pull angle above flat (push: -)", UNK},
    {"mu_k", "",      "mu_k: kinetic friction (none: 0)", UNK},
    {"N",    "N",     "N: normal force = m g - F sin(th)", UNK},
    {"fk",   "N",     "fk: friction force", UNK},
    {"a",    "m/s^2", "a: acceleration (steady speed: 0)", UNK},
};
static const rule_t pull_rules[] = {
    RULE(PA_N,  B(PA_M)|B(PA_F)|B(PA_TH), pa_n, "N = m g - F sin(th)"),
    RULE(PA_FK, B(PA_MU)|B(PA_N), pa_fk, "fk = mu_k N"),
    RULE(PA_A,  B(PA_F)|B(PA_TH)|B(PA_FK)|B(PA_M), pa_a, "F cos(th) - fk = m a"),
    RULE(PA_F,  B(PA_M)|B(PA_A)|B(PA_MU)|B(PA_TH), pa_f, "F cos(th) - mu_k (m g - F sin th) = m a"),
    RULE(PA_FK, B(PA_F)|B(PA_TH)|B(PA_M)|B(PA_A), pa_fka, "F cos(th) - fk = m a"),
    RULE(PA_MU, B(PA_FK)|B(PA_N), pa_mu, "fk = mu_k N"),
    RULE(PA_M,  B(PA_F)|B(PA_TH)|B(PA_MU)|B(PA_A), pa_m, "F cos(th) - mu_k (m g - F sin th) = m a"),
    RULE(PA_TH, B(PA_M)|B(PA_N)|B(PA_F), pa_thn, "N = m g - F sin(th)"),
    {PA_TH, B(PA_M)|B(PA_F)|B(PA_MU)|B(PA_A), pa_th, "F cos(th) - mu_k (m g - F sin th) = m a", R_NOCHK|R_ALT},
};

/* ---------- pushed up a ramp --------------------------------------- */
enum { IP_M, IP_TH, IP_MU, IP_F, IP_A };
#define IPG (GRAV * (sin(v[IP_TH] * DEG) + v[IP_MU] * cos(v[IP_TH] * DEG)))
FN(ip_a,   v[IP_F] / v[IP_M] - IPG)
FN(ip_f,   v[IP_M] * (v[IP_A] + IPG))
FN(ip_m,   v[IP_F] / (v[IP_A] + IPG))
FN(ip_th,  (asin((v[IP_F] / v[IP_M] - v[IP_A]) / GRAV / sqrt(1 + SQ(v[IP_MU]))) - atan(v[IP_MU])) / DEG)
FN(ip_mu,  (v[IP_F] / v[IP_M] - v[IP_A] - GRAV * sin(v[IP_TH] * DEG)) / (GRAV * cos(v[IP_TH] * DEG)))
static const var_t ipush_vars[] = {
    {"m",    "kg",    "m: mass of the block", UNK},
    {"th",   "deg",   "th: angle of the ramp", UNK},
    {"mu_k", "",      "mu_k: kinetic friction (none: 0)", UNK},
    {"F",    "N",     "F: push or pull up along the ramp", UNK},
    {"a",    "m/s^2", "a: accel up the ramp (steady: 0)", UNK},
};
#define IP_EQ "F - m g sin th - mu_k m g cos th = m a"
static const rule_t ipush_rules[] = {
    RULE(IP_A, B(IP_M)|B(IP_TH)|B(IP_MU)|B(IP_F), ip_a, IP_EQ),
    RULE(IP_F, B(IP_M)|B(IP_TH)|B(IP_MU)|B(IP_A), ip_f, IP_EQ),
    RULE(IP_M, B(IP_F)|B(IP_TH)|B(IP_MU)|B(IP_A), ip_m, IP_EQ),
    RULE(IP_MU, B(IP_M)|B(IP_TH)|B(IP_F)|B(IP_A), ip_mu, IP_EQ),
    RULE(IP_TH, B(IP_M)|B(IP_MU)|B(IP_F)|B(IP_A), ip_th, IP_EQ),
};

/* ---------- elevator / apparent weight ----------------------------- */
FN(ev_n,   v[0] * (GRAV + v[1]))
FN(ev_w,   v[0] * GRAV)
FN(ev_a,   v[2] / v[0] - GRAV)
FN(ev_m,   v[2] / (GRAV + v[1]))
FN(ev_mw,  v[3] / GRAV)
static const var_t elev_vars[] = {
    {"m", "kg",    "m: mass", UNK},
    {"a", "m/s^2", "a: accel (up +, down -, free fall -9.8)", UNK},
    {"N", "N",     "N: scale reading or cable tension", UNK},
    {"w", "N",     "w: true weight = m g", UNK},
};
static const rule_t elev_rules[] = {
    RULE(3, B(0), ev_w, "w = m g"),
    RULE(0, B(3), ev_mw, "w = m g"),
    RULE(2, B(0)|B(1), ev_n, "N - m g = m a"),
    RULE(1, B(2)|B(0), ev_a, "N - m g = m a"),
    RULE(0, B(2)|B(1), ev_m, "N - m g = m a"),
};

/* ---------- weight hanging from two ropes -------------------------- */
enum { HR_W, HR_M, HR_T1, HR_A1, HR_T2, HR_A2 };
#define HRS sin((v[HR_A1] + v[HR_A2]) * DEG)
FN(hr_w,   v[HR_M] * GRAV)
FN(hr_m,   v[HR_W] / GRAV)
FN(hr_t1,  v[HR_W] * cos(v[HR_A2] * DEG) / HRS)
FN(hr_t2,  v[HR_W] * cos(v[HR_A1] * DEG) / HRS)
FN(hr_t2h, v[HR_T1] * cos(v[HR_A1] * DEG) / cos(v[HR_A2] * DEG))
FN(hr_t1h, v[HR_T2] * cos(v[HR_A2] * DEG) / cos(v[HR_A1] * DEG))
FN(hr_w1,  v[HR_T1] * HRS / cos(v[HR_A2] * DEG))
FN(hr_a2,  acosz(v[HR_T1] * cos(v[HR_A1] * DEG) / v[HR_T2]) / DEG)
FN(hr_a1,  acosz(v[HR_T2] * cos(v[HR_A2] * DEG) / v[HR_T1]) / DEG)
static const var_t rope_vars[] = {
    {"w",  "N",   "w: weight being held up", UNK},
    {"m",  "kg",  "m: its mass", UNK},
    {"T1", "N",   "T1: tension in rope 1", UNK},
    {"a1", "deg", "a1: rope 1 angle above horizontal", UNK},
    {"T2", "N",   "T2: tension in rope 2", UNK},
    {"a2", "deg", "a2: rope 2 angle (flat rope: 0)", UNK},
};
static const rule_t rope_rules[] = {
    RULE(HR_W,  B(HR_M), hr_w, "w = m g"),
    RULE(HR_T1, B(HR_W)|B(HR_A1)|B(HR_A2), hr_t1, "T1 = w cos(a2) / sin(a1 + a2)"),
    RULE(HR_T2, B(HR_W)|B(HR_A1)|B(HR_A2), hr_t2, "T2 = w cos(a1) / sin(a1 + a2)"),
    RULE(HR_T2, B(HR_T1)|B(HR_A1)|B(HR_A2), hr_t2h, "T1 cos(a1) = T2 cos(a2)"),
    RULE(HR_T1, B(HR_T2)|B(HR_A1)|B(HR_A2), hr_t1h, "T1 cos(a1) = T2 cos(a2)"),
    RULE(HR_W,  B(HR_T1)|B(HR_A1)|B(HR_A2), hr_w1, "T1 sin(a1) + T2 sin(a2) = w"),
    RULE(HR_M,  B(HR_W), hr_m, "w = m g"),
    RULE(HR_A2, B(HR_T1)|B(HR_A1)|B(HR_T2), hr_a2, "T1 cos(a1) = T2 cos(a2)"),
    RULE(HR_A1, B(HR_T2)|B(HR_A2)|B(HR_T1), hr_a1, "T1 cos(a1) = T2 cos(a2)"),
};

/* ---------- work-energy theorem ------------------------------------ */
enum { WE_M, WE_V1, WE_V2, WE_W, WE_F, WE_D };
#define WEK (2 * v[WE_W] / v[WE_M])
FN(we_w,   0.5 * v[WE_M] * (SQ(v[WE_V2]) - SQ(v[WE_V1])))
FN(we_wf,  v[WE_F] * v[WE_D])
FN(we_f,   v[WE_W] / v[WE_D])
FN(we_d,   v[WE_W] / v[WE_F])
FN(we_v2,  rootz(SQ(v[WE_V1]) + WEK, SQ(v[WE_V1]) + fabs(WEK)))
FN(we_v1,  rootz(SQ(v[WE_V2]) - WEK, SQ(v[WE_V2]) + fabs(WEK)))
FN(we_m,   2 * v[WE_W] / (SQ(v[WE_V2]) - SQ(v[WE_V1])))
static const var_t wet_vars[] = {
    {"m",  "kg",  "m: mass", UNK},
    {"v1", "m/s", "v1: speed at the start (rest: 0)", UNK},
    {"v2", "m/s", "v2: speed at the end (stops: 0)", UNK},
    {"W",  "J",   "W: total work done on it", UNK},
    {"F",  "N",   "F: net force along motion (brake: -)", UNK},
    {"d",  "m",   "d: distance it acts over", UNK},
};
static const rule_t wet_rules[] = {
    RULE(WE_W,  B(WE_F)|B(WE_D), we_wf, "W = F d"),
    RULE(WE_W,  B(WE_M)|B(WE_V1)|B(WE_V2), we_w, "W = (1/2) m v2^2 - (1/2) m v1^2"),
    {WE_V2, B(WE_M)|B(WE_V1)|B(WE_W), we_v2, "W = (1/2) m v2^2 - (1/2) m v1^2", R_ABS},
    {WE_V1, B(WE_M)|B(WE_V2)|B(WE_W), we_v1, "W = (1/2) m v2^2 - (1/2) m v1^2", R_ABS},
    RULE(WE_F,  B(WE_W)|B(WE_D), we_f, "W = F d"),
    RULE(WE_D,  B(WE_W)|B(WE_F), we_d, "W = F d"),
    RULE(WE_M,  B(WE_W)|B(WE_V1)|B(WE_V2), we_m, "W = (1/2) m v2^2 - (1/2) m v1^2"),
};

/* ---------- springs: Hooke's law ----------------------------------- */
enum { HK_K, HK_X, HK_F, HK_U, HK_M };
FN(hk_f,   v[HK_K] * v[HK_X])
FN(hk_fm,  v[HK_M] * GRAV)
FN(hk_k,   v[HK_F] / v[HK_X])
FN(hk_x,   v[HK_F] / v[HK_K])
FN(hk_u,   0.5 * v[HK_K] * SQ(v[HK_X]))
FN(hk_m,   v[HK_F] / GRAV)
FN(hk_xu,  root(2 * v[HK_U] / v[HK_K]))
FN(hk_ku,  2 * v[HK_U] / SQ(v[HK_X]))
static const var_t hooke_vars[] = {
    {"k", "N/m", "k: spring constant", UNK},
    {"x", "m",   "x: stretch or squeeze (cm: E-2)", UNK},
    {"F", "N",   "F: force on / from the spring", UNK},
    {"U", "J",   "U: energy stored = (1/2) k x^2", UNK},
    {"m", "kg",  "m: mass hanging at rest on it", UNK},
};
static const rule_t hooke_rules[] = {
    RULE(HK_F, B(HK_M), hk_fm, "F = m g (hanging at rest)"),
    RULE(HK_F, B(HK_K)|B(HK_X), hk_f, "F = k x"),
    RULE(HK_K, B(HK_F)|B(HK_X), hk_k, "F = k x"),
    RULE(HK_X, B(HK_F)|B(HK_K), hk_x, "F = k x"),
    RULE(HK_U, B(HK_K)|B(HK_X), hk_u, "U = (1/2) k x^2"),
    {HK_X, B(HK_U)|B(HK_K), hk_xu, "U = (1/2) k x^2", R_ABS},
    RULE(HK_K, B(HK_U)|B(HK_X), hk_ku, "U = (1/2) k x^2"),
    RULE(HK_M, B(HK_F), hk_m, "F = m g (hanging at rest)"),
};

/* ---------- linear <-> angular ------------------------------------- */
enum { LA_R, LA_W, LA_V, LA_AL, LA_AT, LA_AR, LA_TH, LA_S };
FN(la_v,   v[LA_R] * v[LA_W])
FN(la_w,   v[LA_V] / v[LA_R])
FN(la_rv,  v[LA_V] / v[LA_W])
FN(la_at,  v[LA_R] * v[LA_AL])
FN(la_al,  v[LA_AT] / v[LA_R])
FN(la_ra,  v[LA_AT] / v[LA_AL])
FN(la_ar,  SQ(v[LA_W]) * v[LA_R])
FN(la_war, root(v[LA_AR] / v[LA_R]))
FN(la_rar, v[LA_AR] / SQ(v[LA_W]))
FN(la_s,   v[LA_R] * v[LA_TH])
FN(la_th,  v[LA_S] / v[LA_R])
FN(la_rs,  v[LA_S] / v[LA_TH])
static const var_t linang_vars[] = {
    {"r",    "m",       "r: distance from the axis", UNK},
    {"w",    "rad/s",   "w: angular speed", UNK},
    {"v",    "m/s",     "v: speed of that point", UNK},
    {"alpha","rad/s^2", "alpha: angular acceleration", UNK},
    {"a_tan","m/s^2",   "a_tan: accel along the circle", UNK},
    {"a_rad","m/s^2",   "a_rad: accel toward the center", UNK},
    {"th",   "rad",     "th: angle turned (rad)", UNK},
    {"s",    "m",       "s: distance along the arc", UNK},
};
static const rule_t linang_rules[] = {
    RULE(LA_V, B(LA_R)|B(LA_W), la_v, "v = r w"),
    RULE(LA_W, B(LA_V)|B(LA_R), la_w, "v = r w"),
    RULE(LA_R, B(LA_V)|B(LA_W), la_rv, "v = r w"),
    RULE(LA_AT, B(LA_R)|B(LA_AL), la_at, "a_tan = r alpha"),
    RULE(LA_AL, B(LA_AT)|B(LA_R), la_al, "a_tan = r alpha"),
    RULE(LA_R, B(LA_AT)|B(LA_AL), la_ra, "a_tan = r alpha"),
    RULE(LA_AR, B(LA_W)|B(LA_R), la_ar, "a_rad = w^2 r = v^2 / r"),
    RULE(LA_S, B(LA_R)|B(LA_TH), la_s, "s = r th"),
    RULE(LA_TH, B(LA_S)|B(LA_R), la_th, "s = r th"),
    RULE(LA_R, B(LA_S)|B(LA_TH), la_rs, "s = r th"),
    {LA_W, B(LA_AR)|B(LA_R), la_war, "a_rad = w^2 r", R_ABS},
    RULE(LA_R, B(LA_AR)|B(LA_W), la_rar, "a_rad = w^2 r"),
};

/* ---------- parallel-axis theorem ---------------------------------- */
FN(px_i,   v[0] + v[1] * SQ(v[2]))
FN(px_icm, v[3] - v[1] * SQ(v[2]))
FN(px_m,   (v[3] - v[0]) / SQ(v[2]))
FN(px_d,   root((v[3] - v[0]) / v[1]))
static const var_t pax_vars[] = {
    {"Icm", "kg m^2", "Icm: I about the center of mass", UNK},
    {"M",   "kg",     "M: mass of the object", UNK},
    {"d",   "m",      "d: distance between the two axes", UNK},
    {"I",   "kg m^2", "I: I about the new, parallel axis", UNK},
};
static const rule_t pax_rules[] = {
    RULE(3, B(0)|B(1)|B(2), px_i, "I = Icm + M d^2"),
    RULE(0, B(3)|B(1)|B(2), px_icm, "I = Icm + M d^2"),
    RULE(1, B(3)|B(0)|B(2), px_m, "I = Icm + M d^2"),
    {2, B(3)|B(0)|B(1), px_d, "I = Icm + M d^2", R_ABS},
};

/* ---------- rolling without slipping ------------------------------- */
enum { RL_C, RL_H, RL_V, RL_TH, RL_A, RL_M, RL_KT, RL_KR };
FN(rl_v,   root(2 * GRAV * v[RL_H] / (1 + v[RL_C])))
FN(rl_h,   SQ(v[RL_V]) * (1 + v[RL_C]) / (2 * GRAV))
FN(rl_c,   2 * GRAV * v[RL_H] / SQ(v[RL_V]) - 1)
FN(rl_a,   GRAV * sin(v[RL_TH] * DEG) / (1 + v[RL_C]))
FN(rl_th,  asin(v[RL_A] * (1 + v[RL_C]) / GRAV) / DEG)
FN(rl_ca,  GRAV * sin(v[RL_TH] * DEG) / v[RL_A] - 1)
FN(rl_kt,  0.5 * v[RL_M] * SQ(v[RL_V]))
FN(rl_kr,  0.5 * v[RL_C] * v[RL_M] * SQ(v[RL_V]))
FN(rl_m,   2 * v[RL_KT] / SQ(v[RL_V]))
static const var_t roll_vars[] = {
    {"c",   "",      "c: shape number (see Tips; disk 0.5)", UNK},
    {"h",   "m",     "h: height it rolls down, from rest", UNK},
    {"v",   "m/s",   "v: speed at the bottom", UNK},
    {"th",  "deg",   "th: angle of the ramp", UNK},
    {"a",   "m/s^2", "a: acceleration down the ramp", UNK},
    {"m",   "kg",    "m: mass (only for the energies)", UNK},
    {"Ktr", "J",     "Ktr: moving (translational) KE", UNK},
    {"Krot","J",     "Krot: spinning (rotational) KE", UNK},
};
static const rule_t roll_rules[] = {
    RULE(RL_V, B(RL_H)|B(RL_C), rl_v, "m g h = (1/2)(1 + c) m v^2"),
    RULE(RL_A, B(RL_TH)|B(RL_C), rl_a, "a = g sin(th) / (1 + c)"),
    RULE(RL_H, B(RL_V)|B(RL_C), rl_h, "m g h = (1/2)(1 + c) m v^2"),
    RULE(RL_C, B(RL_H)|B(RL_V), rl_c, "m g h = (1/2)(1 + c) m v^2"),
    RULE(RL_C, B(RL_TH)|B(RL_A), rl_ca, "a = g sin(th) / (1 + c)"),
    RULE(RL_TH, B(RL_A)|B(RL_C), rl_th, "a = g sin(th) / (1 + c)"),
    RULE(RL_KT, B(RL_M)|B(RL_V), rl_kt, "Ktr = (1/2) m v^2"),
    RULE(RL_KR, B(RL_C)|B(RL_M)|B(RL_V), rl_kr, "Krot = (1/2) I w^2 = (1/2) c m v^2"),
    RULE(RL_M, B(RL_KT)|B(RL_V), rl_m, "Ktr = (1/2) m v^2"),
};

/* ---------- conservation of angular momentum ----------------------- */
FN(am_w2,  v[0] * v[1] / v[2])
FN(am_i2,  v[0] * v[1] / v[3])
FN(am_w1,  v[2] * v[3] / v[0])
FN(am_i1,  v[2] * v[3] / v[1])
FN(am_k1,  0.5 * v[0] * SQ(v[1]))
FN(am_k2,  0.5 * v[2] * SQ(v[3]))
static const var_t angm_vars[] = {
    {"I1", "kg m^2", "I1: moment of inertia before", UNK},
    {"w1", "rad/s",  "w1: angular speed before", UNK},
    {"I2", "kg m^2", "I2: moment of inertia after", UNK},
    {"w2", "rad/s",  "w2: angular speed after", UNK},
    {"K1", "J",      "K1: rotational KE before", UNK},
    {"K2", "J",      "K2: rotational KE after", UNK},
};
static const rule_t angm_rules[] = {
    RULE(3, B(0)|B(1)|B(2), am_w2, "I1 w1 = I2 w2"),
    RULE(2, B(0)|B(1)|B(3), am_i2, "I1 w1 = I2 w2"),
    RULE(1, B(0)|B(2)|B(3), am_w1, "I1 w1 = I2 w2"),
    RULE(0, B(1)|B(2)|B(3), am_i1, "I1 w1 = I2 w2"),
    RULE(4, B(0)|B(1), am_k1, "K = (1/2) I w^2"),
    RULE(5, B(2)|B(3), am_k2, "K = (1/2) I w^2"),
};

/* ---------- lever / seesaw balance --------------------------------- */
FN(lv_f2,  v[0] * v[1] / v[3])
FN(lv_d2,  v[0] * v[1] / v[2])
FN(lv_f1,  v[2] * v[3] / v[1])
FN(lv_d1,  v[2] * v[3] / v[0])
static const var_t lever_vars[] = {
    {"F1", "N", "F1: force (or weight m g) on one side", UNK},
    {"d1", "m", "d1: its distance from the pivot", UNK},
    {"F2", "N", "F2: force on the other side", UNK},
    {"d2", "m", "d2: its distance from the pivot", UNK},
};
static const rule_t lever_rules[] = {
    RULE(2, B(0)|B(1)|B(3), lv_f2, "F1 d1 = F2 d2 (torques balance)"),
    RULE(3, B(0)|B(1)|B(2), lv_d2, "F1 d1 = F2 d2 (torques balance)"),
    RULE(0, B(1)|B(2)|B(3), lv_f1, "F1 d1 = F2 d2 (torques balance)"),
    RULE(1, B(0)|B(2)|B(3), lv_d1, "F1 d1 = F2 d2 (torques balance)"),
};

/* ---------- physical pendulum -------------------------------------- */
enum { PP_I, PP_M, PP_D, PP_T, PP_F };
#define PPT2 SQ(v[PP_T] / (2 * PI))
FN(pp_t,   2 * PI * root(v[PP_I] / (v[PP_M] * GRAV * v[PP_D])))
FN(pp_f,   1 / v[PP_T])
FN(pp_tf,  1 / v[PP_F])
FN(pp_i,   v[PP_M] * GRAV * v[PP_D] * PPT2)
FN(pp_d,   v[PP_I] / (v[PP_M] * GRAV * PPT2))
FN(pp_m,   v[PP_I] / (GRAV * v[PP_D] * PPT2))
static const var_t ppend_vars[] = {
    {"I", "kg m^2", "I: moment of inertia about the pivot", UNK},
    {"m", "kg",     "m: mass", UNK},
    {"d", "m",      "d: pivot to center of mass", UNK},
    {"T", "s",      "T: period of the swing", UNK},
    {"f", "Hz",     "f: frequency", UNK},
};
static const rule_t ppend_rules[] = {
    RULE(PP_T, B(PP_I)|B(PP_M)|B(PP_D), pp_t, "T = 2 pi sqrt(I / (m g d))"),
    RULE(PP_T, B(PP_F), pp_tf, "f = 1 / T"),
    RULE(PP_F, B(PP_T), pp_f, "f = 1 / T"),
    RULE(PP_I, B(PP_M)|B(PP_D)|B(PP_T), pp_i, "T = 2 pi sqrt(I / (m g d))"),
    RULE(PP_D, B(PP_I)|B(PP_M)|B(PP_T), pp_d, "T = 2 pi sqrt(I / (m g d))"),
    RULE(PP_M, B(PP_I)|B(PP_D)|B(PP_T), pp_m, "T = 2 pi sqrt(I / (m g d))"),
};

/* ---------- density and pressure ----------------------------------- */
enum { DP_RHO, DP_M, DP_V, DP_H, DP_P0, DP_PG, DP_P, DP_A, DP_F };
FN(dp_rho, v[DP_M] / v[DP_V])
FN(dp_m,   v[DP_RHO] * v[DP_V])
FN(dp_v,   v[DP_M] / v[DP_RHO])
FN(dp_pg,  v[DP_RHO] * GRAV * v[DP_H])
FN(dp_h,   v[DP_PG] / (v[DP_RHO] * GRAV))
FN(dp_rhp, v[DP_PG] / (GRAV * v[DP_H]))
FN(dp_p,   v[DP_P0] + v[DP_PG])
FN(dp_pgp, v[DP_P] - v[DP_P0])
FN(dp_f,   v[DP_P] * v[DP_A])
FN(dp_a,   v[DP_F] / v[DP_P])
FN(dp_pf,  v[DP_F] / v[DP_A])
static const var_t dens_vars[] = {
    {"rho", "kg/m^3", "rho: density (water 1000)", UNK},
    {"m",   "kg",     "m: mass", UNK},
    {"V",   "m^3",    "V: volume (1 L = 1E-3, 1 cm^3 = 1E-6)", UNK},
    {"h",   "m",      "h: depth below the surface", UNK},
    {"p0",  "Pa",     "p0: pressure on top (blank = 1 atm)", 1.013e5},
    {"p_g", "Pa",     "p_g: gauge pressure = rho g h", UNK},
    {"p",   "Pa",     "p: absolute pressure = p0 + rho g h", UNK},
    {"A",   "m^2",    "A: area the pressure pushes on", UNK},
    {"F",   "N",      "F: force = p A", UNK},
};
static const rule_t dens_rules[] = {
    RULE(DP_RHO, B(DP_M)|B(DP_V), dp_rho, "rho = m / V"),
    RULE(DP_M, B(DP_RHO)|B(DP_V), dp_m, "rho = m / V"),
    RULE(DP_V, B(DP_M)|B(DP_RHO), dp_v, "rho = m / V"),
    RULE(DP_PG, B(DP_RHO)|B(DP_H), dp_pg, "p_g = rho g h"),
    RULE(DP_PG, B(DP_P)|B(DP_P0), dp_pgp, "p = p0 + p_g"),
    RULE(DP_P, B(DP_F)|B(DP_A), dp_pf, "F = p A"),
    RULE(DP_P, B(DP_P0)|B(DP_PG), dp_p, "p = p0 + p_g"),
    RULE(DP_H, B(DP_PG)|B(DP_RHO), dp_h, "p_g = rho g h"),
    RULE(DP_RHO, B(DP_PG)|B(DP_H), dp_rhp, "p_g = rho g h"),
    RULE(DP_F, B(DP_P)|B(DP_A), dp_f, "F = p A"),
    RULE(DP_A, B(DP_F)|B(DP_P), dp_a, "F = p A"),
};

/* ---------- buoyancy ----------------------------------------------- */
enum { BY_RF, BY_V, BY_FB, BY_RO, BY_FR, BY_M, BY_WA };
FN(by_fb,  v[BY_RF] * GRAV * v[BY_V])
FN(by_v,   v[BY_FB] / (v[BY_RF] * GRAV))
FN(by_rf,  v[BY_FB] / (GRAV * v[BY_V]))
FN(by_fr,  v[BY_RO] / v[BY_RF])
FN(by_ro,  v[BY_FR] * v[BY_RF])
FN(by_rff, v[BY_RO] / v[BY_FR])
FN(by_wa,  v[BY_M] * GRAV - v[BY_FB])
FN(by_m,   (v[BY_WA] + v[BY_FB]) / GRAV)
FN(by_fbw, v[BY_M] * GRAV - v[BY_WA])
static const var_t buoy_vars[] = {
    {"rho_f","kg/m^3","rho_f: density of the fluid", UNK},
    {"V_sub","m^3",   "V_sub: volume that is under the fluid", UNK},
    {"F_b",  "N",     "F_b: buoyant force = rho_f g V_sub", UNK},
    {"rho_o","kg/m^3","rho_o: density of the object", UNK},
    {"frac", "",      "frac: part under when it floats", UNK},
    {"m",    "kg",    "m: mass of the object", UNK},
    {"w_app","N",     "w_app: apparent weight = m g - F_b", UNK},
};
static const rule_t buoy_rules[] = {
    RULE(BY_FB, B(BY_RF)|B(BY_V), by_fb, "F_b = rho_f g V_sub"),
    RULE(BY_FB, B(BY_M)|B(BY_WA), by_fbw, "w_app = m g - F_b"),
    RULE(BY_FR, B(BY_RO)|B(BY_RF), by_fr, "frac under = rho_o / rho_f"),
    RULE(BY_WA, B(BY_M)|B(BY_FB), by_wa, "w_app = m g - F_b"),
    RULE(BY_V, B(BY_FB)|B(BY_RF), by_v, "F_b = rho_f g V_sub"),
    RULE(BY_RF, B(BY_FB)|B(BY_V), by_rf, "F_b = rho_f g V_sub"),
    RULE(BY_RO, B(BY_FR)|B(BY_RF), by_ro, "frac under = rho_o / rho_f"),
    RULE(BY_RF, B(BY_RO)|B(BY_FR), by_rff, "frac under = rho_o / rho_f"),
    RULE(BY_M, B(BY_WA)|B(BY_FB), by_m, "w_app = m g - F_b"),
};

/* ---------- fluid flow: continuity + Bernoulli --------------------- */
enum { FL_A1, FL_V1, FL_A2, FL_V2, FL_Q, FL_RHO, FL_P1, FL_Y1, FL_P2, FL_Y2 };
#define FLDP (2 * (v[FL_P1] - v[FL_P2]) / v[FL_RHO] + 2 * GRAV * (v[FL_Y1] - v[FL_Y2]))
#define FLSC (2 * (fabs(v[FL_P1]) + fabs(v[FL_P2])) / v[FL_RHO] + 2 * GRAV * (fabs(v[FL_Y1]) + fabs(v[FL_Y2])))
FN(fl_q1,  v[FL_A1] * v[FL_V1])
FN(fl_q2,  v[FL_A2] * v[FL_V2])
FN(fl_v1q, v[FL_Q] / v[FL_A1])
FN(fl_v2q, v[FL_Q] / v[FL_A2])
FN(fl_a1,  v[FL_Q] / v[FL_V1])
FN(fl_a2,  v[FL_Q] / v[FL_V2])
FN(fl_p2,  v[FL_P1] + 0.5 * v[FL_RHO] * (SQ(v[FL_V1]) - SQ(v[FL_V2])) + v[FL_RHO] * GRAV * (v[FL_Y1] - v[FL_Y2]))
FN(fl_p1,  v[FL_P2] - 0.5 * v[FL_RHO] * (SQ(v[FL_V1]) - SQ(v[FL_V2])) - v[FL_RHO] * GRAV * (v[FL_Y1] - v[FL_Y2]))
FN(fl_v2,  rootz(SQ(v[FL_V1]) + FLDP, SQ(v[FL_V1]) + FLSC))
FN(fl_v1,  rootz(SQ(v[FL_V2]) - FLDP, SQ(v[FL_V2]) + FLSC))
FN(fl_y2,  v[FL_Y1] + (v[FL_P1] - v[FL_P2]) / (v[FL_RHO] * GRAV) + (SQ(v[FL_V1]) - SQ(v[FL_V2])) / (2 * GRAV))
FN(fl_y1,  v[FL_Y2] - (v[FL_P1] - v[FL_P2]) / (v[FL_RHO] * GRAV) - (SQ(v[FL_V1]) - SQ(v[FL_V2])) / (2 * GRAV))
FN(fl_rho, (v[FL_P2] - v[FL_P1]) / (0.5 * (SQ(v[FL_V1]) - SQ(v[FL_V2])) + GRAV * (v[FL_Y1] - v[FL_Y2])))
static const var_t flow_vars[] = {
    {"A1",  "m^2",   "A1: pipe area at point 1", UNK},
    {"v1",  "m/s",   "v1: flow speed at point 1", UNK},
    {"A2",  "m^2",   "A2: pipe area at point 2", UNK},
    {"v2",  "m/s",   "v2: flow speed at point 2", UNK},
    {"Q",   "m^3/s", "Q: volume flow rate = A v", UNK},
    {"rho", "kg/m^3","rho: fluid density (blank = water)", 1000},
    {"p1",  "Pa",    "p1: pressure at point 1", UNK},
    {"y1",  "m",     "y1: height of point 1 (blank = 0)", 0},
    {"p2",  "Pa",    "p2: pressure at point 2", UNK},
    {"y2",  "m",     "y2: height of point 2 (blank = 0)", 0},
};
#define FL_BERN "p + (1/2) rho v^2 + rho g y = same"
#define FL_ALL(x) ((B(FL_V1)|B(FL_V2)|B(FL_RHO)|B(FL_P1)|B(FL_Y1)|B(FL_P2)|B(FL_Y2)) & ~B(x))
static const rule_t flow_rules[] = {
    RULE(FL_Q,  B(FL_A1)|B(FL_V1), fl_q1, "Q = A1 v1 = A2 v2"),
    RULE(FL_Q,  B(FL_A2)|B(FL_V2), fl_q2, "Q = A1 v1 = A2 v2"),
    RULE(FL_V1, B(FL_Q)|B(FL_A1), fl_v1q, "Q = A1 v1 = A2 v2"),
    RULE(FL_V2, B(FL_Q)|B(FL_A2), fl_v2q, "Q = A1 v1 = A2 v2"),
    RULE(FL_A1, B(FL_Q)|B(FL_V1), fl_a1, "Q = A1 v1 = A2 v2"),
    RULE(FL_A2, B(FL_Q)|B(FL_V2), fl_a2, "Q = A1 v1 = A2 v2"),
    RULE(FL_P2, FL_ALL(FL_P2), fl_p2, FL_BERN),
    RULE(FL_P1, FL_ALL(FL_P1), fl_p1, FL_BERN),
    {FL_V2, FL_ALL(FL_V2), fl_v2, FL_BERN, R_ABS},
    {FL_V1, FL_ALL(FL_V1), fl_v1, FL_BERN, R_ABS},
    RULE(FL_Y2, FL_ALL(FL_Y2), fl_y2, FL_BERN),
    RULE(FL_Y1, FL_ALL(FL_Y1), fl_y1, FL_BERN),
    RULE(FL_RHO, FL_ALL(FL_RHO), fl_rho, FL_BERN),
};

/* ---------- stress, strain, Young's modulus ------------------------ */
enum { YS_F, YS_A, YS_L0, YS_DL, YS_SS, YS_SN, YS_Y };
FN(ys_ss,  v[YS_F] / v[YS_A])
FN(ys_f,   v[YS_SS] * v[YS_A])
FN(ys_a,   v[YS_F] / v[YS_SS])
FN(ys_sn,  v[YS_DL] / v[YS_L0])
FN(ys_dl,  v[YS_SN] * v[YS_L0])
FN(ys_l0,  v[YS_DL] / v[YS_SN])
FN(ys_y,   v[YS_SS] / v[YS_SN])
FN(ys_ssy, v[YS_Y] * v[YS_SN])
FN(ys_sny, v[YS_SS] / v[YS_Y])
static const var_t young_vars[] = {
    {"F",   "N",   "F: stretching (or squeezing) force", UNK},
    {"A",   "m^2", "A: cross-section area (wire: pi r^2)", UNK},
    {"L0",  "m",   "L0: original length", UNK},
    {"dL",  "m",   "dL: change in length", UNK},
    {"stress","Pa","stress = F / A", UNK},
    {"strain","",  "strain = dL / L0 (no units)", UNK},
    {"Y",   "Pa",  "Y: Young's modulus (steel 2.0E11)", UNK},
};
static const rule_t young_rules[] = {
    RULE(YS_SS, B(YS_F)|B(YS_A), ys_ss, "stress = F / A"),
    RULE(YS_SN, B(YS_DL)|B(YS_L0), ys_sn, "strain = dL / L0"),
    RULE(YS_Y,  B(YS_SS)|B(YS_SN), ys_y, "Y = stress / strain"),
    RULE(YS_SS, B(YS_Y)|B(YS_SN), ys_ssy, "Y = stress / strain"),
    RULE(YS_SN, B(YS_SS)|B(YS_Y), ys_sny, "Y = stress / strain"),
    RULE(YS_F,  B(YS_SS)|B(YS_A), ys_f, "stress = F / A"),
    RULE(YS_A,  B(YS_F)|B(YS_SS), ys_a, "stress = F / A"),
    RULE(YS_DL, B(YS_SN)|B(YS_L0), ys_dl, "strain = dL / L0"),
    RULE(YS_L0, B(YS_DL)|B(YS_SN), ys_l0, "strain = dL / L0"),
};

/* ---------- topic table -------------------------------------------- */
#define N(a) ((uint8_t)(sizeof(a) / sizeof((a)[0])))
#define TOPIC(name, vars, rules, tips) { name, vars, N(vars), rules, N(rules), tips }
static const topic_t topics[] = {
#if PART == 0 || PART == 1
    /* Motion */
    TOPIC("Straight-line motion", lin_vars, lin_rules,
        "starts from rest: v0 = 0\ncomes to a stop: v = 0\ndropped: v0 = 0, a = -9.8\n"
        "thrown straight up: a = -9.8\nat the highest point: v = 0\n"
        "falls 20 m: dx = -20\nconstant speed: a = 0"),
    TOPIC("Projectile", proj_vars, proj_rules,
        "thrown horizontally: th = 0\nlands 20 m below: y = -20\n"
        "at the top: vy = 0\nH and R: just read them off\n"
        "R only works if it lands at the\nsame height it was launched from.\n"
        "Off a cliff: enter y, read x and t."),
    TOPIC("Circular motion", circ_vars, circ_rules,
        "one lap every 4 s: T = 4\nspins at 300 rpm: rpm = 300\n"
        "string tension (flat circle): Fc\ncar on a flat curve: mu_s, r, vmax\n"
        "diameter given: r = half of it"),
    TOPIC("Vertical circle + banked", vcirc_vars, vcirc_rules,
        "bucket, loop, ball on a string:\n  top: N_top   bottom: N_bot\n"
        "car over the top of a hill: N_hill\nfeels weightless: N = 0, so v = v_min\n"
        "banked road, no friction: th, r\n  give v_bank"),
    TOPIC("Vector components", vec_vars, vec_rules,
        "30 deg above +x: th = 30\n30 deg north of west: th = 150\n"
        "straight down: th = 270 (or -90)\nEnter Ax and Ay to get size + angle."),
    TOPIC("Add two vectors", addv_vars, addv_rules,
        "east 0, north 90, west 180, south 270\n30 north of east: 30\n30 west of north: 120\n"
        "boat 5 m/s north, current 3 m/s east:\n  A = 5, thA = 90, B = 3, thB = 0\n"
        "plane in wind: same idea\nsubtract: add 180 to thB\n"
        "Know the result and one vector:\n  type R, thR, A, thA to get B"),
    TOPIC("Convert units", conv_vars, conv_rules,
        "Type a number next to its unit;\nthe matching units fill in.\n"
        "rpm to rad/s: times 0.1047\ncm to m: E-2    mm to m: E-3\n"
        "g to kg: E-3    km to m: E3\nL to m^3: E-3   cm^2 to m^2: E-4"),
    /* Forces */
    TOPIC("F = ma + friction", force_vars, force_rules,
        "frictionless / smooth: mu_k = 0\nconstant speed: a = 0\n"
        "pulled with 50 N: F = 50\nweighs 98 N: w = 98\n"
        "For a = from motion, use\nStraight-line motion first."),
    TOPIC("Pull at an angle", pull_vars, pull_rules,
        "rope pulls 30 above horizontal: 30\npushing DOWN at 30: th = -30\n"
        "constant speed: a = 0\nfrictionless: mu_k = 0\n"
        "F needed to keep it moving:\n  type m, th, mu_k, a = 0"),
    TOPIC("Incline (ramp)", incl_vars, incl_rules,
        "frictionless ramp: mu_k = 0\nslides down: read a\n"
        "pushed up and slows: read a_up\njust about to slip: read th_s\n"
        "mass is only needed for forces."),
    TOPIC("Incline + push", ipush_vars, ipush_rules,
        "Block moving UP the ramp with a\npush or rope along the ramp.\n"
        "steady speed: a = 0\nfrictionless: mu_k = 0\n"
        "sliding with no push: use\nIncline (ramp) instead"),
    TOPIC("Atwood (2 hanging)", atw_vars, atw_rules,
        "two masses over a pulley,\nboth hanging. m2 is the heavy one.\n"
        "massless string + pulley assumed."),
    TOPIC("Table + hanging mass", tab_vars, tab_rules,
        "m1 slides on the table,\nm2 hangs off the edge.\n"
        "frictionless table: mu_k = 0"),
    TOPIC("Elevator (apparent wt)", elev_vars, elev_rules,
        "speeding up going up: a is +\nslowing down going up: a is -\n"
        "speeding up going down: a is -\nconstant speed: a = 0\n"
        "cable tension lifting a load:\n  same thing, read N\ncable snaps: a = -9.8, N = 0"),
    TOPIC("Hanging from 2 ropes", rope_vars, rope_rules,
        "Angles are measured up from\nhorizontal.\none rope horizontal: its angle = 0\n"
        "angle given from the ceiling: same\nangle given from vertical: 90 - it\n"
        "traffic light, picture on wires"),
    /* Energy */
    TOPIC("Work and power", work_vars, work_rules,
        "force along the motion: phi = 0\n(blank phi counts as 0)\n"
        "friction: phi = 180 (negative W)\nlifting at steady speed: F = m g\n"
        "2 hp motor: P = 1492"),
    TOPIC("Work-energy theorem", wet_vars, wet_rules,
        "from rest: v1 = 0\ncomes to a stop: v2 = 0\n"
        "braking / friction: F is negative\nstopping distance: type F, read d\n"
        "friction force = mu_k m g on a flat\nroad (type it as -mu_k m g)"),
    TOPIC("Energy conservation", en_vars, en_rules,
        "released from rest: v1 = 0\nreaches the ground: h2 = 0\n"
        "highest point: v2 = 0\nno friction: leave W blank\n"
        "friction does 30 J: W = -30\nno spring: leave x1, x2, k blank\n"
        "mass is not needed unless there\nis a spring or W."),
    TOPIC("Spring (Hooke's law)", hooke_vars, hooke_rules,
        "stretches 5 cm: x = 0.05\nmass hung on it: type m, read F\n"
        "  then x or k\n1.85 N/cm: k = 185\nenergy stored: read U"),
#endif
#if PART == 0 || PART == 2
    /* Momentum */
    TOPIC("Collision: stick together", stick_vars, stick_rules,
        "at rest / stationary: v = 0\nmoving left / west: negative v\n"
        "lock together, stick, embed, catch,\ndropped onto a moving cart: here\n"
        "dropped straight down: its v = 0\n4.2 g bullet: m = 4.2E-3\n"
        "weighs 625 N: m = 625/9.8 = 63.8\n"
        "spring bumpers, max energy stored:\n  it is -dK, at speed vf"),
    TOPIC("Collision: elastic", elas_vars, elas_rules,
        "elastic / bounce with no energy\nlost: use this one\n"
        "at rest: v = 0\nmoving left: negative v\n"
        "head-on = 1-D, signs matter\n"
        "masses as a ratio (deuteron =\n  2 x neutron): type 1 and 2\n"
        "speed as a fraction: type v1 = 1\n"
        "unknown mass: enter both speeds\n  of m1, leave m2 blank"),
    TOPIC("Recoil / any collision", gen_vars, gen_rules,
        "Fill in 5 of the 6, get the 6th.\nthrow, push off, recoil, expel:\n"
        "  both start at rest: v1 = v2 = 0\nbounces back: negative after\n"
        "east / right = +, west / left = -\ndK shows energy lost or created\n"
        "weighs 625 N: m = 625/9.8 = 63.8\nexpels 1.65 kg: m1 = what is left\n"
        "fired at an angle on ice: use\n  v cos(angle) for its v"),
    TOPIC("Impulse", imp_vars, imp_rules,
        "starts at rest: v1 = 0\nbounces back: v2 is negative\n"
        "brought to rest: v2 = 0\nin contact 2.4 ms: t = 2.4E-3\n"
        "force to the left: negative F\nforce-time graph: J = area under\n"
        "  it, then type J\naverage force: read F\n"
        "dropped and bounces: type h_drop\n  and h_up, not v1 and v2\n"
        "2-D with Fx, Fy given: do x,\n  then do y, one at a time"),
    TOPIC("Momentum and KE", mk_vars, mk_rules,
        "One object: enter any two of\nm, v, p, K to get the others.\n"
        "weighs 625 N: m = 625/9.8 = 63.8"),
    TOPIC("Stick, then what happens", then_vars, then_rules,
        "Something sticks in a block that\nwas at rest, THEN the pair:\n"
        "  swings up (ballistic pendulum),\n  or climbs a hill: h\n"
        "  slides to a stop: mu_k and d\n  squeezes a spring: k and x\n"
        "  flies off a table: Ht and xr\n"
        "Type what you know from either\nend; it works backward too.\n"
        "5 g bullet: m1 = 5E-3"),
    TOPIC("Spring between 2 blocks", spr2_vars, spr2_rules,
        "Blocks pushed apart from rest.\nv1 and v2 are speeds (opposite\ndirections).\n"
        "1.85 N/cm: k = 185\nsqueezed 20 cm: x = 0.2\n"
        "energy stored in spring: read U"),
    TOPIC("Center of mass", com_vars, com_rules,
        "Put one object at 0 (x1 blank).\n40 m ahead: x2 = 40\n"
        "stacked: use heights of centers\n"
        "same density, sizes given:\n  use m = diameter^3\n"
        "Sun 1.99E30, Jupiter 1.90E27,\n  apart 7.78E11"),
    TOPIC("Rocket thrust", rock_vars, rock_rules,
        "burns 0.054 kg each second: R\ngas speed relative to rocket: vex\n"
        "astronaut + pack: m = both masses\ngas used in 5 s: t = 5, read dm"),
    /* Momentum in 2-D */
    TOPIC("Impulse in 2-D", imp2_vars, imp2_rules,
        "Angles go counterclockwise from\nthe +x (right) direction.\n"
        "moving right: 0 (or leave blank)\nmoving left: 180\n"
        "left, 30 above horizontal: 150\nleft, 30 below horizontal: 210\n"
        "Read Fx and Fy for components."),
    TOPIC("Collision in 2-D", col2_vars, col2_rules,
        "m1 comes in along +x (a1 blank).\nm2 at rest: leave v2, a2 blank.\n"
        "Angles after: above the original\nline = +, below = -\n"
        "Know one object after: get the\nother's speed and direction.\n"
        "Know both angles after: get both\nspeeds.\n"
        "equal masses, not given: m = 1"),
    TOPIC("Stick together in 2-D", stk2_vars, stk2_rules,
        "Angles: east 0, north 90,\nwest 180, south 270 (or -90)\n"
        "60 south of east: -60\ntackle, enmeshed, move as one,\nlocked together: use this one\n"
        "Know the angle after but not one\nspeed before: type af, leave that\nspeed blank."),
    /* Changing force (calc) */
    TOPIC("Force F = A t^n", pow_vars, pow_rules,
        "F = 3.00 t:  A = 3, n = 1\nF = A t^2: n = 2; type F2 and t2\n  to get A\n"
        "impulse from 2 s to 3.5 s:\n  t1 = 2, t2 = 3.5, read J\n"
        "change in velocity: leave v1\n  blank, read v2\n"
        "two terms or i and j parts: do\n  each one, add the J's\n"
        "Other formulas: in FlowCE press\nMATH for the integral template,\nthen type the answer as J in\nImpulse."),
    TOPIC("Rate of change (d/dt)", der_vars, der_rules,
        "p = B t^n  ->  F = dp/dt\nx = B t^n  ->  v = dx/dt\nv = B t^n  ->  a = dv/dt\n"
        "constant term: its rate is 0\nB t (n = 1): rate is just B\n"
        "two terms: do each, then add\n"
        "Other formulas: FlowCE's MATH\nkey has a derivative template."),
#endif
#if PART == 0 || PART == 3
    /* Rotation */
    TOPIC("Rotational motion", rot_vars, rot_rules,
        "starts from rest: w0 = 0\nstops: w = 0\n"
        "5 revolutions: dth = 31.42\nrpm to rad/s: times 0.1047"),
    TOPIC("Linear <-> angular", linang_vars, linang_rules,
        "point on the rim: r = radius\nrpm to rad/s: times 0.1047\n"
        "1 rev = 6.283 rad\nrolling wheel: v of the car = r w\n"
        "total accel = sqrt(a_tan^2 +\n  a_rad^2)"),
    TOPIC("Torque and rotation", torq_vars, torq_rules,
        "force at a right angle: leave\nphi blank (counts as 90)\n"
        "need I? use Moment of inertia\nwrench 0.3 m long: r = 0.3"),
    TOPIC("Moment of inertia", moi_vars, moi_rules,
        "Enter M and R (or L for a rod).\nEvery shape's I shows up;\npick the one that matches."),
    TOPIC("Parallel-axis theorem", pax_vars, pax_rules,
        "Get Icm from Moment of inertia.\nrod about its end: d = L / 2\n"
        "disk about its rim: d = R"),
    TOPIC("Rolling", roll_vars, roll_rules,
        "c is the number in I = c M R^2:\n  hoop / ring         1\n  solid disk/cylinder 0.5\n"
        "  solid sphere        0.4\n  hollow sphere       0.667\n"
        "sliding, no rolling: c = 0\nstarts from rest at height h"),
    TOPIC("Angular momentum", angm_vars, angm_rules,
        "skater pulls arms in: I drops,\n  w goes up\n"
        "something dropped on a turntable:\n  I2 = I1 + I of the new thing\n"
        "point mass at radius r: I = m r^2\nrpm is fine for w if you only\n  want w (not the K's)"),
    TOPIC("Lever balance", lever_vars, lever_rules,
        "seesaw, wrench, balanced beam\na weight: F = m g\n"
        "the beam's own weight acts at\n  its middle\nthree forces: Ask AI"),
    /* Gravity + springs */
    TOPIC("Gravitation", grav_vars, grav_rules,
        "Earth: M = 5.97E24, R = 6.37E6\nMoon: M = 7.35E22, R = 1.74E6\n"
        "Sun: M = 1.99E30\n300 km up: r = 6.37E6 + 3E5\n   = 6.67E6\n^ key types the E"),
    TOPIC("Springs and pendulums", shm_vars, shm_rules,
        "spring: enter m and k\npendulum: enter L only\n"
        "10 cycles in 5 s: T = 0.5\npulled back 4 cm: A = 0.04"),
    TOPIC("Physical pendulum", ppend_vars, ppend_rules,
        "A swinging object that is not a\npoint mass (rod, hoop, sign).\n"
        "rod pivoted at its end:\n  I = (1/3) m L^2, d = L / 2\n"
        "use Parallel-axis for other pivots"),
    /* Fluids + solids */
    TOPIC("Density + pressure", dens_vars, dens_rules,
        "water 1000, seawater 1030,\nmercury 13600, air 1.2 (kg/m^3)\n"
        "1 atm = 1.013E5 Pa\n1 L = 1E-3 m^3, 1 cm^3 = 1E-6 m^3\n"
        "pressure 10 m down in a lake:\n  rho = 1000, h = 10, read p\n"
        "sealed tank: type its own p0"),
    TOPIC("Buoyancy", buoy_vars, buoy_rules,
        "fully under: V_sub = whole volume\n  = m / rho_o\n"
        "floating: F_b = m g, and frac\n  tells how much is under\n"
        "ice in water: 917 / 1000\nscale reading under water: w_app"),
    TOPIC("Flow (Bernoulli)", flow_vars, flow_rules,
        "pipe narrows: A1 v1 = A2 v2\nround pipe: A = pi r^2\n"
        "level pipe: leave y1, y2 blank\nwater: leave rho blank\n"
        "hole in a tank, h below surface:\n  v1 = 0, y1 = h, p1 = p2 (type\n  the same number), read v2"),
    TOPIC("Stress and strain", young_vars, young_rules,
        "wire of radius r: A = pi r^2\n  (2 mm diameter: A = 3.14E-6)\n"
        "steel 2.0E11, copper 1.1E11,\naluminum 7.0E10 (Pa)\n"
        "stretch 1 mm: dL = 1E-3"),
#endif
};
#define NUM_TOPICS N(topics)

typedef struct { const char *name; uint8_t count; } cat_t;
static const cat_t cats[] = {
#if PART == 0 || PART == 1
    {"Motion", 7},
    {"Forces", 8},
    {"Energy", 4},
#endif
#if PART == 0 || PART == 2
    {"Momentum", 9},
    {"Momentum in 2-D", 3},
    {"Changing force (calc)", 2},
#endif
#if PART == 0 || PART == 3
    {"Rotation", 8},
    {"Gravity + springs", 3},
    {"Fluids + solids", 4},
#endif
};
#define NUM_CATS N(cats)
static uint8_t cat_first(uint8_t c) {
    uint8_t f = 0;
    for (uint8_t i = 0; i < c; i++) f = (uint8_t)(f + cats[i].count);
    return f;
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
static uint8_t st[MAXV];
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
            if (fabs(x) < 1e-9) { pick = r; px = 0; break; }
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
        tol = 5e-3 * (fabs(x) > fabs(y) ? fabs(x) : fabs(y)) + 1e-6;
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
    tput("g = 9.80 m/s^2   G = 6.674E-11\nangles in degrees unless it says rad\n");
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
#if PART == 0 || PART == 1
    if (tp->vars == lin_vars) {
        if (st[3] == S_SOLVED && val[3] < 0) return "t < 0: check signs of v0, a, dx";
    } else if (tp->vars == proj_vars) {
        if (st[P_T] == S_SOLVED && tp->rules[by[P_T]].fn == p_ty && KN(P_Y) && val[P_Y] > 0)
            return "t is the 2nd time (on the way down)";
        if (st[P_R] == S_SOLVED && KN(P_Y) && fabs(val[P_Y]) > 1e-6) return "R only if it lands at launch height";
    } else if (tp->vars == circ_vars) {
        if (KN(C_V) && KN(C_VMAX) && val[C_V] > val[C_VMAX]) return "v > vmax: too fast, it skids";
    } else if (tp->vars == force_vars) {
        if (KN(F_A) && KN(F_APP) && KN(F_FK) && val[F_A] < 0) return "a < 0: friction beats the push";
    } else if (tp->vars == incl_vars) {
        if (KN(I_A) && val[I_A] < 0) return "a < 0: friction holds it in place";
        if (KN(I_TH) && KN(I_THS) && val[I_TH] < val[I_THS]) return "th < th_s: it will not start to slip";
    } else if (tp->vars == atw_vars) {
        if (KN(2) && val[2] < 0) return "a < 0: m1 is the one going down";
    } else if (tp->vars == tab_vars) {
        if (KN(3) && val[3] < 0) return "a < 0: friction holds, nothing moves";
    } else if (tp->vars == vcirc_vars) {
        if (KN(VC_NT) && val[VC_NT] < 0) return "N_top < 0: too slow, it falls away";
    } else if (tp->vars == pull_vars) {
        if (KN(PA_N) && val[PA_N] < 0) return "N < 0: the pull lifts it off the floor";
    } else if (tp->vars == elev_vars) {
        if (KN(2) && val[2] < 0) return "N < 0: falling faster than free fall?";
    }
#endif
#if PART == 0 || PART == 2
    if (tp->vars == stick_vars) {
        if (KN(5) && val[5] > 1e-6) return "dK > 0 cannot happen: check signs";
    } else if (tp->vars == gen_vars) {
        if (KN(6)) return val[6] > 1e-6 ? "KE went up: push-off or explosion"
                        : val[6] < -1e-6 ? "KE was lost: not an elastic collision" : "KE unchanged: elastic";
    } else if (tp->vars == col2_vars) {
        if (KN(T_DK) && val[T_DK] > 1e-6) return "dK > 0: check speeds and angles";
    }
#endif
#if PART == 0 || PART == 3
    if (tp->vars == shm_vars) {
        if (KN(H_X) && KN(H_A) && fabs(val[H_X]) > fabs(val[H_A])) return "x is beyond the amplitude A";
    } else if (tp->vars == buoy_vars) {
        if (KN(BY_FR) && val[BY_FR] > 1) return "frac > 1: denser than the fluid, sinks";
    }
#endif
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
        bool good = st[i] != S_UNK && fabs(val[i] - want) <= 2e-3 * fabs(want) + 1e-6;
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

int main(int argc, char **argv) {
    char b[32];
    if (argc > 1 && !strcmp(argv[1], "--solve")) return solve_lines();
    if (argc > 1 && !strcmp(argv[1], "--topics")) return list_topics();
    if (argc > 1 && !strcmp(argv[1], "--sheet")) { tp = &topics[0]; build_formula_sheet(); printf("%s\n[%u chars]\n", text, text_len); return 0; }
    /* number formatting */
    struct { double x; const char *s; } f[] = { {0, "0"}, {5.1, "5.1"}, {-19.799, "-19.8"}, {35.3477, "35.35"},
        {123456, "123456"}, {1234567, "1.235E6"}, {0.000123456, "1.235E-4"}, {0.05, "0.05"}, {9.99996, "10"},
        {5.97e24, "5.97E24"}, {1000, "1000"}, {0.001, "0.001"}, {-0.5, "-0.5"}, {99999.6, "100000"}, {6.674e-11, "6.674E-11"} };
    for (unsigned i = 0; i < sizeof f / sizeof f[0]; i++) {
        fmt_num(f[i].x, b);
        if (strcmp(b, f[i].s)) { fails++; printf("FMT %g -> '%s' want '%s'\n", f[i].x, b, f[i].s); }
    }
    if (fabs(parse_num("-2.5E-3") + 0.0025) > 1e-12 || parse_num("3E8") != 3e8 || ok(parse_num("-")) || ok(parse_num("1E")) || parse_num(".5") != 0.5) { fails++; printf("PARSE fail\n"); }

    /* every rule consistent (all variables known) */
    full("Straight-line motion", "v0=3 a=2 t=4");
    full("Rotational motion", "w0=1.5 alpha=0.8 t=6");
    full("Projectile", "v0=20 th=35 t=1.2");
    full("Circular motion", "r=2 T=3 m=1.5 mu_s=0.6");
    full("Vector components", "A=10 th=30");
    full("F = ma + friction", "m=4 mu_k=0.25 F=30");
    full("Incline (ramp)", "m=3 th=35 mu_k=0.2 mu_s=0.5");
    full("Atwood (2 hanging)", "m1=2 m2=3");
    full("Table + hanging mass", "m1=4 m2=3 mu_k=0.3");
    full("Work and power", "F=50 d=8 phi=25 t=4");
    full("Energy conservation", "m=2 v1=3 h1=5 h2=1 W=-12 x1=0.2 x2=0.05 k=400");
    full("Collision: stick together", "m1=2 v1=5 m2=3 v2=-1");
    full("Collision: elastic", "m1=2 v1=5 m2=3 v2=-1");
    full("Recoil / any collision", "m1=2 v1=5 m2=3 v2=-1 v1f=-0.5");
    full("Impulse", "m=0.5 h_drop=2 h_up=1.2 t=0.004");
    full("Momentum and KE", "m=3 v=4");
    full("Impulse in 2-D", "m=0.2 v1=10 a1=20 v2=14 a2=160 t=0.01");
    full("Collision in 2-D", "m1=2 v1=6 a1=10 m2=3 v2=1.5 a2=200 v1f=3 b1=70");
    full("Spring between 2 blocks", "m1=2 m2=3 k=500 x=0.1");
    full("Stick together in 2-D", "m1=2 v1=6 a1=10 m2=3 v2=1.5 a2=200");
    full("Stick, then what happens", "m1=0.5 v1=24 m2=8 mu_k=0.3 k=200 Ht=2.2 L=1.4");
    full("Center of mass", "m1=3 x1=1 m2=5 x2=9 v1=2 v2=-4");
    full("Rocket thrust", "vex=1500 R=0.05 m=40 t=6");
    full("Force F = A t^n", "A=4 n=2 t1=1 t2=3 m=5 v1=2 p1=-3");
    full("Rate of change (d/dt)", "B=-0.75 n=2 t=3");
    full("Torque and rotation", "F=20 r=0.4 phi=60 I=0.5 w=6");
    full("Moment of inertia", "M=3 R=0.2 L=1.5");
    full("Gravitation", "M=5.97E24 m=1000 r=6.67E6");
    full("Springs and pendulums", "m=0.5 k=200 A=0.04 x=0.01");

    /* word problems with known answers */
    check("Straight-line motion", "v0=0 a=-9.8 dx=-20 | t=2.020 v=-19.80");                 /* dropped 20 m */
    check("Straight-line motion", "v0=25 v=0 a=-5 | dx=62.5 t=5");                         /* braking car */
    check("Straight-line motion", "v0=15 a=-9.8 v=0 | dx=11.48 t=1.531");                  /* thrown up, top */
    check("Straight-line motion", "v0=0 t=6 dx=90 | a=5 v=30");
    check("Projectile", "v0=20 th=30 | H=5.102 R=35.35 v0x=17.32 v0y=10");
    check("Projectile", "v0=12 th=0 y=-20 | t=2.020 x=24.24 vy=-19.80 v=23.15");           /* off a cliff */
    check("Projectile", "th=45 R=40 | v0=19.80");
    check("Projectile", "v0x=15 x=45 v0y=0 | t=3 y=-44.1");
    check("Circular motion", "r=50 mu_s=0.8 | vmax=19.80");
    check("Circular motion", "rpm=300 r=0.1 | w=31.42 v=3.142 ac=98.70");
    check("F = ma + friction", "m=10 F=50 mu_k=0.2 | fk=19.6 Fnet=30.4 a=3.04");
    check("F = ma + friction", "m=10 F=50 a=3.04 | mu_k=0.2");
    check("Incline (ramp)", "th=30 mu_k=0 | a=4.9");
    check("Incline (ramp)", "th=30 a=2.5 | mu_k=0.2828");
    check("Incline (ramp)", "mu_s=0.5 | th_s=26.57");
    check("Atwood (2 hanging)", "m1=2 m2=3 | a=1.96 T=23.52");
    check("Table + hanging mass", "m1=4 m2=2 mu_k=0 | a=3.267 T=13.07");
    check("Work and power", "F=100 d=5 | W=500");                                           /* phi blank = 0 */
    check("Work and power", "F=100 d=5 phi=60 t=10 | W=250 P=25");
    check("Energy conservation", "v1=0 h1=20 h2=0 | v2=19.80");                             /* no mass needed */
    check("Energy conservation", "v1=10 h1=0 v2=0 | h2=5.102");
    check("Energy conservation", "m=2 v1=3 h1=0 v2=0 h2=0 x2=0.1 | k=1800");                /* block hits spring */
    check("Energy conservation", "m=2 v1=0 h1=5 v2=8 h2=0 | W=-34");                        /* work by friction */
    check("Energy conservation", "m=0.5 k=200 x1=0.1 v1=0 h1=0 h2=0 | v2=2");               /* spring launch */
    check("Energy conservation", "m=2 k=1800 v1=3 h1=0 v2=0 h2=0 | x2=0.1");                /* how far it compresses */
    check("Energy conservation", "m=0.5 k=200 v1=0 h1=0 v2=2 h2=0 x2=0 | x1=0.1");          /* how far was it pulled back */
    check("Energy conservation", "v1=0 h1=20 v2=19.799 h2=0 | W=0");
    check("Collision: stick together", "m1=1500 v1=20 m2=1000 v2=0 | vf=12 dK=-120000");
    check("Collision: stick together", "m1=0.005 m2=2 v2=0 vf=1.5 | v1=601.5");             /* bullet in block */
    check("Collision: elastic", "m1=2 v1=3 m2=1 v2=0 | v1f=1 v2f=4");
    check("Collision: elastic", "m1=1 v1=4 m2=1 v2=-2 | v1f=-2 v2f=4");
    check("Recoil / any collision", "m1=4 v1=0 m2=0.01 v2=0 v2f=300 | v1f=-0.75");          /* rifle recoil */
    check("Recoil / any collision", "m1=0.2 v1=10 m2=0.3 v2=0 v1f=-2 | v2f=8");
    check("Impulse", "m=0.145 v1=-40 v2=50 t=0.002 | J=13.05 F=6525");
    check("Impulse", "m=0.06 v1=0 F=300 t=0.005 | v2=25");
    check("Rotational motion", "w0=0 alpha=2 t=5 | w=10 dth=25");
    check("Torque and rotation", "F=40 r=0.3 I=2 | tau=12 alpha=6");                        /* phi blank = 90 */
    check("Moment of inertia", "M=2 R=0.5 | disk=0.25 hoop=0.5 sphere=0.2");
    check("Gravitation", "M=5.97E24 r=6.37E6 | g=9.820 v_esc=11180");
    check("Gravitation", "M=5.97E24 T=86400 | r=4.224E7");
    check("Springs and pendulums", "L=1 | T=2.007 f=0.4983");
    check("Springs and pendulums", "m=0.5 k=200 A=0.04 | w=20 T=0.3142 vmax=0.8 amax=16");
    check("Vector components", "Ax=3 Ay=4 | A=5 th=53.13");

    /* Chapter 8 exercises (numbers only) */
    check("Impulse", "m=0.0450 v1=0 v2=24.9 t=2.40E-3 | F=466.9");                          /* 8.7 */
    check("Impulse", "m=0.143 v1=-40.5 v2=54.5 t=1.50E-3 | J=13.585 F=9057");               /* 8.8 */
    check("Impulse", "m=0.160 v1=3.00 F=25.0 t=0.050 | v2=10.81");                          /* 8.9a */
    check("Impulse", "m=0.160 v1=3.00 F=-12.0 t=0.050 | v2=-0.75");                         /* 8.9b */
    check("Impulse in 2-D", "m=0.145 v1=40.0 v2=52.0 a2=150 t=1.75E-3 | Fx=-7046 Fy=2154"); /* 8.10 */
    check("Impulse", "m=80 v1=0 J=1200 | v2=15");                                           /* 8.12: J = area */
    check("Impulse", "m=2.00 v1=6.00 J=2.50 | v2=7.25");                                    /* 8.13 right */
    check("Impulse", "m=2.00 v1=6.00 J=-2.50 | v2=4.75");                                   /* 8.13 left */
    check("Momentum and KE", "m=40.0 v=20.0 | p=800 K=8000");                               /* 8.15ab */
    check("Impulse", "m=40.0 v1=20.0 v2=0 t=5.00 | F=-160");                                /* 8.15c */
    check("Recoil / any collision", "m1=65.5 v1=0 m2=2.50 v2=0 v2f=3.10 | v1f=-0.1183");    /* 8.16 */
    check("Recoil / any collision", "m1=63.78 v1=0 m2=73.98 v2=0 v2f=1.50 | v1f=-1.740 dK=179.8"); /* 8.18 */
    check("Recoil / any collision", "m1=5.85 v1=0 m2=1.65 v2=0 v1f=2.70 | v2f=-9.573 dK=96.9");    /* 8.19 */
    check("Collision: stick together", "m1=0.400 v1=10.0 m2=70.0 v2=0 | vf=0.05682");       /* 8.20a */
    check("Recoil / any collision", "m1=0.400 v1=10.0 m2=70.0 v2=0 v1f=-7.0 | v2f=0.09714"); /* 8.20b */
    check("Recoil / any collision", "m1=0.245 m2=0.360 v2=0 v1f=-0.118 v2f=0.660 | v1=0.8518 dK=-0.00877"); /* 8.21 */
    check("Recoil / any collision", "m1=1800 v1=1.60 m2=1450 v2=-1.10 v1f=0.270 | v2f=0.5510 dK=-2895.5"); /* 8.22 */
    check("Spring between 2 blocks", "m1=1.60 m2=1.60 k=185 x=0.200 | U=3.7 v1=1.521 v2=1.521");  /* 8.23 */
    check("Spring between 2 blocks", "m1=1.00 m2=3.00 v2=1.10 | v1=3.3 U=7.26");            /* 8.24 */
    check("Recoil / any collision", "m1=69.5 v1=0 m2=4.20E-3 v2=0 v2f=975 | v1f=-0.05892"); /* 8.25a */
    check("Collision in 2-D", "m1=45.0 v1=14.0 m2=70.0 v1f=8.00 b1=55.1 | v2f=7.381 b2=-34.86 dK=-1063"); /* 8.27 */
    check("Collision: stick together", "m1=60 v1=4.5 m2=2.5 v2=0 | vf=4.32");               /* 8.29a */
    check("Collision in 2-D", "m1=1 v1=40.0 m2=1 b1=30.0 b2=-45.0 | v1f=29.28 v2f=20.71 dK=-156.9"); /* 8.30 */
    check("Recoil / any collision", "m1=4.00 v1=2.00 m2=6.00 v2=-2.50 v2f=0.50 | v1f=-2.5 dK=-13.5"); /* 8.35 */
    check("Collision: stick together", "m1=1050 v1=-13.0 m2=6320 v2=11.0 | vf=7.581");      /* 8.36 */
    check("Stick together in 2-D", "m1=75 v1=7.2 m2=90 v2=8.8 a2=90 | vf=5.810 af=55.71");  /* 8.37 */
    check("Stick together in 2-D", "m1=1900 m2=1500 v2=16 a2=-90 af=-60 | v1=7.293 vf=8.151"); /* 8.38 */

    check("Collision: stick together", "m1=20.0 v1=7.5 m2=4.40 v2=0 | vf=6.148");          /* fish eats fish */
    check("Stick, then what happens", "m1=5.00E-3 m2=1.20 mu_k=0.20 d=0.310 | vf=1.1024 v1=265.7");   /* 8.42 */
    check("Stick, then what happens", "m1=12.0E-3 v1=380 m2=9.00 | h=0.013063 K1=866.4 Kf=1.1537");  /* 8.43 */
    check("Recoil / any collision", "m1=3.00 v1=8.00 m2=15.0 v2=0 v1f=-2.00 | v2f=2");      /* 8.44 step 1 */
    check("Energy conservation", "m=15.0 v1=2 h1=0 v2=0 h2=0 k=575 | x2=0.3230");           /* 8.44 step 2 */
    check("Stick, then what happens", "m1=3.00 v1=14.0 m2=4.00 L=1.40 | vf=6 T=248.6");     /* 8.45 */
    check("Collision: elastic", "m1=0.160 v1=0.90 m2=0.291 v2=-2.27 | v1f=-3.191 v2f=-0.02078"); /* 8.46 */
    check("Collision: stick together", "m1=6.00 v1=7.00 m2=10.50 v2=0 | vf=2.545 dK=-93.55"); /* 8.47a */
    check("Collision: elastic", "m1=6.00 v1=7.00 m2=10.50 v2=0 | v1f=-1.909 v2f=5.091");    /* 8.47b */
    check("Collision: elastic", "m1=0.0100 v1=-0.400 m2=0.0300 v2=0.200 | v1f=0.5 v2f=-0.1 dp1=0.009 dK1=4.5E-4"); /* 8.48 */
    check("Collision: elastic", "m1=1 v1=1 m2=2 v2=0 | v1f=-0.3333");                       /* 8.49a */
    check("Collision: elastic", "m1=1 v1=3.60E7 v2=0 v1f=-3.30E7 | m2=23 v2f=3.0E6");       /* 8.50 */
    check("Collision: elastic", "m2=1 v1=6.0 v2=0 v2f=3.0 | m1=0.3333 v1f=-3");             /* 8.51 */
    check("Center of mass", "m1=1.99E30 m2=1.90E27 x2=7.78E11 | xcm=7.421E8");              /* 8.52 */
    check("Center of mass", "m1=1.3312E10 m2=1.9531E9 x2=19700 | xcm=2520.5");              /* 8.53, m = d^3 */
    check("Center of mass", "m1=1200 m2=1800 x2=40.0 v1=12.0 v2=20.0 | xcm=24 vcm=16.8 P=50400"); /* 8.54 */
    check("Center of mass", "m1=0.500 x1=0.150 m2=0.800 x2=0.700 | xcm=0.4885");            /* 8.55 */
    check("Center of mass", "x1=0 m2=0.10 x2=8.0 xcm=2.0 v2=0 vcm=5.0 | m1=0.3 P=2 v1=6.667"); /* 8.56 */
    check("Rocket thrust", "m=180 a=0.029 vex=490 t=5.0 | F=5.22 dm=0.05327");              /* 8.61 */
    check("Rocket thrust", "R=5.40E-2 vex=1550 | F=83.7");                                  /* 8.62 */
    check("Impulse", "m=0.0360 h_drop=2.07 h_up=1.54 t=2.50E-3 | J=0.4271 F=170.8");        /* 8.64 */
    check("Impulse", "m=0.05714 v1=20.0 F=-380 t=3.00E-3 | J=-1.14 v2=0.049");              /* 8.65 x */
    check("Impulse", "m=0.05714 v1=-4.0 F=110 t=3.00E-3 | J=0.33 v2=1.775");                /* 8.65 y */
    check("Collision: stick together", "m1=2.00 v1=2.00 m2=10.00 v2=-0.500 | vf=-0.08333 dK=-5.208"); /* 8.67a */
    check("Collision: elastic", "m1=2.00 v1=2.00 m2=10.00 v2=-0.500 | v1f=-2.167 v2f=0.3333"); /* 8.67b */
    check("Recoil / any collision", "m1=130 v1=4.70 m2=20.0 v2=4.70 v2f=0 | v1f=5.423");    /* 8.68b */
    check("Collision: stick together", "m1=150 v1=4.70 m2=20.0 v2=-6.30 | vf=3.406");       /* 8.68c */
    check("Impulse", "m=3.00 v1=0 J=25 | v2=8.333");                                        /* 8.70: J = area */
    check("Stick, then what happens", "m1=0.500 v1=24.0 m2=8.00 Ht=2.20 | vf=1.4118 xr=0.9460"); /* 8.71 */
    check("Stick, then what happens", "m1=5.00 v1=12.0 m2=5.00 | vf=6 h=1.8367");           /* 8.73 */
    check("Force F = A t^n", "n=2 t2=1.25 F2=781.25 | A=500");                              /* 8.11a */
    check("Force F = A t^n", "A=500 n=2 t1=2.00 t2=3.50 m=2150 | J=5812.5 v2=2.7035");      /* 8.11bc */
    check("Force F = A t^n", "A=0.280 n=1 t2=2.00 p1=-3.00 | p2=-2.44");                    /* 8.14 x */
    check("Force F = A t^n", "A=-0.450 n=2 t2=2.00 p1=4.00 | p2=2.8");                      /* 8.14 y */
    check("Force F = A t^n", "A=3.00 n=1 m=2.00 v2=9.00 | t2=3.464 F2=10.39");              /* 8.66 */
    check("Rate of change (d/dt)", "B=0.9375 n=2 t=3.0 | dq/dt=5.625");                     /* 8.58c */
    check("Rate of change (d/dt)", "B=-0.75 n=2 t=2 | dq/dt=-3");                           /* 8.59 x */

    /* conflict detection */
    tp = &topics[find_topic("Straight-line motion")];
    memset(txt, 0, sizeof txt);
    strcpy(txt[0], "0"); strcpy(txt[1], "10"); strcpy(txt[2], "2"); strcpy(txt[3], "9");
    solve();
    if (conflict < 0) { fails++; printf("conflict not detected\n"); }

    if (argc > 1) {     /* show a steps page */
        check("Projectile", "v0=12 th=0 y=-20 |");
        build_steps(P_V); printf("%s\n", text);
        check("Energy conservation", "m=2 v1=3 h1=0 v2=0 h2=0 x2=0.1 |");
        build_steps(E_K); printf("%s\n", text);
    }
    printf(fails ? "\n%d FAILURES\n" : "\nALL TESTS PASSED\n", fails);
    return fails != 0;
}

/* ================================================================== */
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
enum { ST_HOME, ST_TOPICS, ST_FORM, ST_READ, ST_WAIT, ST_COMPOSE };
static uint8_t state = ST_HOME, read_back = ST_FORM;
static bool dirty = true;
static const char *notice = NULL;
static const char *read_title = "";
static bool read_is_ai = false;

static uint8_t home_sel, topic_sel, cur_cat, cur_topic, row, row_top;
static bool pending;            /* typed something; re-solve once the keys go quiet */
static uint32_t solve_at;

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
    print_right("PHYS1", SCREEN_W - TEXT_X, 5);
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

#if PART == 1
#define APP_TITLE "PHYS1 motion, forces, energy"
#elif PART == 2
#define APP_TITLE "PHYS2 momentum, collisions"
#elif PART == 3
#define APP_TITLE "PHYS3 rotation, gravity, fluids"
#else
#define APP_TITLE "PHYS all topics"
#endif

static void draw_home(void) {
#if USE_CLAUDE
    static const char *const soft[5] = { "Formula", "Ask AI", "", "Theme", "Quit" };
#else
    static const char *const soft[5] = { "Formula", "", "", "Theme", "Quit" };
#endif
    const char *items[NUM_CATS];
    for (uint8_t i = 0; i < NUM_CATS; i++) items[i] = cats[i].name;
    draw_list(APP_TITLE, "What is the problem about?", items, NUM_CATS, home_sel, soft);
}
static void draw_topics(void) {
    static const char *const soft[5] = { "Back", "", "", "", "" };
    const char *items[10];
    const cat_t *c = &cats[cur_cat];
    for (uint8_t i = 0; i < c->count; i++) items[i] = topics[cat_first(cur_cat) + i].name;
    draw_list(c->name, "Pick the kind of problem", items, c->count, topic_sel, soft);
}

static void draw_form(void) {
#if USE_CLAUDE
    static const char *const soft[5] = { "Menu", "Reset", "Tips", "Steps", "Ask AI" };
#else
    static const char *const soft[5] = { "Menu", "Reset", "Tips", "Steps", "" };
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
    } else {
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

static void open_topic(uint8_t t) {
    cur_topic = t;
    tp = &topics[t];
    memset(txt, 0, sizeof txt);
    row = row_top = 0;
    pending = false;
    solve();
    state = ST_FORM;
    dirty = true;
}
static void edited(void) { pending = true; solve_at = now_ms() + 250; dirty = true; }
static void settle(void) { if (pending) { solve(); pending = false; dirty = true; } }

static void show_steps(void) {
    settle();
    build_steps(row);
    open_reader("Steps", ST_FORM, false);
}
static void show_tips(void) {
    text_len = 0;
    tput("WORDS IN THE PROBLEM -> WHAT TO TYPE\n\n");
    tput(tp->tips);
    tput("\n\nFORMULAS USED HERE\n");
    put_formulas(tp);
    tput("\nKEYS\n(-) negative   ^ types E (3E8)\nCLEAR empties a value\nENTER on a solved value = steps\n");
    open_reader("Tips", ST_FORM, false);
}

static void form_keys(void) {
    int8_t d;
    char c = 0;
    uint8_t n = (uint8_t)strlen(txt[row]);

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
#endif
    else if (PRESSED(kb_KeyClear)) {
        if (n) { txt[row][0] = 0; edited(); }
        else state = ST_TOPICS;
        dirty = true;
    }
    else if (pressed_repeat(kb_KeyDel)) {
        if (n) { txt[row][n - 1] = 0; edited(); }
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
#endif
    if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyGraph)) return true;
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
