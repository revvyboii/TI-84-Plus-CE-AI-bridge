/*
 * PHYS1 - Physics 1 problem solver for the TI-84 Plus CE
 *
 * Pick a topic, type the numbers a word problem gives you, and every other
 * quantity fills in by itself. ENTER on a solved value shows the steps.
 *
 * Two builds from this one file:
 *   USE_CLAUDE=1  PHYS1  local solver + "Ask AI" over USB (needs bridge.js)
 *   USE_CLAUDE=0  PHYSL  local solver only, no USB code at all
 *
 * Keys (form): digits . (-) type a value   ^ or , = E (3E8)
 *   UP/DOWN field   DEL backspace   CLEAR erase field / go back
 *   ENTER steps for this value      MODE light/dark
 *   Y= menu  WINDOW reset  ZOOM tips  TRACE steps  GRAPH Ask AI
 */

#ifndef USE_CLAUDE
#define USE_CLAUDE 1
#endif

#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <string.h>

/* ================================================================== */
/* Solver engine (plain C, also compiled on a PC for testing)          */
/* ================================================================== */

#define MAXV    13
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
enum { R_ABS = 1, R_NOCHK = 2 };     /* rule flags: sign-ambiguous / skip check */
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

#define FN(name, expr) static double name(const double *v) { return (expr); }
#define SQ(x) ((x) * (x))
static double root(double x) { return x < 0 ? UNK : sqrt(x); }

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
    s = root(b * b + 2 * a * dx);
    if (!ok(s)) return UNK;
    t1 = (-b + s) / a; t2 = (-b - s) / a;
    if (t1 > t2) { double x = t1; t1 = t2; t2 = x; }
    return t1 > 1e-9 ? t1 : t2;           /* first time it gets there */
}
static double k_vv2(const double *v) {
    double s = root(SQ(v[0]) + 2 * v[2] * v[4]);
    return v[4] < 0 ? -s : s;
}
static double k_v0v2(const double *v) {
    double s = root(SQ(v[1]) - 2 * v[2] * v[4]);
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
    {3, B(0)|B(2)|B(4), k_tq,   E2, R_NOCHK}, \
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
FN(p_ty,    (v[P_V0Y] + root(SQ(v[P_V0Y]) - 2 * GRAV * v[P_Y])) / GRAV)
FN(p_vyy,   -root(SQ(v[P_V0Y]) - 2 * GRAV * v[P_Y]))
FN(p_v0yvy, root(SQ(v[P_VY]) + 2 * GRAV * v[P_Y]))
FN(p_v,     sqrt(SQ(v[P_VX]) + SQ(v[P_VY])))
FN(p_h,     SQ(v[P_V0Y]) / (2 * GRAV))
FN(p_v0yh,  root(2 * GRAV * v[P_H]))
FN(p_r,     SQ(v[P_V0]) * sin(2 * v[P_TH] * DEG) / GRAV)
FN(p_v0r,   root(v[P_R] * GRAV / sin(2 * v[P_TH] * DEG)))

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
    {P_T,   B(P_V0Y)|B(P_Y), p_ty, "y = v0y t - (1/2) g t^2 (later t)", R_NOCHK},
    {P_VY,  B(P_V0Y)|B(P_Y), p_vyy, "vy^2 = v0y^2 - 2 g y (going down)", R_ABS},
    {P_V0Y, B(P_VY)|B(P_Y), p_v0yvy, "vy^2 = v0y^2 - 2 g y", R_ABS},
    {P_V,   B(P_VX)|B(P_VY), p_v, "v = sqrt(vx^2 + vy^2)", 0},
    {P_H,   B(P_V0Y), p_h, "H = v0y^2 / (2 g)", 0},
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
};

/* ---------- 10. conservation of energy ----------------------------- */
enum { E_M, E_V1, E_H1, E_V2, E_H2, E_X1, E_X2, E_K, E_W };
#define E_EQ "K1 + Ug1 + Us1 + W = K2 + Ug2 + Us2"
/* extra energy per kg from spring + other work; UNK if m is needed but absent */
#define ESPR (v[E_K] * (SQ(v[E_X1]) - SQ(v[E_X2])) + 2 * v[E_W])
#define EALL (v[E_M] * (SQ(v[E_V1]) - SQ(v[E_V2])) + 2 * v[E_M] * GRAV * (v[E_H1] - v[E_H2]) + 2 * v[E_W])
FN(e_v2,  root(SQ(v[E_V1]) + 2 * GRAV * (v[E_H1] - v[E_H2]) + ESPR / v[E_M]))
FN(e_v1,  root(SQ(v[E_V2]) - 2 * GRAV * (v[E_H1] - v[E_H2]) - ESPR / v[E_M]))
FN(e_h2,  v[E_H1] + (SQ(v[E_V1]) - SQ(v[E_V2]) + ESPR / v[E_M]) / (2 * GRAV))
FN(e_h1,  v[E_H2] - (SQ(v[E_V1]) - SQ(v[E_V2]) + ESPR / v[E_M]) / (2 * GRAV))
FN(e_v2n, ESPR != 0 ? UNK : root(SQ(v[E_V1]) + 2 * GRAV * (v[E_H1] - v[E_H2])))
FN(e_v1n, ESPR != 0 ? UNK : root(SQ(v[E_V2]) - 2 * GRAV * (v[E_H1] - v[E_H2])))
FN(e_h2n, ESPR != 0 ? UNK : v[E_H1] + (SQ(v[E_V1]) - SQ(v[E_V2])) / (2 * GRAV))
FN(e_h1n, ESPR != 0 ? UNK : v[E_H2] - (SQ(v[E_V1]) - SQ(v[E_V2])) / (2 * GRAV))
FN(e_w,   0.5 * v[E_M] * (SQ(v[E_V2]) - SQ(v[E_V1])) + v[E_M] * GRAV * (v[E_H2] - v[E_H1])
          + 0.5 * v[E_K] * (SQ(v[E_X2]) - SQ(v[E_X1])))
FN(e_m,   ESPR == 0 ? UNK : -ESPR / (SQ(v[E_V1]) - SQ(v[E_V2]) + 2 * GRAV * (v[E_H1] - v[E_H2])))
FN(e_k,   EALL / (SQ(v[E_X2]) - SQ(v[E_X1])))
FN(e_x2,  root(SQ(v[E_X1]) + EALL / v[E_K]))
FN(e_x1,  root(SQ(v[E_X2]) - EALL / v[E_K]))
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
FN(g_m1,   v[2] * (v[5] - v[3]) / (v[1] - v[4]))
FN(g_m2,   v[0] * (v[1] - v[4]) / (v[5] - v[3]))
static const var_t elas_vars[] = {
    {"m1",  "kg",  "m1: first object", UNK},
    {"v1",  "m/s", "v1: its velocity before (left = -)", UNK},
    {"m2",  "kg",  "m2: second object", UNK},
    {"v2",  "m/s", "v2: its velocity before (rest: 0)", UNK},
    {"v1f", "m/s", "v1f: velocity of m1 after", UNK},
    {"v2f", "m/s", "v2f: velocity of m2 after", UNK},
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
};
static const rule_t gen_rules[] = {
    {5, B(0)|B(1)|B(2)|B(3)|B(4), g_v2f, G_EQ, 0},
    {4, B(0)|B(1)|B(2)|B(3)|B(5), g_v1f, G_EQ, 0},
    {1, B(0)|B(2)|B(3)|B(4)|B(5), g_v1, G_EQ, 0},
    {3, B(0)|B(1)|B(2)|B(4)|B(5), g_v2, G_EQ, 0},
    {0, B(1)|B(2)|B(3)|B(4)|B(5), g_m1, G_EQ, 0},
    {2, B(0)|B(1)|B(3)|B(4)|B(5), g_m2, G_EQ, 0},
};

/* ---------- 14. impulse -------------------------------------------- */
enum { J_M, J_V1, J_V2, J_J, J_F, J_T };
FN(j_j,   v[J_M] * (v[J_V2] - v[J_V1]))
FN(j_v2,  v[J_V1] + v[J_J] / v[J_M])
FN(j_v1,  v[J_V2] - v[J_J] / v[J_M])
FN(j_m,   v[J_J] / (v[J_V2] - v[J_V1]))
FN(j_jf,  v[J_F] * v[J_T])
FN(j_f,   v[J_J] / v[J_T])
FN(j_t,   v[J_J] / v[J_F])
static const var_t imp_vars[] = {
    {"m",  "kg",   "m: mass", UNK},
    {"v1", "m/s",  "v1: velocity before (left = -)", UNK},
    {"v2", "m/s",  "v2: velocity after (bounce back = -)", UNK},
    {"J",  "kg m/s", "J: impulse = change in momentum", UNK},
    {"F",  "N",    "F: average force", UNK},
    {"t",  "s",    "t: contact time (ms: use E-3)", UNK},
};
static const rule_t imp_rules[] = {
    {J_J, B(J_M)|B(J_V1)|B(J_V2), j_j, "J = m (v2 - v1)", 0},
    {J_J, B(J_F)|B(J_T), j_jf, "J = F t", 0},
    {J_F, B(J_J)|B(J_T), j_f, "J = F t", 0},
    {J_T, B(J_J)|B(J_F), j_t, "J = F t", 0},
    {J_V2, B(J_M)|B(J_V1)|B(J_J), j_v2, "J = m (v2 - v1)", 0},
    {J_V1, B(J_M)|B(J_V2)|B(J_J), j_v1, "J = m (v2 - v1)", 0},
    {J_M, B(J_V1)|B(J_V2)|B(J_J), j_m, "J = m (v2 - v1)", 0},
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
};

/* ---------- 17. moment of inertia of common shapes ----------------- */
FN(mi_hoop,  v[0] * SQ(v[1]))
FN(mi_disk,  0.5 * v[0] * SQ(v[1]))
FN(mi_sph,   0.4 * v[0] * SQ(v[1]))
FN(mi_shell, 2.0 / 3.0 * v[0] * SQ(v[1]))
FN(mi_rodc,  v[0] * SQ(v[2]) / 12.0)
FN(mi_rode,  v[0] * SQ(v[2]) / 3.0)
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
};

/* ---------- 18. gravitation ---------------------------------------- */
enum { G_M1, G_M2, G_R, G_F, G_G, G_VO, G_T, G_VE };
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
static const var_t grav_vars[] = {
    {"M",    "kg",    "M: planet (Earth: 5.97E24)", UNK},
    {"m",    "kg",    "m: small object or satellite", UNK},
    {"r",    "m",     "r: center to center (Earth R 6.37E6)", UNK},
    {"F",    "N",     "F: gravitational force", UNK},
    {"g",    "m/s^2", "g: gravity strength at r", UNK},
    {"v_orb","m/s",   "v_orb: speed of a circular orbit", UNK},
    {"T",    "s",     "T: time for one orbit", UNK},
    {"v_esc","m/s",   "v_esc: escape speed from r", UNK},
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
};

/* ---------- 19. springs and pendulums ------------------------------ */
enum { H_M, H_K, H_L, H_T, H_F, H_W, H_A, H_VM, H_AM };
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
};

/* ---------- topic table -------------------------------------------- */
#define N(a) ((uint8_t)(sizeof(a) / sizeof((a)[0])))
#define TOPIC(name, vars, rules, tips) { name, vars, N(vars), rules, N(rules), tips }
static const topic_t topics[] = {
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
    TOPIC("Vector components", vec_vars, vec_rules,
        "30 deg above +x: th = 30\n30 deg north of west: th = 150\n"
        "straight down: th = 270 (or -90)\nEnter Ax and Ay to get size + angle."),
    /* Forces */
    TOPIC("F = ma + friction", force_vars, force_rules,
        "frictionless / smooth: mu_k = 0\nconstant speed: a = 0\n"
        "pulled with 50 N: F = 50\nweighs 98 N: w = 98\n"
        "For a = from motion, use\nStraight-line motion first."),
    TOPIC("Incline (ramp)", incl_vars, incl_rules,
        "frictionless ramp: mu_k = 0\nslides down: read a\n"
        "pushed up and slows: read a_up\njust about to slip: read th_s\n"
        "mass is only needed for forces."),
    TOPIC("Atwood (2 hanging)", atw_vars, atw_rules,
        "two masses over a pulley,\nboth hanging. m2 is the heavy one.\n"
        "massless string + pulley assumed."),
    TOPIC("Table + hanging mass", tab_vars, tab_rules,
        "m1 slides on the table,\nm2 hangs off the edge.\n"
        "frictionless table: mu_k = 0"),
    /* Energy */
    TOPIC("Work and power", work_vars, work_rules,
        "force along the motion: phi = 0\n(blank phi counts as 0)\n"
        "friction: phi = 180 (negative W)\nlifting at steady speed: F = m g\n"
        "2 hp motor: P = 1492"),
    TOPIC("Energy conservation", en_vars, en_rules,
        "released from rest: v1 = 0\nreaches the ground: h2 = 0\n"
        "highest point: v2 = 0\nno friction: leave W blank\n"
        "friction does 30 J: W = -30\nno spring: leave x1, x2, k blank\n"
        "mass is not needed unless there\nis a spring or W."),
    /* Momentum */
    TOPIC("Collision: stick together", stick_vars, stick_rules,
        "at rest / stationary: v = 0\nmoving left / toward: negative v\n"
        "couple, stick, embed, catch,\njump onto: use this one\n"
        "g to kg: 5 g = 5E-3"),
    TOPIC("Collision: elastic", elas_vars, elas_rules,
        "elastic / bounce with no energy\nlost: use this one\n"
        "at rest: v = 0\nmoving left: negative v\n"
        "head-on = 1-D, signs matter"),
    TOPIC("Recoil / any collision", gen_vars, gen_rules,
        "Fill in 5 of the 6, get the 6th.\nrecoil, explosion, push off:\n"
        "  both start at rest: v1 = v2 = 0\nbounces back: negative after\n"
        "passes through: enter its v after"),
    TOPIC("Impulse", imp_vars, imp_rules,
        "bounces back: v2 is negative\nstops: v2 = 0\n"
        "in contact 2 ms: t = 2E-3\naverage force: read F"),
    /* Rotation */
    TOPIC("Rotational motion", rot_vars, rot_rules,
        "starts from rest: w0 = 0\nstops: w = 0\n"
        "5 revolutions: dth = 31.42\nrpm to rad/s: times 0.1047"),
    TOPIC("Torque and rotation", torq_vars, torq_rules,
        "force at a right angle: leave\nphi blank (counts as 90)\n"
        "need I? use Moment of inertia\nwrench 0.3 m long: r = 0.3"),
    TOPIC("Moment of inertia", moi_vars, moi_rules,
        "Enter M and R (or L for a rod).\nEvery shape's I shows up;\npick the one that matches."),
    /* Gravity + oscillations */
    TOPIC("Gravitation", grav_vars, grav_rules,
        "Earth: M = 5.97E24, R = 6.37E6\nMoon: M = 7.35E22, R = 1.74E6\n"
        "Sun: M = 1.99E30\n300 km up: r = 6.37E6 + 3E5\n   = 6.67E6\n^ key types the E"),
    TOPIC("Springs and pendulums", shm_vars, shm_rules,
        "spring: enter m and k\npendulum: enter L only\n"
        "10 cycles in 5 s: T = 0.5\npulled back 4 cm: A = 0.04"),
};
#define NUM_TOPICS N(topics)

typedef struct { const char *name; uint8_t first, count; } cat_t;
static const cat_t cats[] = {
    {"Motion", 0, 4}, {"Forces", 4, 4}, {"Energy", 8, 2},
    {"Momentum", 10, 4}, {"Rotation", 14, 3}, {"Gravity + springs", 17, 2},
};
#define NUM_CATS N(cats)

/* ---------- form state + solver ------------------------------------ */
enum { S_UNK, S_USER, S_SOLVED, S_DEF };
#define VAL_LEN 11

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

static void propagate(void) {
    bool changed = true;
    while (changed) {
        uint16_t k = known_mask();
        changed = false;
        for (uint8_t r = 0; r < tp->nrules; r++) {
            const rule_t *R = &tp->rules[r];
            double x;
            if (st[R->target] != S_UNK || (R->need & k) != R->need) continue;
            x = R->fn(val);
            if (!ok(x)) continue;
            if (fabs(x) < 1e-9) x = 0;
            val[R->target] = x;
            st[R->target] = S_SOLVED;
            by[R->target] = r;
            k |= B(R->target);
            changed = true;
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

static void solve(void) {
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
        if (R->flags & R_ABS) { x = x * x; y = y * y; }     /* square-root rules: compare squares */
        tol = 2e-3 * (fabs(x) > fabs(y) ? fabs(x) : fabs(y)) + 1e-6;
        for (uint8_t j = 0; j < tp->nvars; j++) {    /* allow for rounding in what was typed */
            double sc = fabs(val[j]);
            if (R->flags & R_ABS) sc *= sc;
            if ((R->need & B(j)) && 1e-4 * sc > tol) tol = 1e-4 * sc;
        }
        if (fabs(x - y) > tol) { conflict = r; break; }
    }
}

/* ---------- text builder (steps, tips, prompts) -------------------- */
#define MAX_TEXT 4096
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
        tput("\nNot enough information yet for ");
        tput(tp->vars[target].sym);
        tput(".\nEnter more of what the problem gives you, or check Tips (ZOOM).\n");
    } else {
        tput("\n"); tput(tp->vars[target].sym);
        tput(" is a number you typed in.\nMove to a solved (colored) value and press ENTER to see its steps.\n");
    }
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

int main(int argc, char **argv) {
    char b[32];
    (void)argv; (void)cats;
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
    full("Impulse", "m=0.5 v1=8 v2=-6 t=0.004");
    full("Torque and rotation", "F=20 r=0.4 phi=60 I=0.5 w=6");
    full("Moment of inertia", "M=3 R=0.2 L=1.5");
    full("Gravitation", "M=5.97E24 m=1000 r=6.67E6");
    full("Springs and pendulums", "m=0.5 k=200 A=0.04");

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
#define ROWS       11
#define HELP_Y     194
#define SOFT_Y     222
#define MAX_LINES  480

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
        int y = BODY_Y + 20 + i * 22;
        char num[4] = { (char)('1' + i), 0, 0, 0 };
        if (i == sel) {
            gfx_SetColor(C_ACCENT);
            gfx_FillRectangle_NoClip(TEXT_X - 2, (uint8_t)(y - 6), SCREEN_W - 2 * TEXT_X + 4, 20);
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

#if USE_CLAUDE
#define HOME_ITEMS (NUM_CATS + 1)
#else
#define HOME_ITEMS NUM_CATS
#endif

static void draw_home(void) {
    static const char *const soft[5] = { "", "", "", "Theme", "Quit" };
    const char *items[NUM_CATS + 1];
    for (uint8_t i = 0; i < NUM_CATS; i++) items[i] = cats[i].name;
    items[NUM_CATS] = "Ask Claude (type a question)";
    draw_list("PHYS1  Physics 1 solver", "What is the problem about?", items, HOME_ITEMS, home_sel, soft);
}
static void draw_topics(void) {
    static const char *const soft[5] = { "Back", "", "", "", "" };
    const char *items[6];
    const cat_t *c = &cats[cur_cat];
    for (uint8_t i = 0; i < c->count; i++) items[i] = topics[c->first + i].name;
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
            gfx_SetTextFGColor(st[i] == S_USER || !txt[i][0] ? C_TEXT : C_WARN);
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
    if (row_top > 0) gfx_PrintStringXY("^", 236, BODY_Y - 5);
    if (row_top + ROWS < tp->nvars) gfx_PrintStringXY("v more", 236, HELP_Y - 8);

    /* help line for the selected value, or a warning */
    gfx_SetColor(C_RULE);
    gfx_HorizLine_NoClip(TEXT_X, HELP_Y - 2, SCREEN_W - 2 * TEXT_X);
    if (conflict >= 0) {
        gfx_SetTextFGColor(C_WARN);
        gfx_PrintStringXY("Numbers disagree. Check:", TEXT_X, HELP_Y + 2);
        print_n(tp->rules[conflict].eq, (uint8_t)strlen(tp->rules[conflict].eq), TEXT_X, HELP_Y + 13);
    } else {
        gfx_SetTextFGColor(C_TEXT);
        gfx_PrintStringXY(tp->vars[row].desc, TEXT_X, HELP_Y + 2);
        gfx_SetTextFGColor(C_DIM);
        gfx_PrintStringXY(st[row] == S_SOLVED ? "ENTER: show the steps for this"
                          : "type the givens; blanks get solved", TEXT_X, HELP_Y + 13);
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
    solve();
    state = ST_FORM;
    dirty = true;
}
static void show_steps(void) {
    build_steps(row);
    open_reader("Steps", ST_FORM, false);
}
static void show_tips(void) {
    text_len = 0;
    tput("WORDS IN THE PROBLEM -> WHAT TO TYPE\n\n");
    tput(tp->tips);
    tput("\n\nKEYS\n(-) negative   ^ types E (3E8)\nCLEAR empties a value\nENTER on a solved value = steps\n");
    open_reader("Tips", ST_FORM, false);
}

static void form_keys(void) {
    int8_t d;
    char c = 0;
    uint8_t n = (uint8_t)strlen(txt[row]);

    if (pressed_repeat(kb_KeyDown)) { row = (uint8_t)((row + 1) % tp->nvars); dirty = true; }
    else if (pressed_repeat(kb_KeyUp)) { row = (uint8_t)((row + tp->nvars - 1) % tp->nvars); dirty = true; }
    else if (PRESSED(kb_KeyEnter)) {
        if (st[row] == S_SOLVED) show_steps();
        else { row = (uint8_t)((row + 1) % tp->nvars); dirty = true; }
    }
    else if (PRESSED(kb_KeyTrace)) show_steps();
    else if (PRESSED(kb_KeyZoom)) show_tips();
    else if (PRESSED(kb_KeyWindow)) { memset(txt, 0, sizeof txt); solve(); dirty = true; }
    else if (PRESSED(kb_KeyYequ)) { state = ST_TOPICS; dirty = true; }
#if USE_CLAUDE
    else if (PRESSED(kb_KeyGraph)) ask_ai_form();
#endif
    else if (PRESSED(kb_KeyClear)) {
        if (n) { txt[row][0] = 0; solve(); }
        else state = ST_TOPICS;
        dirty = true;
    }
    else if (pressed_repeat(kb_KeyDel)) {
        if (n) { txt[row][n - 1] = 0; solve(); dirty = true; }
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
        solve();
        dirty = true;
    }
    if (row < row_top) row_top = row;
    if (row >= row_top + ROWS) row_top = (uint8_t)(row - ROWS + 1);
}

static bool home_keys(void) {       /* true = quit */
    int8_t d = digit_pressed(), pick = -1;
    if (pressed_repeat(kb_KeyDown)) { home_sel = (uint8_t)((home_sel + 1) % HOME_ITEMS); dirty = true; }
    if (pressed_repeat(kb_KeyUp))   { home_sel = (uint8_t)((home_sel + HOME_ITEMS - 1) % HOME_ITEMS); dirty = true; }
    if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyRight)) pick = (int8_t)home_sel;
    if (d >= 1 && d <= (int8_t)HOME_ITEMS) pick = (int8_t)(d - 1);
    if (PRESSED(kb_KeyTrace) || PRESSED(kb_KeyMode)) { night = !night; apply_theme(); dirty = true; }
    if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyGraph)) return true;
    if (pick < 0) return false;
    home_sel = (uint8_t)pick;
    if (pick < (int8_t)NUM_CATS) { cur_cat = (uint8_t)pick; topic_sel = 0; state = ST_TOPICS; }
#if USE_CLAUDE
    else { read_back = ST_HOME; state = ST_COMPOSE; }
#endif
    dirty = true;
    return false;
}
static void topics_keys(void) {
    const cat_t *c = &cats[cur_cat];
    int8_t d = digit_pressed();
    if (pressed_repeat(kb_KeyDown)) { topic_sel = (uint8_t)((topic_sel + 1) % c->count); dirty = true; }
    else if (pressed_repeat(kb_KeyUp)) { topic_sel = (uint8_t)((topic_sel + c->count - 1) % c->count); dirty = true; }
    else if (PRESSED(kb_KeyEnter) || PRESSED(kb_KeyRight)) open_topic((uint8_t)(c->first + topic_sel));
    else if (PRESSED(kb_KeyClear) || PRESSED(kb_KeyYequ) || PRESSED(kb_KeyLeft)) { state = ST_HOME; dirty = true; }
    else if (d >= 1 && d <= (int8_t)c->count) { topic_sel = (uint8_t)(d - 1); open_topic((uint8_t)(c->first + topic_sel)); }
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
        case ST_FORM:   form_keys(); break;
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
