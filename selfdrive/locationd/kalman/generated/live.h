/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_3702839190652464607);
void inv_err_fun(double *nom_x, double *true_x, double *out_1567270458521939487);
void H_mod_fun(double *state, double *out_6209011822632597713);
void f_fun(double *state, double dt, double *out_688236372372532204);
void F_fun(double *state, double dt, double *out_1958849141401370296);
void h_3(double *state, double *unused, double *out_4037526428729242589);
void H_3(double *state, double *unused, double *out_2514133031888983281);
void h_4(double *state, double *unused, double *out_8199363292054432930);
void H_4(double *state, double *unused, double *out_1365864478765828093);
void h_9(double *state, double *unused, double *out_7747252342084683849);
void H_9(double *state, double *unused, double *out_7661839025497878621);
void h_10(double *state, double *unused, double *out_248830372333436934);
void H_10(double *state, double *unused, double *out_3756852166495392008);
void h_12(double *state, double *unused, double *out_8297322391528677603);
void H_12(double *state, double *unused, double *out_7558968725669421741);
void h_13(double *state, double *unused, double *out_4968379810315873373);
void H_13(double *state, double *unused, double *out_7617782100233727186);
void h_14(double *state, double *unused, double *out_7747252342084683849);
void H_14(double *state, double *unused, double *out_7661839025497878621);
void h_19(double *state, double *unused, double *out_7524281701416122087);
void H_19(double *state, double *unused, double *out_569876710731556009);
#define DIM 23
#define EDIM 22
#define MEDIM 22
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_3 = 3.841459;
void update_3(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_4 = 7.814728;
void update_4(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_9 = 7.814728;
void update_9(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_10 = 7.814728;
void update_10(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_12 = 7.814728;
void update_12(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_13 = 7.814728;
void update_13(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_14 = 7.814728;
void update_14(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_19 = 7.814728;
void update_19(double *, double *, double *, double *, double *);