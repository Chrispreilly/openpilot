/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_862951853940638262);
void inv_err_fun(double *nom_x, double *true_x, double *out_4196267415493599366);
void H_mod_fun(double *state, double *out_4972867834322587385);
void f_fun(double *state, double dt, double *out_6206102444964281666);
void F_fun(double *state, double dt, double *out_2039902045151723830);
void h_25(double *state, double *unused, double *out_3210887464065636601);
void H_25(double *state, double *unused, double *out_7660322170446165983);
void h_24(double *state, double *unused, double *out_838322440835917861);
void H_24(double *state, double *unused, double *out_589892170851736179);
void h_30(double *state, double *unused, double *out_6752297139885956667);
void H_30(double *state, double *unused, double *out_1455446809497690677);
void h_26(double *state, double *unused, double *out_4875757984572891982);
void H_26(double *state, double *unused, double *out_3055641902892824233);
void h_27(double *state, double *unused, double *out_4930031343789200054);
void H_27(double *state, double *unused, double *out_3944297225607655941);
void h_29(double *state, double *unused, double *out_5961981834362151480);
void H_29(double *state, double *unused, double *out_2969857451529201489);
void h_28(double *state, double *unused, double *out_8205368988943805089);
void H_28(double *state, double *unused, double *out_3045396279621735753);
#define DIM 8
#define EDIM 8
#define MEDIM 8
typedef void (*Hfun)(double *, double *, double *);

void predict(double *x, double *P, double *Q, double dt);
const static double MAHA_THRESH_25 = 3.841459;
void update_25(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_24 = 5.991465;
void update_24(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_30 = 3.841459;
void update_30(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_26 = 3.841459;
void update_26(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_27 = 3.841459;
void update_27(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_29 = 3.841459;
void update_29(double *, double *, double *, double *, double *);
const static double MAHA_THRESH_28 = 5.991465;
void update_28(double *, double *, double *, double *, double *);
void set_mass(double x);

void set_rotational_inertia(double x);

void set_center_to_front(double x);

void set_center_to_rear(double x);

void set_stiffness_front(double x);

void set_stiffness_rear(double x);
