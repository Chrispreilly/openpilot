
extern "C"{

double mass;

void set_mass(double x){ mass = x;}

double rotational_inertia;

void set_rotational_inertia(double x){ rotational_inertia = x;}

double center_to_front;

void set_center_to_front(double x){ center_to_front = x;}

double center_to_rear;

void set_center_to_rear(double x){ center_to_rear = x;}

double stiffness_front;

void set_stiffness_front(double x){ stiffness_front = x;}

double stiffness_rear;

void set_stiffness_rear(double x){ stiffness_rear = x;}

}
extern "C" {
#include <math.h>
/******************************************************************************
 *                       Code generated with sympy 1.4                        *
 *                                                                            *
 *              See http://www.sympy.org/ for more information.               *
 *                                                                            *
 *                         This file is part of 'ekf'                         *
 ******************************************************************************/
void err_fun(double *nom_x, double *delta_x, double *out_862951853940638262) {
   out_862951853940638262[0] = delta_x[0] + nom_x[0];
   out_862951853940638262[1] = delta_x[1] + nom_x[1];
   out_862951853940638262[2] = delta_x[2] + nom_x[2];
   out_862951853940638262[3] = delta_x[3] + nom_x[3];
   out_862951853940638262[4] = delta_x[4] + nom_x[4];
   out_862951853940638262[5] = delta_x[5] + nom_x[5];
   out_862951853940638262[6] = delta_x[6] + nom_x[6];
   out_862951853940638262[7] = delta_x[7] + nom_x[7];
}
void inv_err_fun(double *nom_x, double *true_x, double *out_4196267415493599366) {
   out_4196267415493599366[0] = -nom_x[0] + true_x[0];
   out_4196267415493599366[1] = -nom_x[1] + true_x[1];
   out_4196267415493599366[2] = -nom_x[2] + true_x[2];
   out_4196267415493599366[3] = -nom_x[3] + true_x[3];
   out_4196267415493599366[4] = -nom_x[4] + true_x[4];
   out_4196267415493599366[5] = -nom_x[5] + true_x[5];
   out_4196267415493599366[6] = -nom_x[6] + true_x[6];
   out_4196267415493599366[7] = -nom_x[7] + true_x[7];
}
void H_mod_fun(double *state, double *out_4972867834322587385) {
   out_4972867834322587385[0] = 1.0;
   out_4972867834322587385[1] = 0.0;
   out_4972867834322587385[2] = 0.0;
   out_4972867834322587385[3] = 0.0;
   out_4972867834322587385[4] = 0.0;
   out_4972867834322587385[5] = 0.0;
   out_4972867834322587385[6] = 0.0;
   out_4972867834322587385[7] = 0.0;
   out_4972867834322587385[8] = 0.0;
   out_4972867834322587385[9] = 1.0;
   out_4972867834322587385[10] = 0.0;
   out_4972867834322587385[11] = 0.0;
   out_4972867834322587385[12] = 0.0;
   out_4972867834322587385[13] = 0.0;
   out_4972867834322587385[14] = 0.0;
   out_4972867834322587385[15] = 0.0;
   out_4972867834322587385[16] = 0.0;
   out_4972867834322587385[17] = 0.0;
   out_4972867834322587385[18] = 1.0;
   out_4972867834322587385[19] = 0.0;
   out_4972867834322587385[20] = 0.0;
   out_4972867834322587385[21] = 0.0;
   out_4972867834322587385[22] = 0.0;
   out_4972867834322587385[23] = 0.0;
   out_4972867834322587385[24] = 0.0;
   out_4972867834322587385[25] = 0.0;
   out_4972867834322587385[26] = 0.0;
   out_4972867834322587385[27] = 1.0;
   out_4972867834322587385[28] = 0.0;
   out_4972867834322587385[29] = 0.0;
   out_4972867834322587385[30] = 0.0;
   out_4972867834322587385[31] = 0.0;
   out_4972867834322587385[32] = 0.0;
   out_4972867834322587385[33] = 0.0;
   out_4972867834322587385[34] = 0.0;
   out_4972867834322587385[35] = 0.0;
   out_4972867834322587385[36] = 1.0;
   out_4972867834322587385[37] = 0.0;
   out_4972867834322587385[38] = 0.0;
   out_4972867834322587385[39] = 0.0;
   out_4972867834322587385[40] = 0.0;
   out_4972867834322587385[41] = 0.0;
   out_4972867834322587385[42] = 0.0;
   out_4972867834322587385[43] = 0.0;
   out_4972867834322587385[44] = 0.0;
   out_4972867834322587385[45] = 1.0;
   out_4972867834322587385[46] = 0.0;
   out_4972867834322587385[47] = 0.0;
   out_4972867834322587385[48] = 0.0;
   out_4972867834322587385[49] = 0.0;
   out_4972867834322587385[50] = 0.0;
   out_4972867834322587385[51] = 0.0;
   out_4972867834322587385[52] = 0.0;
   out_4972867834322587385[53] = 0.0;
   out_4972867834322587385[54] = 1.0;
   out_4972867834322587385[55] = 0.0;
   out_4972867834322587385[56] = 0.0;
   out_4972867834322587385[57] = 0.0;
   out_4972867834322587385[58] = 0.0;
   out_4972867834322587385[59] = 0.0;
   out_4972867834322587385[60] = 0.0;
   out_4972867834322587385[61] = 0.0;
   out_4972867834322587385[62] = 0.0;
   out_4972867834322587385[63] = 1.0;
}
void f_fun(double *state, double dt, double *out_6206102444964281666) {
   out_6206102444964281666[0] = state[0];
   out_6206102444964281666[1] = state[1];
   out_6206102444964281666[2] = state[2];
   out_6206102444964281666[3] = state[3];
   out_6206102444964281666[4] = state[4];
   out_6206102444964281666[5] = dt*((-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]))*state[6] + stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*state[1]) + (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*state[4])) + state[5];
   out_6206102444964281666[6] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*state[4])) + state[6];
   out_6206102444964281666[7] = state[7];
}
void F_fun(double *state, double dt, double *out_2039902045151723830) {
   out_2039902045151723830[0] = 1;
   out_2039902045151723830[1] = 0;
   out_2039902045151723830[2] = 0;
   out_2039902045151723830[3] = 0;
   out_2039902045151723830[4] = 0;
   out_2039902045151723830[5] = 0;
   out_2039902045151723830[6] = 0;
   out_2039902045151723830[7] = 0;
   out_2039902045151723830[8] = 0;
   out_2039902045151723830[9] = 1;
   out_2039902045151723830[10] = 0;
   out_2039902045151723830[11] = 0;
   out_2039902045151723830[12] = 0;
   out_2039902045151723830[13] = 0;
   out_2039902045151723830[14] = 0;
   out_2039902045151723830[15] = 0;
   out_2039902045151723830[16] = 0;
   out_2039902045151723830[17] = 0;
   out_2039902045151723830[18] = 1;
   out_2039902045151723830[19] = 0;
   out_2039902045151723830[20] = 0;
   out_2039902045151723830[21] = 0;
   out_2039902045151723830[22] = 0;
   out_2039902045151723830[23] = 0;
   out_2039902045151723830[24] = 0;
   out_2039902045151723830[25] = 0;
   out_2039902045151723830[26] = 0;
   out_2039902045151723830[27] = 1;
   out_2039902045151723830[28] = 0;
   out_2039902045151723830[29] = 0;
   out_2039902045151723830[30] = 0;
   out_2039902045151723830[31] = 0;
   out_2039902045151723830[32] = 0;
   out_2039902045151723830[33] = 0;
   out_2039902045151723830[34] = 0;
   out_2039902045151723830[35] = 0;
   out_2039902045151723830[36] = 1;
   out_2039902045151723830[37] = 0;
   out_2039902045151723830[38] = 0;
   out_2039902045151723830[39] = 0;
   out_2039902045151723830[40] = dt*(stiffness_front*(-state[2] - state[3] + state[7])/(mass*state[1]) + (-stiffness_front - stiffness_rear)*state[5]/(mass*state[4]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[6]/(mass*state[4]));
   out_2039902045151723830[41] = -dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(mass*pow(state[1], 2));
   out_2039902045151723830[42] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2039902045151723830[43] = -dt*stiffness_front*state[0]/(mass*state[1]);
   out_2039902045151723830[44] = dt*((-1 - (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*pow(state[4], 2)))*state[6] - (-stiffness_front*state[0] - stiffness_rear*state[0])*state[5]/(mass*pow(state[4], 2)));
   out_2039902045151723830[45] = dt*(-stiffness_front*state[0] - stiffness_rear*state[0])/(mass*state[4]) + 1;
   out_2039902045151723830[46] = dt*(-state[4] + (-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(mass*state[4]));
   out_2039902045151723830[47] = dt*stiffness_front*state[0]/(mass*state[1]);
   out_2039902045151723830[48] = dt*(center_to_front*stiffness_front*(-state[2] - state[3] + state[7])/(rotational_inertia*state[1]) + (-center_to_front*stiffness_front + center_to_rear*stiffness_rear)*state[5]/(rotational_inertia*state[4]) + (-pow(center_to_front, 2)*stiffness_front - pow(center_to_rear, 2)*stiffness_rear)*state[6]/(rotational_inertia*state[4]));
   out_2039902045151723830[49] = -center_to_front*dt*stiffness_front*(-state[2] - state[3] + state[7])*state[0]/(rotational_inertia*pow(state[1], 2));
   out_2039902045151723830[50] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2039902045151723830[51] = -center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2039902045151723830[52] = dt*(-(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])*state[5]/(rotational_inertia*pow(state[4], 2)) - (-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])*state[6]/(rotational_inertia*pow(state[4], 2)));
   out_2039902045151723830[53] = dt*(-center_to_front*stiffness_front*state[0] + center_to_rear*stiffness_rear*state[0])/(rotational_inertia*state[4]);
   out_2039902045151723830[54] = dt*(-pow(center_to_front, 2)*stiffness_front*state[0] - pow(center_to_rear, 2)*stiffness_rear*state[0])/(rotational_inertia*state[4]) + 1;
   out_2039902045151723830[55] = center_to_front*dt*stiffness_front*state[0]/(rotational_inertia*state[1]);
   out_2039902045151723830[56] = 0;
   out_2039902045151723830[57] = 0;
   out_2039902045151723830[58] = 0;
   out_2039902045151723830[59] = 0;
   out_2039902045151723830[60] = 0;
   out_2039902045151723830[61] = 0;
   out_2039902045151723830[62] = 0;
   out_2039902045151723830[63] = 1;
}
void h_25(double *state, double *unused, double *out_3210887464065636601) {
   out_3210887464065636601[0] = state[6];
}
void H_25(double *state, double *unused, double *out_7660322170446165983) {
   out_7660322170446165983[0] = 0;
   out_7660322170446165983[1] = 0;
   out_7660322170446165983[2] = 0;
   out_7660322170446165983[3] = 0;
   out_7660322170446165983[4] = 0;
   out_7660322170446165983[5] = 0;
   out_7660322170446165983[6] = 1;
   out_7660322170446165983[7] = 0;
}
void h_24(double *state, double *unused, double *out_838322440835917861) {
   out_838322440835917861[0] = state[4];
   out_838322440835917861[1] = state[5];
}
void H_24(double *state, double *unused, double *out_589892170851736179) {
   out_589892170851736179[0] = 0;
   out_589892170851736179[1] = 0;
   out_589892170851736179[2] = 0;
   out_589892170851736179[3] = 0;
   out_589892170851736179[4] = 1;
   out_589892170851736179[5] = 0;
   out_589892170851736179[6] = 0;
   out_589892170851736179[7] = 0;
   out_589892170851736179[8] = 0;
   out_589892170851736179[9] = 0;
   out_589892170851736179[10] = 0;
   out_589892170851736179[11] = 0;
   out_589892170851736179[12] = 0;
   out_589892170851736179[13] = 1;
   out_589892170851736179[14] = 0;
   out_589892170851736179[15] = 0;
}
void h_30(double *state, double *unused, double *out_6752297139885956667) {
   out_6752297139885956667[0] = state[4];
}
void H_30(double *state, double *unused, double *out_1455446809497690677) {
   out_1455446809497690677[0] = 0;
   out_1455446809497690677[1] = 0;
   out_1455446809497690677[2] = 0;
   out_1455446809497690677[3] = 0;
   out_1455446809497690677[4] = 1;
   out_1455446809497690677[5] = 0;
   out_1455446809497690677[6] = 0;
   out_1455446809497690677[7] = 0;
}
void h_26(double *state, double *unused, double *out_4875757984572891982) {
   out_4875757984572891982[0] = state[7];
}
void H_26(double *state, double *unused, double *out_3055641902892824233) {
   out_3055641902892824233[0] = 0;
   out_3055641902892824233[1] = 0;
   out_3055641902892824233[2] = 0;
   out_3055641902892824233[3] = 0;
   out_3055641902892824233[4] = 0;
   out_3055641902892824233[5] = 0;
   out_3055641902892824233[6] = 0;
   out_3055641902892824233[7] = 1;
}
void h_27(double *state, double *unused, double *out_4930031343789200054) {
   out_4930031343789200054[0] = state[3];
}
void H_27(double *state, double *unused, double *out_3944297225607655941) {
   out_3944297225607655941[0] = 0;
   out_3944297225607655941[1] = 0;
   out_3944297225607655941[2] = 0;
   out_3944297225607655941[3] = 1;
   out_3944297225607655941[4] = 0;
   out_3944297225607655941[5] = 0;
   out_3944297225607655941[6] = 0;
   out_3944297225607655941[7] = 0;
}
void h_29(double *state, double *unused, double *out_5961981834362151480) {
   out_5961981834362151480[0] = state[1];
}
void H_29(double *state, double *unused, double *out_2969857451529201489) {
   out_2969857451529201489[0] = 0;
   out_2969857451529201489[1] = 1;
   out_2969857451529201489[2] = 0;
   out_2969857451529201489[3] = 0;
   out_2969857451529201489[4] = 0;
   out_2969857451529201489[5] = 0;
   out_2969857451529201489[6] = 0;
   out_2969857451529201489[7] = 0;
}
void h_28(double *state, double *unused, double *out_8205368988943805089) {
   out_8205368988943805089[0] = state[5];
   out_8205368988943805089[1] = state[6];
}
void H_28(double *state, double *unused, double *out_3045396279621735753) {
   out_3045396279621735753[0] = 0;
   out_3045396279621735753[1] = 0;
   out_3045396279621735753[2] = 0;
   out_3045396279621735753[3] = 0;
   out_3045396279621735753[4] = 0;
   out_3045396279621735753[5] = 1;
   out_3045396279621735753[6] = 0;
   out_3045396279621735753[7] = 0;
   out_3045396279621735753[8] = 0;
   out_3045396279621735753[9] = 0;
   out_3045396279621735753[10] = 0;
   out_3045396279621735753[11] = 0;
   out_3045396279621735753[12] = 0;
   out_3045396279621735753[13] = 0;
   out_3045396279621735753[14] = 1;
   out_3045396279621735753[15] = 0;
}
}

extern "C"{
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
}

#include <eigen3/Eigen/Dense>
#include <iostream>

typedef Eigen::Matrix<double, DIM, DIM, Eigen::RowMajor> DDM;
typedef Eigen::Matrix<double, EDIM, EDIM, Eigen::RowMajor> EEM;
typedef Eigen::Matrix<double, DIM, EDIM, Eigen::RowMajor> DEM;

void predict(double *in_x, double *in_P, double *in_Q, double dt) {
  typedef Eigen::Matrix<double, MEDIM, MEDIM, Eigen::RowMajor> RRM;
  
  double nx[DIM] = {0};
  double in_F[EDIM*EDIM] = {0};

  // functions from sympy
  f_fun(in_x, dt, nx);
  F_fun(in_x, dt, in_F);


  EEM F(in_F);
  EEM P(in_P);
  EEM Q(in_Q);

  RRM F_main = F.topLeftCorner(MEDIM, MEDIM);
  P.topLeftCorner(MEDIM, MEDIM) = (F_main * P.topLeftCorner(MEDIM, MEDIM)) * F_main.transpose();
  P.topRightCorner(MEDIM, EDIM - MEDIM) = F_main * P.topRightCorner(MEDIM, EDIM - MEDIM);
  P.bottomLeftCorner(EDIM - MEDIM, MEDIM) = P.bottomLeftCorner(EDIM - MEDIM, MEDIM) * F_main.transpose();

  P = P + dt*Q;

  // copy out state
  memcpy(in_x, nx, DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
}

// note: extra_args dim only correct when null space projecting
// otherwise 1
template <int ZDIM, int EADIM, bool MAHA_TEST>
void update(double *in_x, double *in_P, Hfun h_fun, Hfun H_fun, Hfun Hea_fun, double *in_z, double *in_R, double *in_ea, double MAHA_THRESHOLD) {
  typedef Eigen::Matrix<double, ZDIM, ZDIM, Eigen::RowMajor> ZZM;
  typedef Eigen::Matrix<double, ZDIM, DIM, Eigen::RowMajor> ZDM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, EDIM, Eigen::RowMajor> XEM;
  //typedef Eigen::Matrix<double, EDIM, ZDIM, Eigen::RowMajor> EZM;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1> X1M;
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> XXM;

  double in_hx[ZDIM] = {0};
  double in_H[ZDIM * DIM] = {0};
  double in_H_mod[EDIM * DIM] = {0};
  double delta_x[EDIM] = {0};
  double x_new[DIM] = {0};


  // state x, P
  Eigen::Matrix<double, ZDIM, 1> z(in_z);
  EEM P(in_P);
  ZZM pre_R(in_R);
  
  // functions from sympy
  h_fun(in_x, in_ea, in_hx);
  H_fun(in_x, in_ea, in_H);
  ZDM pre_H(in_H); 
  
  // get y (y = z - hx)
  Eigen::Matrix<double, ZDIM, 1> pre_y(in_hx); pre_y = z - pre_y;
  X1M y; XXM H; XXM R;
  if (Hea_fun){
    typedef Eigen::Matrix<double, ZDIM, EADIM, Eigen::RowMajor> ZAM;
    double in_Hea[ZDIM * EADIM] = {0};
    Hea_fun(in_x, in_ea, in_Hea);
    ZAM Hea(in_Hea);
    XXM A = Hea.transpose().fullPivLu().kernel();
   

    y = A.transpose() * pre_y;
    H = A.transpose() * pre_H;
    R = A.transpose() * pre_R * A;
  } else {
    y = pre_y;
    H = pre_H;
    R = pre_R;
  }
  // get modified H
  H_mod_fun(in_x, in_H_mod);
  DEM H_mod(in_H_mod);
  XEM H_err = H * H_mod;
  
  // Do mahalobis distance test
  if (MAHA_TEST){
    XXM a = (H_err * P * H_err.transpose() + R).inverse();
    double maha_dist = y.transpose() * a * y;
    if (maha_dist > MAHA_THRESHOLD){
      R = 1.0e16 * R;
    }
  }

  // Outlier resilient weighting
  double weight = 1;//(1.5)/(1 + y.squaredNorm()/R.sum());

  // kalman gains and I_KH
  XXM S = ((H_err * P) * H_err.transpose()) + R/weight;
  XEM KT = S.fullPivLu().solve(H_err * P.transpose());
  //EZM K = KT.transpose(); TODO: WHY DOES THIS NOT COMPILE?
  //EZM K = S.fullPivLu().solve(H_err * P.transpose()).transpose();
  //std::cout << "Here is the matrix rot:\n" << K << std::endl;
  EEM I_KH = Eigen::Matrix<double, EDIM, EDIM>::Identity() - (KT.transpose() * H_err);

  // update state by injecting dx
  Eigen::Matrix<double, EDIM, 1> dx(delta_x);
  dx  = (KT.transpose() * y);
  memcpy(delta_x, dx.data(), EDIM * sizeof(double));
  err_fun(in_x, delta_x, x_new);
  Eigen::Matrix<double, DIM, 1> x(x_new);
 
  // update cov 
  P = ((I_KH * P) * I_KH.transpose()) + ((KT.transpose() * R) * KT);

  // copy out state
  memcpy(in_x, x.data(), DIM * sizeof(double));
  memcpy(in_P, P.data(), EDIM * EDIM * sizeof(double));
  memcpy(in_z, y.data(), y.rows() * sizeof(double));
}



extern "C"{

      void update_25(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_25, H_25, NULL, in_z, in_R, in_ea, MAHA_THRESH_25);
      }
    
      void update_24(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_24, H_24, NULL, in_z, in_R, in_ea, MAHA_THRESH_24);
      }
    
      void update_30(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_30, H_30, NULL, in_z, in_R, in_ea, MAHA_THRESH_30);
      }
    
      void update_26(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_26, H_26, NULL, in_z, in_R, in_ea, MAHA_THRESH_26);
      }
    
      void update_27(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_27, H_27, NULL, in_z, in_R, in_ea, MAHA_THRESH_27);
      }
    
      void update_29(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<1,3,0>(in_x, in_P, h_29, H_29, NULL, in_z, in_R, in_ea, MAHA_THRESH_29);
      }
    
      void update_28(double *in_x, double *in_P, double *in_z, double *in_R, double *in_ea) {
        update<2,3,0>(in_x, in_P, h_28, H_28, NULL, in_z, in_R, in_ea, MAHA_THRESH_28);
      }
    
}
