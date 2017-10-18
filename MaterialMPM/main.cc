#include <cmath>

#include <boost/array.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/vector.hpp>

int main() {

  typedef boost::numeric::ublas::bounded_matrix<double, 6, 6> Mat6x6_;
  typedef boost::numeric::ublas::bounded_vector<double, 6> Vec6x1_;

  // Input
  double youngModulus_ = 10000000;
  double poissonRatio_ = 0.2;
  double phi_ = 30;
  double coh_ = 1000;
  double psi_ = 0;
  double sigt_ = 0;

  // \param[in]  F dstrain vector
  // \param[out] S stress vector
  Vec6x1_ F;
  F(0) = 0.00005;
  F(1) = -0.0001;
  F(2) = 0.00005;
  F(3) = 0;
  F(4) = 0;
  F(5) = 0;

  Vec6x1_ S;
  S(0) = 0;
  S(1) = -3464.1016;
  S(2) = 0;
  S(3) = 0;
  S(4) = 0;
  S(5) = 0;

  // Elastic stiffness matrix
  Mat6x6_ De;
  Vec6x1_ dS;
  // Bulk and shear modulus
  double K, G;
  double a1, a2;

  K = youngModulus_ / (3.0 * (1. - 2. * poissonRatio_));
  G = youngModulus_ / (2.0 * (1. + poissonRatio_));

  a1 = K + (4.0 / 3.0) * G;
  a2 = K - (2.0 / 3.0) * G;

  // compute elasticityTensor
  De(0,0)=a1;    De(0,1)=a2;    De(0,2)=a2;    De(0,3)=0;    De(0,4)=0;    De(0,5)=0;
  De(1,0)=a2;    De(1,1)=a1;    De(1,2)=a2;    De(1,3)=0;    De(1,4)=0;    De(1,5)=0;
  De(2,0)=a2;    De(2,1)=a2;    De(2,2)=a1;    De(2,3)=0;    De(2,4)=0;    De(2,5)=0;
  De(3,0)= 0;    De(3,1)= 0;    De(3,2)= 0;    De(3,3)=G;    De(3,4)=0;    De(3,5)=0;
  De(4,0)= 0;    De(4,1)= 0;    De(4,2)= 0;    De(4,3)=0;    De(4,4)=G;    De(4,5)=0;
  De(5,0)= 0;    De(5,1)= 0;    De(5,2)= 0;    De(5,3)=0;    De(5,4)=0;    De(5,5)=G;

  double fx, fy, fz, fxy;
  double mohr_c, mohr_r, f1, f2;
  double sigma1, sigma2, sigma3;
  unsigned mohr_flag;
  double nphi, npsi, f_s, f_t, alphap, heichi;
  double lamda_s, lamda_t;
  double cs2, si2, dc2, dss;
  double PI = std::atan(1.0) * 4.;
  double phi, psi;

  phi = phi_ * PI / 180.;
  psi = psi_ * PI / 180.;

  boost::numeric::ublas::noalias(dS) = boost::numeric::ublas::prod(De, F);

  S += dS;

  fx = S(0);
  fy = S(1);
  fz = S(2);
  fxy = S(3);

  mohr_c = 0.5 * (fx + fy);
  mohr_r = 0.5 * sqrt((fx - fy) * (fx - fy) + 4.0 * fxy * fxy);

  f1 = mohr_c - mohr_r;
  f2 = mohr_c + mohr_r;

  if (f1 > fz) {
    sigma1 = fz;
    sigma2 = f1;
    sigma3 = f2;
    mohr_flag = 2;
  } else if (f2 < fz) {
    sigma1 = f1;
    sigma2 = f2;
    sigma3 = fz;
    mohr_flag = 3;
  } else {
    sigma1 = f1;
    sigma2 = fz;
    sigma3 = f2;
    mohr_flag = 1;
  }

  // Mohr-Coulomb failure criteria
  nphi = (1.0 + sin(phi)) / (1.0 - sin(phi));
  npsi = (1.0 + sin(psi)) / (1.0 - sin(psi));
  f_s = sigma1 - sigma3 * nphi + 2.0 * coh_ * sqrt(nphi);
  f_t = sigt_ - sigma3;
  alphap = sqrt(1.0 + nphi * nphi) + nphi;
  heichi = sigma3 - sigt_ +
           alphap * (sigma1 - nphi * sigt_ + 2.0 * coh_ * sqrt(nphi));

  K = youngModulus_ / (3.0 * (1 - 2 * poissonRatio_));
  G = youngModulus_ / (2.0 * (1 + poissonRatio_));
  a1 = K + (4.0 / 3.0) * G;
  a2 = K - (2.0 / 3.0) * G;

  if (heichi < 0.) {
    if (f_s < 0.) {
      // tens + shear failure
      lamda_s = f_s / ((a1 - a2 * npsi) - (a2 - a1 * npsi) * nphi);
      sigma1 -= lamda_s * (a1 - a2 * npsi);
      sigma2 -= lamda_s * a2 * (1.0 - npsi);
      sigma3 -= lamda_s * (a2 - a1 * npsi);
    }
  } else if (heichi > 0.) {
    if (f_t < 0.) {
      // tens failure
      lamda_t = f_t / a1;
      sigma1 += lamda_t * a2;
      sigma2 += lamda_t * a2;
      sigma3 += lamda_t * a1;
      sigt_ = 0.0;
    }
  }

  if (sigma1 == sigma3) {
    cs2 = 1.0;
    si2 = 0.0;
  } else {
    cs2 = (fx - fy) / (f1 - f2);
    si2 = 2.0 * fxy / (f1 - f2);
  
}
  if (mohr_flag == 1) {
    dc2 = (sigma1 - sigma3) * cs2;
    dss = sigma1 + sigma3;
    S(0) = 0.5 * (dss + dc2);
    S(1) = 0.5 * (dss - dc2);
    S(2) = sigma2;
    S(3) = 0.5 * (sigma1 - sigma3) * si2;
    S(4) = 0.0;
    S(5) = 0.0;
  } else if (mohr_flag == 2) {
    dc2 = (sigma2 - sigma3) * cs2;
    dss = sigma2 + sigma3;
    S(0) = 0.5 * (dss + dc2);
    S(1) = 0.5 * (dss - dc2);
    S(2) = sigma1;
    S(3) = 0.5 * (sigma2 - sigma3) * si2;
    S(4) = 0.0;
    S(5) = 0.0;
  } else {
    dc2 = (sigma1 - sigma2) * cs2;
    dss = sigma1 + sigma2;
    S(0) = 0.5 * (dss + dc2);
    S(1) = 0.5 * (dss - dc2);
    S(2) = sigma3;
    S(3) = 0.5 * (sigma1 - sigma2) * si2;
    S(4) = 0.0;
    S(5) = 0.0;
  }

  std::cout << S(0) << '\t' << S(1) << '\t' << S(2) << '\t' 
            << S(3) << '\t' << S(4) << '\t' << S(5) << '\t'<< "\n";

  return 0;

}