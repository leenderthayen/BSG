#ifndef MATRIXELEMENTS
#define MATRIXELEMENTS

#include "gsl/gsl_sf_coupling.h"
#include "NilssonOrbits.h"

#include <cmath>

namespace MatrixElements {

inline int gL(int k) {
  return (k > 0) ? k : std::abs(k)-1;
}

inline int kL(int l, int s) {
  return (s == 1) ? -(l+1) : l;
}

inline double sign(double x) {
  return (x == 0) ? 0 : x/std::abs(x);
}

inline double getGKLs(int kf, int ki, double Ji, int K, int L, int s, int li, int lf, int si, int sf) {
  //cout << "Calc GKLs " << kf << " " << ki << endl;
  double first = std::sqrt((2*s+1)*(2*K+1)*(2*gL(kf)+1)*(2*gL(ki)+1)*(2*(lf+sf/2.)+1)*(2*(li+si/2.)+1));
  //double second = std::pow(sqrt(-1.), gL(ki)+gL(kf)+L)*std::pow(-1, li+si/2.-lf-sf/2.);
  double second = 1.;
  double third = gsl_sf_coupling_3j(2*gL(kf), 2*gL(ki), 2*L, 0, 0, 0)*std::sqrt(2*L+1.)*std::pow(-1, gL(kf)-gL(ki)-1);
  double fourth = gsl_sf_coupling_3j(2*gL(kf), 2*gL(ki), 2*L, 0, 0, 0)*std::sqrt(2*L+1.)*std::pow(-1, gL(kf)-gL(ki)-1);

  //cout << first << " " << second << " " << third << " " << fourth << endl;

  return first*second*third*fourth;
  /*return std::sqrt((2*s+1)*(2*K+1)*(2*gL(kf)+1)*(2*gL(ki)+1)*
                (2*(lf+sf/2.)+1)*(2*(li+si/2.)+1))*std::pow(sqrt(-1.), gL(ki)+gL(kf)+L)*std::pow(-1, li+si/2.-lf-sf/2.)
                *gsl_sf_coupling_3j(2*gL(kf), 2*gL(ki), 2*L, 0, 0, 0)*std::sqrt(2*L+1.)*std::pow(-1, gL(kf)-gL(ki)-1)*
                gsl_sf_coupling_9j(2*K, 2*s, 2*L, 2*lf+sf, 1, gL(kf), 2*li+si, 1, gL(ki));*/
}

inline double GetSingleParticleMatrixElement(bool V, double Ji, int K, int L, int s, int li, int lf, int si, int sf, double R) {
  //cout << "Calculating SP ME " << K << " " << L << " " << s << " state: " << li << " " << si << " " << lf << " " << sf << endl;
  double Mn = (protonMasskeV+neutronMasskeV)/2./electronMasskeV;
  double result = std::sqrt(2/(2*Ji+1));
  if (V) {
    if (s == 0) {
      result *= getGKLs(kL(lf, sf), kL(li, si), Ji, K, L, s, li, lf, si, sf);
    }
    else if (s == 1) {
      if (li == lf) {
        if (si == sf) {
          if (si == 1) {
            result *= -(li+1)*std::sqrt(6.*(li+1)*(2.*li+3)/(2.*li+1.))/(2.*Mn*R);
          }
          else {
            result *= -li*std::sqrt(6.*li*(2*li-1.)/(2.*li+1.))/(2.*Mn*R);
          }
        }
        else {
          result *= sf*std::sqrt(3.*li*(li+1.)/2./(2.*li+1))/(Mn*R);
        }
      }
      else {
        result = 0.;
      }
      /*result *= 1./((protonMasskeV+neutronMasskeV)/electronMasskeV*R)*(sign(kL(li, si))*getGKLs(kL(lf, sf), -kL(li, si), Ji, K, L, s, li, lf, si, sf)
                + sign(kL(lf, sf))*getGKLs(-kL(lf, sf), kL(li, si), Ji, K, L, s, li, lf, si, sf));
      result *= 3./4.;*/
    }
  }
  else {
    //result *= getGKLs(kL(lf, sf), kL(li, si), Ji, K, L, s, li, lf, si, sf);
    if (li == lf) {
      if (L == 0) {
        if (si == sf) {
          if (si == 1) {
            result *= std::sqrt((li+1.)*(2.*li+3.)/(2.*li+1));
          }
          else {
            result *= -std::sqrt(li*(2.*li-1.)/(2.*li+1.));
          }
        }
        else {
          result *= si*2.*std::sqrt(li*(li+1.)/(2.*li+1.));
        }
      }
      else if (L == 2) {
        if (si == sf) {
          if (si == 1) {
            result *= -li*std::sqrt(2.*(li+1.)/(2.*li+1)/(2.*li+3.));
          }
          else {
            result *= (li+1.)*std::sqrt(li/(2.*li+1)/(2.*li-1.));
          }
        }
        else {
          result *= si*std::sqrt(li*(li+1.)/2./(2.*li+1.));
        }
      }
    }
    else {
      result = 0.;
    }
  }
  //cout << "Result: " << result << endl;
  return result;
}

inline int delta(double x, double y) {
  return (x == y) ? 1 : 0;
}

inline double CalculateDeformedSPMatrixElement(nilsson::SingleParticleState spsi, nilsson::SingleParticleState spsf, bool V, int K, int L, int s, double Ji, double Jf, double Ki, double Kf, double R) {
  double result = 0.;

  //cout << "Calculating deformed SP ME " << K << " " << L << " " << s << endl;
  std::vector<nilsson::WFComp> finalStates = spsf.componentsHO;
  std::vector<nilsson::WFComp> initStates = spsi.componentsHO;
  int inO = 2 * spsi.O;
  int fO = 2 * spsf.O;
  for (std::vector<nilsson::WFComp>::iterator fIt = finalStates.begin(); fIt != finalStates.end(); ++fIt) {
    for(std::vector<nilsson::WFComp>::iterator inIt = initStates.begin(); inIt != initStates.end(); ++inIt) {
      result += (*fIt).C*(*inIt).C*(pow(-1, Jf-Kf+(*fIt).l+(*fIt).s/2.-fO/2.)*
	gsl_sf_coupling_3j(2*Jf, 2*K, 2*Ji, -2*Kf, fO-inO, 2*Ki)*
	gsl_sf_coupling_3j(2*(*fIt).l+(*fIt).s, 2*K, 2*(*inIt).l+(*inIt).s, -fO, fO-inO, inO)
	+ gsl_sf_coupling_3j(2*Jf, 2*K, 2*Ji, 2*Kf, -fO-inO, 2*Ki)*
	gsl_sf_coupling_3j(2*(*fIt).l+(*fIt).s, 2*K, 2*(*inIt).l+(*inIt).s, fO, -fO-inO, inO))*
	GetSingleParticleMatrixElement(V, Ji, K, L, s, (*inIt).l, (*fIt).l, (*inIt).s, (*fIt).s, R);
    }
  }

  result *= std::sqrt((2*Ji+1)*(2*Jf+1)/(1.+delta(Kf, 0.))/(1.+delta(Ki, 0.)));

  //cout << "Result: " << result << endl;

  return result;
}
}
#endif
