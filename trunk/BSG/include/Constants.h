#ifndef CONSTANTS_H
#define CONSTANTS_H

const double electronMasskeV = 510.998902;  // in keV
const double neutronMasskeV = 939565.378;   // in keV
const double protonMasskeV = 938271.998;    // in keV
const double nucleonMasskeV = (neutronMasskeV + protonMasskeV) / 2.;
const double hbar = 6.58211889e-16;       // in ev*s
const double speedOfLight = 299792458.0;  // in m/s
const double pionMasskeV = 134976.6;

const double NATLENGTH = hbar * speedOfLight / 1e3 / electronMasskeV;

const double alpha = 0.0072973525664;
const double EulMasConst = 0.577215664901532;  // Euler-Mascheroni constant

// GSL_CONST_NUM_FINE_STRUCTURE

#define sqr(x) ((x) * (x))

#endif  // CONSTANTS_H
