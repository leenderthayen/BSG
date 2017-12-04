#ifndef CONSTANTS_H
#define CONSTANTS_H

const double electronMasskeV = 510.998902;  /**< the electron rest mass in keV */
const double neutronMasskeV = 939565.378;   /**< the neutron rest mass in keV */
const double protonMasskeV = 938271.998;    /**< the proton rest mass in keV */
const double nucleonMasskeV = (neutronMasskeV + protonMasskeV) / 2.; /**< the average nucleus mass in keV */
const double hbar = 6.58211889e-16;       /**< reduced Planck's constant in units eV*s */
const double speedOfLight = 299792458.0;  /**< the speed of light in m/s */
const double pionMasskeV = 134976.6; /**< the pion rest mass in unuts in keV */

const double NATLENGTH = hbar * speedOfLight / 1e3 / electronMasskeV; /**< the beta decay natural length scale using hbar=c=m_e=1 */

const double alpha = 0.0072973525664; /**< the fine-structure constant */
const double EulMasConst = 0.577215664901532;  /**< the Euler-Mascheroni constant */

#define sqr(x) ((x) * (x))

#endif  // CONSTANTS_H
