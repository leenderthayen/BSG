#ifndef CONSTANTS_H
#define CONSTANTS_H

namespace bsg {

const double ELECTRON_MASS_KEV = 510.998902;  /**< the electron rest mass in keV */
const double NEUTRON_MASS_KEV = 939565.378;   /**< the neutron rest mass in keV */
const double PROTON_MASS_KEV = 938271.998;    /**< the proton rest mass in keV */
const double NUCLEON_MASS_KEV = (NEUTRON_MASS_KEV + PROTON_MASS_KEV) / 2.; /**< the average nucleus mass in keV */
const double HBAR = 6.58211889e-16;       /**< reduced Planck's constant in units eV*s */
const double SPEED_OF_LIGHT = 299792458.0;  /**< the speed of light in m/s */
const double PION_MASS_KEV = 134976.6; /**< the pion rest mass in unuts in keV */

const double NATURAL_LENGTH = HBAR * SPEED_OF_LIGHT / 1e3 / ELECTRON_MASS_KEV; /**< the beta decay natural length scale using HBAR=c=m_e=1 3.861592643659598e-13 m*/

const double ALPHA = 0.0072973525664; /**< the fine-structure constant */
const double EULER_MASCHERONI_CONSTANT = 0.577215664901532;  /**< the Euler-Mascheroni constant */

#define sqr(x) ((x) * (x))

}

#endif  // CONSTANTS_H
