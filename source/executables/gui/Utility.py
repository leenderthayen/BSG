# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:14:14 2017

@author: leendert
"""

import numpy as np
import logging
import ConfigParser

''' General variables and utility functions '''

atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
         'Sc', 'Ti',
         'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb',
         'Mo', 'Tc', 'Ru', 'Rh',
         'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu',
         'Gd', 'Tb', 'Dy', 'Ho', 'Er',
         'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
         'Fr', 'Ra', 'Ac', 'Th', 'Pa',
         'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt',
         'Ds', 'Rg', 'Cn', 'Uut', 'Fl',
         'Uup', 'Lv', 'Uus', 'Uuo']


def convertFloat(s):
    """

    :param s: string representing a number

    """
    if s.strip() == '':
        return 0.
    else:
        try:
            return float(s)
        except ValueError:
            return s

def interpretENSDFSpinParities(mJpi, dJpi, preference=0):
    """Determine what to do with spin parities when it's not crystal clear

    :param mJpi: string, mother spin parities
    :param dJpi: string, daughter spin parities
    :param preference: int, when faced with multiple choices, which one to pick (index)

    Look at ensdf-manual on p52, section on J to properly interpret spins (Default value = 0)

    """

    """If either one is unknown, assume same spin as other one, i.e. allowed"""
    newMJpi = mJpi
    newDJpi = dJpi

    replacementList = ('(', ')', '[', ']', 'LE', 'GT', 'GE', ' ', '|>', 'NOT', '>', 'HIGHJ')

    for r in replacementList:
        newMJpi = newMJpi.replace(r, '')
        newDJpi = newDJpi.replace(r, '')

    #print(newMJpi, newDJpi)

    if newMJpi in ('', '-', '+'):
        newMJpi = newDJpi
    elif newDJpi in ('', '-', '+'):
        newDJpi = newMJpi

    newMJpi = newMJpi.replace(':', 'TO')
    newDJpi = newDJpi.replace(':', 'TO')

    newMJpi = newMJpi.replace('AND', 'TO')
    newDJpi = newDJpi.replace('AND', 'TO')

    indexM = min(len(newMJpi.split('TO'))-1, preference)
    indexD = min(len(newDJpi.split('TO'))-1, preference)

    newMJpi = newMJpi.split('TO')[indexM]
    newDJpi = newDJpi.split('TO')[indexD]

    indexM = min(len(newMJpi.split(','))-1, preference)
    indexD = min(len(newDJpi.split(','))-1, preference)
    newMJpi = newMJpi.split(',')[indexM]
    newDJpi = newDJpi.split(',')[indexD]

    if newMJpi == '' and newDJpi == '':
        newMJpi = '1+'
        newDJpi = '0+'
    try:
        if newMJpi[-1] not in ('-', '+'):
            if newDJpi[-1] in ('-', '+'):
                newMJpi += newDJpi[-1]
            else:
                newMJpi += '+'
        if newDJpi[-1] not in ('-', '+'):
            if newMJpi[-1] in ('-', '+'):
                newDJpi += newMJpi[-1]
            else:
                newDJpi += '+'
    except TypeError:
        print(newMJpi, newDJpi)
        print(mJpi, dJpi)
        raise

    replaced = False

    if newMJpi != mJpi or newDJpi != dJpi:
        #debug('Spin replacement: initial %s to %s; now %s to %s' % (mJpi, dJpi, newMJpi, newDJpi))
        replaced = True

    return (newMJpi, newDJpi, replaced)

class Level:
    """Class representing a nuclear level"""
    def __init__(self, E, dE, Jpi, deformation=None):
        """Initialize

        :param E: energy of nuclear level relative to ground state in keV
        :param dE: uncertainty on energy in keV
        :param Jpi: string, spin and parity
        :param deformation: dictionary containing beta deformations (Default value: None)

        """
        self.E = E
        self.dE = dE
        self.Jpi = Jpi
        self.deformation = deformation

class BetaBranch:
    """Independent container class for one particular beta branch,
    able to retrieve the spectrum shape


    """
    def __init__(self, process, E, dE, IB, dIB, logft, motherLevel, daughterLevel, partialHalflife):
        """Initialize a new beta branch

        :param process: Beta decay process (B-, B+)
        :param E: Endpoint energy in keV
        :param dE: Error on energy in keV
        :param IB: Intensity in %
        :param dIB: Error in intensity
        :param motherLevel: Level object of mother state
        :param daughterLevel: Level object of daughter state

        """
        self.process = process
        self.E = E
        self.dE = dE
        self.IB = IB
        self.dIB = dIB

        try:
            motherLevel.Jpi, daughterLevel.Jpi, self.spinReplacement = interpretENSDFSpinParities(motherLevel.Jpi, daughterLevel.Jpi)
        except:
            #TODO Change this!
            motherLevel = Level(0., 0., '0+')
            daughterLevel = Level(0., 0., '0+')

        self.motherLevel = motherLevel
        self.daughterLevel = daughterLevel

        self.logft = logft

        self.partialHalflife = partialHalflife


''' Data loaders '''


def removeDuplicateLevels(levels):
    """

    :param levels: list of levels from pyne.ensdf

    """
    removals = []
    for i in range(len(levels)):
        for j in range(i + 1, len(levels)):
            if levels[i][0] == levels[j][0]:
                removals.append(j)
    newLevels = [levels[i] for i in range(len(levels)) if i not in removals]
    return newLevels


def getENSDFBetaBranches(filename, Z, A, state):
    """

    :param filename: ensdf filename
    :param Z: proton number of mother
    :param A: mass number of other
    :param state: number designating metastability (0 ground state, 1 1st metastable state)

    """
    from pyne import ensdf

    logger = logging.getLogger(__name__)

    daugtherNameBetaM = '%3s' % str(A) + atoms[Z].upper()
    daughterNameBetaP = '%3s' % str(A) + atoms[Z-2].upper()
    parentName = '%3s' % str(A) + atoms[Z - 1].upper()

    decays = ensdf.decays(filename)
    levels = ensdf.levels(filename)
    qvals = ensdf.qvalues(filename)

    levels = removeDuplicateLevels(levels)

    ZA = Z * 1000 + A
    parentLevel = None
    for i in range(len(levels)):
        '''ZA's in PyNE ensdf are multiplied by 1E4'''
        if levels[i][0] / 10000 == ZA:
            if state == 0:
                if levels[i][5] == 0:
                    parentLevel = levels[i]
                    break
            else:
                if levels[i][-2] == state:
                    parentLevel = levels[i]
                    break

    if parentLevel == None:
        parentLevel = [ZA * 10000 + state, ]
        motherLevel = Level(0., 0., '')
    else:
        motherLevel = Level(parentLevel[5], 0. if parentLevel[6] is None else parentLevel[6], parentLevel[2])

    try:
        qBetaM, dqBetaM = [q[1:] for q in qvals if q[0] / 10000 == ZA][0]
    except IndexError:
        print("IndexError on q value for %s Z%d A%d S%d %s" % (atoms[Z-1], Z, A, state, filename))
        return 0., list()

    if dqBetaM is None:
        dqBetaM = 0.

    qBetaP = 0.
    dqBetaP = 0.

    halfLife = 0.

    bbs = list()
    correctDecays = []
    for d in decays:
        if d[0] == parentLevel[0]:
            halfLife = d[3]
            correctDecays += d[15] + d[16]
            #endpoint energy
            if d[5]:
                if d[1] > d[0]:
                    qBetaM = d[5]
                else:
                    qBetaP = d[5]

    for cd in correctDecays:
        for l in levels:
            try:
                # Level is equal to daughter level from correct decays
                if l[0] == cd[1]:
                    q, dq = (qBetaM, dqBetaM) if l[0] > cd[0] else (qBetaP, dqBetaP)
                    name = daugtherNameBetaM if l[0] > cd[0] else daughterNameBetaP
                    BR, dBR, logft = cd[4:7]

                    if BR == None:
                        BR = 0.
                    partialHalflife = 0.
                    if BR > 0:
                        partialHalflife = halfLife * 100. / BR

                    daughterLevel = Level(l[5], 0. if l[6] is None else l[6], l[2])
                    try:
                        # Endpoint energy
                        if cd[2] is not None:
                            if not (abs(motherLevel.E + q - daughterLevel.E - cd[2]) / cd[2] < 0.2 or
                                    abs(motherLevel.E + q - daughterLevel.E - cd[2]) < max(60, 3 * dq)):
                                print("Endpoint doesn't match up for %s (%.1f KeV) -> %s (%.1f keV) Q: %.1f B: %.1f"
                                      % (parentName, motherLevel.E, name, daughterLevel.E,
                                         motherLevel.E + q - daughterLevel.E, cd[2]))
                    except TypeError:
                        print(filename, Z, A, state)
                        print(motherLevel.E, daughterLevel.E, q, cd)
                        raise
                    # Logft
                    if logft is None:
                        logft = 6.
                    bb = BetaBranch('B-' if l[0] > cd[0] else 'B+', motherLevel.E + q - daughterLevel.E,
                                    np.sqrt(motherLevel.dE ** 2.0 + dq ** 2.0 + daughterLevel.dE ** 2.0),
                                    BR, dBR, logft, motherLevel, daughterLevel, partialHalflife)
                    if BR is None:
                        print(
                                '%s (%.1f KeV) -> %s (%.1f keV) branching ratio not known.' % (
                            parentName, motherLevel.E, name, daughterLevel.E))
                    bbs.append(bb)
            except IndexError:
                print(filename, Z, A, state)
                print(l)
                print(cd)
                raise
    return bbs


def loadDeformationData(filename):
    """

    :param filename: Default value = 'FRDM2012.dat')

    """
    return np.genfromtxt(filename, usecols=(0, 2, 7, 8, 9, 10))


#
def loadDeformation(Z, A, deformations):
    """Load deformation data from Moeller et al., arXiv:1508.06294

    :param Z: Atomic number
    :param A: Mass number
    :param deformations:

    """
    deformation = np.zeros(4)
    try:
        index = np.where((deformations[:, 0] == Z) & (deformations[:, 1] == A))[0][0]
        deformation = deformations[index, 2:]
    except:
        pass
    return deformation


def loadChargeRadiiData(filename):
    """

    :param filename: Default value = 'nuclear_charge_radii.txt')

    """
    return np.genfromtxt(filename, skip_header=12, names=True, missing_values='-',
                         filling_values=0.0)


def loadAtomicMassData(filename='mass16.txt'):
    """

    :param filename: Default value = 'mass16.txt')

    """
    return np.genfromtxt(filename, skip_header=39)


#
def loadChargeRadius(Z, A, radii):
    """Load experimental or parametrically deduced charge radius
    Parametrisation from Bao et al., PHYSICAL REVIEW C 94, 064315 (2016)

    :param Z: Atomic number
    :param A: Mass number
    :param radii: dict of loaded charge radii

    """
    radius = 0.
    try:
        if 'A' in radii.dtype.names:
            possibleRadii = list(radii[(radii['Z'] == Z) & (radii['A'] == A)][0])[3:]
        else:
            possibleRadii = list(radii[(radii['Z'] == Z) & (radii['N'] == A - Z)][0])[2:]
        if possibleRadii[0] > 0:
            radius = possibleRadii[0]
        else:
            radius = possibleRadii[2]
    except:
        # print("IndexError: Z: {0} A: {1}".format(Z, A))
        radius = UF.getEltonNuclearRadius(A)
    return radius

def writeIniFile(Zm, Zd, A, Q, process, decayType, beta2m, beta4m, beta6m, beta2d, beta4d, \
beta6d, mRad, dRad, mJpi, dJpi, phl=None, logft=None, dE=0.0, mE=0.0, name=None, prefix='', **kwargs):

    if not name:
        name = '%sZ%d_A%d_Q%.0f.ini' % (prefix, Zm, A, Q)
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.optionxform = str
    config.add_section('Transition')
    config.add_section('Mother')
    config.add_section('Daughter')
    config.set('Transition', 'Process', process)
    config.set('Transition', 'Type', decayType)
    config.set('Transition', 'MixingRatio', 0.0)
    config.set('Transition', 'QValue', Q)
    if phl:
        config.set('Transition', 'PartialHalflife', phl)
    if logft:
        config.set('Transition', 'Logft', logft)
    config.set('Mother', 'Z', Zm)
    config.set('Mother', 'A', A)
    config.set('Mother', 'Radius', mRad)
    config.set('Mother', 'Beta2', beta2m)
    config.set('Mother', 'Beta4', beta4m)
    config.set('Mother', 'Beta6', beta6m)
    config.set('Mother', 'SpinParity', mJpi)
    config.set('Mother', 'ExcitationEnergy', mE)
    config.set('Daughter', 'Z', Zd)
    config.set('Daughter', 'A', A)
    config.set('Daughter', 'Radius', dRad)
    config.set('Daughter', 'Beta2', beta2d)
    config.set('Daughter', 'Beta4', beta4d)
    config.set('Daughter', 'Beta6', beta6d)
    config.set('Daughter', 'SpinParity', dJpi)
    config.set('Daughter', 'ExcitationEnergy', dE)
    with open(name, 'wb') as configFile:
        config.write(configFile)

def writeConfigFile(name, directory, computational, constants, spectrumCheckBoxes,\
spectrumDSB, spectrumNME, enforceNME, spectrumComboBoxes):
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.optionxform = str
    config.add_section('Spectrum')
    config.add_section('General')
    config.add_section('Computational')
    config.add_section('Constants')
    config.set('General', 'Folder', directory)
    for key in computational:
        try:
            config.set('Computational', computational[key], key.value())
        except:
            try:
                config.set('Computational', computational[key], key.isChecked())
            except:
                try:
                    config.set('Computational', computational[key], key.currentText())
                except:
                    pass
    for key in constants:
        config.set('Constants', constants[key], key.value())

    for key in spectrumCheckBoxes:
        config.set('Spectrum', spectrumCheckBoxes[key], key.isChecked())
    for key in spectrumDSB:
        config.set('Spectrum', spectrumDSB[key], key.value())
    if enforceNME:
        for key in spectrumNME:
            config.set('Spectrum', spectrumNME[key], key.value())
    for key in spectrumComboBoxes:
        config.set('Spectrum', spectrumComboBoxes[key], key.currentText())
    with open(name, 'wb') as configFile:
        config.write(configFile)


#
# def writeIniFile(Z, A, Q, process, decayType, beta2m, beta4m, beta6m, beta2d, beta4d,
#                  beta6d, mRad, dRad, mJpi, dJpi, mE=0.0, dE=0.0, name=None, prefix='', **kwargs):
#     """
#
#     :param Z: proton number of mother
#     :param A: mass number of mother
#     :param Q: Q value in keV
#     :param process: B-/B+/EC
#     :param decayType: Fermi/Gamow-Teller
#     :param beta2m: mother quadrupole deformation
#     :param beta4m: mother beta4 deformation
#     :param beta6m: mother beta6 deformation
#     :param beta2d: daughter quadrupole deformation
#     :param beta4d: daughter beta4 deformation
#     :param beta6d: daughter beta6 deformation
#     :param mRad: mother charge radius (fm)
#     :param dRad: daughter charge radius (fm)
#     :param mJpi: mother double of spin*parity
#     :param dJpi: daughter double of spin*parity
#     :param mE: mother excitation energy in keV (Default value = 0.0)
#     :param dE: daughter excitation energy in keV (Default value = 0.0)
#     :param name: optional filename (Default value = None)
#     :param prefix: optional prefix to filename (Default value = '')
#     :param **kwargs:
#
#     """
#
#     if not name:
#         name = '%sZ%d_A%d_Q%.0f.ini' % (prefix, Z, A, Q)
#     with open(name, 'w') as f:
#         text = '[Transition]\nProcess=%s\nType=%s\nMixingRatio=0.0\nQValue=%f\n\n' \
#                % (process, decayType, Q)
#         f.write(text)
#         text = '[Mother]\nZ=%d\nA=%d\nRadius=%f\nBeta2=%f\nBeta4=%f\nBeta6=%f\nSpinParity=%d\nExcitationEnergy=%f\n' \
#                % (Z, A, mRad, beta2m, beta4m, beta6m, mJpi, mE)
#         f.write(text)
#         text = '[Daughter]\nZ=%d\nA=%d\nRadius=%f\nBeta2=%f\nBeta4=%f\nBeta6=%f\nSpinParity=%d\nExcitationEnergy=%f\n' \
#                % (Z + 1 if process == 'B-' else Z - 1, A, dRad, beta2d, beta4d, beta6d, dJpi, dE)
#         f.write(text)
#     return name
#
#
# # @profile
# def createIniFileFromBetaInfo(betaInfo, radii, deformations, prefix='fission_'):
#     """
#
#     :param betaInfo: dict containing transition info
#     :param radii: dict containing pre-loaded charge radii
#     :param deformations: dict containing pre-loaded nuclear deformation
#     :param prefix: filename prefix (Default value = 'fission_')
#
#     """
#
#     if not 'beta2m' in betaInfo:
#         mDef = loadDeformation(betaInfo['Z'], betaInfo['A'], deformations)
#         betaInfo['beta2m'] = mDef[0]
#         betaInfo['beta4m'] = mDef[1]
#         betaInfo['beta6m'] = mDef[2]
#     if not 'beta2d' in betaInfo:
#         dDef = loadDeformation(betaInfo['Z'] + 1 if betaInfo['process'] == 'B-' else betaInfo['Z'] - 1,
#                                betaInfo['A'], deformations)
#         betaInfo['beta2d'] = dDef[0]
#         betaInfo['beta4d'] = dDef[1]
#         betaInfo['beta6d'] = dDef[2]
#     if not 'mJpi' in betaInfo:
#         betaInfo['mJpi'] = 2
#     if not 'dJpi' in betaInfo:
#         betaInfo['dJpi'] = 0
#     betaInfo['mRad'] = loadChargeRadius(betaInfo['Z'], betaInfo['A'], radii)
#     betaInfo['dRad'] = loadChargeRadius(betaInfo['Z'] + 1 if betaInfo['process'] == 'B-' else betaInfo['Z'] - 1,
#                                         betaInfo['A'], radii)
#
#     betaInfo['prefix'] = prefix
#
#     properBetaInfo = dict(betaInfo)
#     properBetaInfo['Q'] = betaInfo['E0'] + betaInfo['dE'] - betaInfo['mE']
#     del properBetaInfo['E0']
#
#     return writeIniFile(**properBetaInfo)
#
#
# def createIniFile(Z, A, Q, process='B-', decayType='Gamow-Teller', prefix='fission_'):
#     """Generate INI file for beta spectrum generator program
#
#     :param Z: Atomic number of parent
#     :param A: Mass number of parent
#     :param Q: Q value
#     :param process: string, distinction between beta+/- or electron capture (Default value = 'B-')
#     :param decayType: string, Fermi or Gamow-Teller (Default value = 'Gamow-Teller')
#     :param prefix: prefix for INI filename (Default value = 'fission_')
#
#     """
#     Zd = Z + 1
#     if process == "B+" or process == "EC":
#         Zd = Z - 1
#     motherDeformation = loadDeformation(Z, A)
#     daughterDeformation = loadDeformation(Zd, A)
#     motherRadius = loadChargeRadius(Z, A)
#     daughterRadius = loadChargeRadius(Zd, A)
#
#     return writeIniFile(Z, A, Q, process, decayType, motherDeformation[0], motherDeformation[1],
#                         motherDeformation[2], daughterDeformation[0], daughterDeformation[1], daughterDeformation[2],
#                         motherRadius, daughterRadius, 2, 0, prefix=prefix)


class DataManager:
    """Responsible for pre-loading nuclear databases for charge radii, deformation & masses"""

    def __init__(self, name, iniPrefix, chargeRadiiFile, deformFile):
        self.name = name
        self.iniPrefix = iniPrefix
        self.radiiData = None
        self.deformationData = None

        logger = logging.getLogger(__name__)

        try:
            self.radiiData = loadChargeRadiiData(chargeRadiiFile)
            self.deformationData = loadDeformationData(deformFile)
        except OSError:
            logger.warn("Radii or deformation not found.")
            pass

        # self.massData = loadAtomicMassData()

    def createIniFileFromBetaInfo(self, betaInfo):
        """

        :param betaInfo:

        """
        return createIniFileFromBetaInfo(betaInfo, self.radiiData, self.deformationData, self.iniPrefix)

    def loadDeformation(self, Z, A):
        """

        :param Z: param A:
        :param A:

        """
        return loadDeformation(Z, A, self.deformationData)

    def loadChargeRadii(self, Z, A):
        """

        :param Z: param A:
        :param A:

        """
        return loadChargeRadius(Z, A, self.radiiData)
