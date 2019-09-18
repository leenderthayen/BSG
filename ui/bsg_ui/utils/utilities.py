# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:14:14 2017

@author: leendert
"""

import numpy as np
import logging
import configparser

''' General variables and utility functions '''

atoms = ['Nn', 'H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca',
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


def getEltonNuclearRadius(A, natUnits=False):
    """Return nuclear charge radius in fm according to Elton formula

    :param A: mass number
    :param natUnits: return result in units where m_e=hbar=c=1 (Default value = False)

    """
    R = 1.121 * A ** (1. / 3.) + 2.426 * A ** (-1. / 3.) - 6.614 / A

    if natUnits:
        """ hbar * c / m_e, in fm"""
        natLength = 3.861592643659598e2
        R /= natLength
    return R

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
    import ensdf
    logger = logging.getLogger(__name__)

    daughterNameBetaM = '%3s' % str(A) + atoms[Z+1].upper()
    daughterNameBetaP = '%3s' % str(A) + atoms[Z-1].upper()
    parentName = '%3s' % str(A) + atoms[Z].upper()

    decays = ensdf.decays(filename)
    levels = ensdf.levels(filename)
    qvals = ensdf.qvalues(filename)


    levels = removeDuplicateLevels(levels)

    ZA = Z * 1000 + A
    parentLevel = None
    for i in range(len(levels)):
        '''Level listing allows for 999 excited states'''
        if levels[i][0] / 1000 == ZA:
            if state == 0:
                if levels[i][5] == 0:
                    parentLevel = levels[i]
                    break
            else:
                if levels[i][-2] == state:
                    parentLevel = levels[i]
                    break

    if parentLevel == None:
        parentLevel = [ZA * 1000 + state, ]
        motherLevel = Level(0., 0., '')
    else:
        motherLevel = Level(parentLevel[5], 0. if parentLevel[6] is None else parentLevel[6], parentLevel[2])

    try:
        qBetaM, dqBetaM = [q[1:] for q in qvals if q[0] / 1000 == ZA][0]
    except IndexError:
        print("IndexError on q value for %s Z%d A%d S%d %s" % (atoms[Z], Z, A, state, filename))
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
                    name = daughterNameBetaM if l[0] > cd[0] else daughterNameBetaP
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


def loadDeformation(Z, A, filename = 'FRDM2012/FRDM2012.dat'):
    """Load deformation data from Moeller et al., arXiv:1508.06294
    :param Z: Atomic number
    :param A: Mass number
    :param filename: location and name to data file (Default value = 'FRDM2012/FRDM2012.dat')
    """

    deformation = np.zeros(4)
    with open(filename, 'r') as f:
        for line in f:
            d = np.array([float(s) for s in line.split()])
            if d[0] == Z and d[2] == A:
                deformation = d[7:11]
                break
    return deformation

def loadChargeRadius(Z, A, filename = 'ChargeRadii/nuclear_charge_radii.txt'):
    """Load experimental or parametrically deduced charge radius
    Parametrisation from Bao et al., PHYSICAL REVIEW C 94, 064315 (2016)
    :param Z: Atomic number
    :param A: Mass number
    :param filename: filename for data (Default value = 'ChargeRadii/nuclear_charge_radii.txt')
    """
    radius = 0.
    radii = np.array(np.genfromtxt(filename, dtype=None, skip_header=10,
    missing_values='-', filling_values=0.0, names=True, autostrip=True))
    try:
        if 'A' in radii.dtype.names:
            possibleRadii = list(radii[(radii['Z'] == Z) & (radii['A'] == A)][0])[3:]
        else:
            possibleRadii = list(radii[(radii['Z'] == Z) & (radii['N'] == A-Z)][0])[2:]
        if possibleRadii[0] > 0:
            radius = possibleRadii[0]
        else:
            radius = possibleRadii[2]
    except IndexError:
        radius = getEltonNuclearRadius(A)
    return radius

def writeIniFile(Zm, Zd, A, Q, process, decayType, beta2m, beta4m, beta6m, beta2d, beta4d,
beta6d, mRad, dRad, mJpi, dJpi, phl=None, logft=None, dE=0.0, mE=0.0, name=None, prefix='',
robtdFile=None, **kwargs):

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
        config.set('Transition', 'LogFt', logft)
    if robtdFile:
        config.set('Transition', 'ROBTDFile', robtdFile)
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

def writeConfigFile(name, config_settings):
    config = configparser.ConfigParser()
    config.optionxform = str
    for conf_key in self.config_settings:
        config.add_section(conf_key)
        for key in self.config_settings[conf_key]:
            try:
                if isinstance(self.config_settings[conf_key][key],tuple):
                    value = self.config_settings[conf_key][key][0]
                else:
                    value = self.config_settings[conf_key][key]
                if isinstance(value,bool):
                    config.set(conf_key, key,'True' if value==True else 'False')
                elif isinstance(value,int):
                    config.set(conf_key, key, str(value))
                elif isinstance(value,float):
                    config.set(conf_key, key, str(value))
                else:
                    config.set(conf_key, key, value)
            except (configparser.NoSectionError, configparser.NoOptionError):
                continue

    with open(name, 'wb') as configFile:
        config.write(configFile)
