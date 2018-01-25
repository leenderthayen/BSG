# -*- coding: utf-8 -*-
"""
Created on Fri Sep  1 11:14:14 2017

@author: leendert
"""

import numpy as np
import ConfigParser

''' General variables and utility functions '''

atoms = ['H', 'He', 'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne', 'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar', 'K', 'Ca', 'Sc', 'Ti',\
'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh',\
'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe', 'Cs', 'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', \
'Tm', 'Yb', 'Lu', 'Hf', 'Ta', 'W', 'Re', 'Os',  'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn', 'Fr', 'Ra', 'Ac', 'Th', 'Pa',\
'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr', 'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Fl',\
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
        
dbColumnNames = ['name', 'Z', 'A', 'E0', 'IB', 'mJpi', 'dJpi', 'dE', 'deltaJ', 'f', 'u', 'fissionYield',\
    'beta2m', 'beta3m', 'beta4m', 'beta6m', 'beta2d', 'beta3d', 'beta4d', 'beta6d']
    
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
    radii = np.array(np.genfromtxt(filename, dtype=None, skip_header=12, \
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
        print("IndexError: Z: {0} A: {1}".format(Z, A))
        #radius = UF.getEltonNuclearRadius(A)
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
        
def interpretSpinParities(mJpi, dJpi, preference=0):
    """Determine what to do with spin parities when it's not crystal clear

    :param mJpi: string, mother spin parities
    :param dJpi: string, daughter spin parities
    :param preference: int, when faced with multiple choices, which one to pick (index)
    
    Look at ensdf-manual op p52, section on J to properly interpret spins (Default value = 0)

    """
    
    """If either one is unknown, assume same spin as other one, i.e. allowed"""
    newMJpi = mJpi
    newDJpi = dJpi
    
    replacementList = ('(', ')', '[', ']', 'LE', 'GT', 'GE', ' ')
    
    for r in replacementList:
        newMJpi = newMJpi.replace(r, '')
        newDJpi = newDJpi.replace(r, '')
        
    print(newMJpi, newDJpi)
    
    if newMJpi in ('', '-', '+'):
        newMJpi = newDJpi
    elif newDJpi in ('', '-', '+'):
        newDJpi = newMJpi;
        
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
    def __init__(self, E, dE, IB, dIB, motherLevel, daughterLevel, process = 'B-'):
        """Initialize a new beta branch
        
        :param E: Endpoint energy in keV
        :param dE: Error on energy in keV
        :param IB: Intensity in %
        :param dIB: Error in intensity
        :param motherLevel: Level object of mother state
        :param daughterLevel: Level object of daughter state
        :param process: Beta process (B-, B+, EC)
        
        """
        self.E = E
        self.dE = dE
        self.IB = IB
        self.dIB = dIB
        
        motherLevel.Jpi, daughterLevel.Jpi, self.spinReplacement = interpretSpinParities(motherLevel.Jpi, daughterLevel.Jpi)
        
        self.motherLevel = motherLevel
        self.daughterLevel = daughterLevel
        
        self.process = process
        
class BetaCollection:
    """Container of BetaBranch objects for one specific (Z, A) isotope"""
    def __init__(self, name, Z, A):
        """Initialize
        
        :param name: name of mother isotope
        :param Z: Atomic number of mother
        :param A: Mass number
        
        """
        self.name = name
        self.Z = Z
        self.A = A
        self.branches = []
        
    def addBetaBranch(self, bb):
        """Add BetaBranch object to container

        :param bb: BetaBranch object

        """
        self.branches.append(bb)
        
def buildBetaBranches(filename, Z, A):
    """Parse ENSDF file to construct all individual beta- branches
    Returns a BetaCollection object containing all BetaBranch objects

    :param filename: ENSDF file name
    :param Z: Atomic number of the mother
    :param A: Mass number

    """
    name = '%3s' % str(A) + atoms[Z].upper()
    parentName = '%3s' % str(A) + atoms[Z-1].upper()
    
    #debug('Building beta branches from ' + parentName + ' to ' + name)
    
    bc = BetaCollection(str(A) + atoms[Z-1], Z, A)
    
    Q = 0.
    dQ = 0.
    
    foundDecay = False
    
    motherLevel = None
    daughterLevel = None
    
    # Have to loop over the entire file to find the B- branch
    with open(filename, 'r') as f:
        for line in f:
            if line.startswith(name) or (line.startswith(parentName) and line[7] == 'P'):
                mod = line[6]
                recType = line[7]
                if recType == ' ' and 'B- DECAY' in line:
                    foundDecay = True
                # If position of B- decay info is found, start parsing
                if foundDecay:
                    # Parent level
                    if mod == ' ' and recType == 'P':
                        pE, dpE, pJpi, _, _, Q, dQ = processENSDFLine(line)
                        pJpi = pJpi.strip()
                        motherLevel = Level(pE, dpE, pJpi)
                    # Daughter level
                    elif line[5] == ' ' and mod == ' ' and recType == 'L':
                        data = processENSDFLine(line.replace('+X', '  '))
                        dJpi = data[2].strip()
                        daughterLevel = Level(data[0], data[1], dJpi)
                    # Beta decay information level: Branching ratio, ft value
                    # Signals end of beta branch
                    elif mod == ' ' and recType == 'B' and line[5] != 'S':
                        data = processENSDFLine(line)
                        
                        try:
                            bb = BetaBranch(motherLevel.E+Q-daughterLevel.E, np.sqrt(motherLevel.dE**2.0+dQ**2.0+daughterLevel.dE**2.0), \
                            data[2], data[3], motherLevel, daughterLevel)
                            bc.addBetaBranch(bb)
                            
                            #debug('Creating beta branch with E0 = %.1f keV (%.5f' % (bb.E, bb.IB) + '%)')
                        except TypeError as e:
                            print "TypeError ({0}): {1}".format(e.args, e.message)        
                            #debug('ERROR: %d %d' % (Z, A))
            elif foundDecay:
                break
    return bc

''' Data loaders '''

''' Convert the FORTRAN style of decay branch files to proper output '''
def processENSDFLine(line):
    """Process a single line from an ENSDF data file, returns array of relevant data

    :param line: string, single line from ENSDF data

    """
    name = line[:5]
    recType = line[7]
    mod = line[6]
    
    data = []
    
    if mod == 'C' or mod == 'c':
        pass
    else:
        if recType == 'Q':
            Q = convertFloat(line[9:19])
            DQ = convertFloat(line[19:21])
            SN = convertFloat(line[21:29])
            DSN = convertFloat(line[29:31])
            SP = convertFloat(line[31:39])
            DSP = convertFloat(line[39:41])
            QA = convertFloat(line[41:49])
            DQA = convertFloat(line[49:55])
            data = [Q, DQ, SN, DSN, SP, DSP, QA, DQA]
        elif recType == 'P':
            E = convertFloat(line[9:19])
            DE = convertFloat(line[19:21])
            J = line[21:39]
            T = convertFloat(line[39:49])
            DT = convertFloat(line[49:55])
            QP = convertFloat(line[64:74])
            DQP = convertFloat(line[74:76])
            if line[74:76] == 'SY':
                DQP = 100.
            data = [E, DE, J, T, DT, QP, DQP]
        elif recType == 'N' and mod == '':
            NR = convertFloat(line[9:19])
            DNR = convertFloat(line[19:21])
            NT = convertFloat(line[21:29])
            DNT = convertFloat(line[29:31])
            BR = convertFloat(line[31:39])
            DBR = convertFloat(line[39:41])
            NB = convertFloat(line[41:49])
            DNB = convertFloat(line[49:55])
            NP = convertFloat(line[55:62])
            DNP = convertFloat(line[62:64])
            data = [NR, DNR, NT, DNT, BR, DBR, NB, DNB, NP, DNP]
        elif recType == 'L':
            E = convertFloat(line[9:19])
            DE = convertFloat(line[19:21])
            J = line[21:39]
            T = convertFloat(line[39:49])
            DT = convertFloat(line[49:55])
            L = convertFloat(line[55:64])
            S = convertFloat(line[64:74])
            DS = convertFloat(line[74:76])
            MS = line[77:79]
            Q = line[79]
            data = [E, DE, J, T, DT, L, S, DS, MS, Q]
        elif recType == 'B':
            E = convertFloat(line[9:19])
            DE = convertFloat(line[19:21])
            IB = convertFloat(line[21:29])
            DIB = convertFloat(line[29:31])
            LOGFT = convertFloat(line[41:49])
            DFT = convertFloat(line[49:55])
            UN = line[77:79]
            Q = line[79]
            data = [E, DE, IB, DIB, LOGFT, DFT, UN, Q]
        elif recType == 'B' and line[5] == 'S':
            s = (line[13:]).split()
            EAV = convertFloat(s[0])
            DEAV = convertFloat(s[1])
            data = [EAV, DEAV]
        elif recType == 'G':
            E = convertFloat(line[9:19])
            DE = convertFloat(line[19:21])
            RI = convertFloat(line[21:29])
            DRI = convertFloat(line[29:31])
            M = convertFloat(line[31:41])
            MR = convertFloat(line[41:49])
            DMR = convertFloat(line[49:55])
            CC = convertFloat(line[55:62])
            DCC = convertFloat(line[62:64])
            TI = convertFloat(line[64:74])
            DTI = convertFloat(line[74:76])
            COIN = line[77]
            Q = line[79]
            data = [E, DE, RI, DRI, M, MR, DMR, CC, DCC, TI, DTI, COIN, Q]
            
    return data
    