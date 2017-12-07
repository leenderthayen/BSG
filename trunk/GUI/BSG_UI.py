# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:39 2014

@author: leendert
"""
import sys

from PySide import QtCore, QtGui, QtUiTools

from MainWindowGUI import Ui_MainWindow

import numpy as np

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
        radius = UF.getEltonNuclearRadius(A)
    return radius
    
def writeIniFile(Zm, Zd, A, Q, process, decayType, beta2m, beta4m, beta6m, beta2d, beta4d, \
beta6d, mRad, dRad, mJpi, dJpi, dE=0.0, name=None, prefix='', **kwargs):
    if not name:
        name = '%sZ%d_A%d_Q%.0f.ini' % (prefix, Z, A, Q)
    with open(name, 'w') as f:
        text = '[Transition]\nProcess=%s\nType=%s\nMixingRatio=0.0\nQValue=%f\n\n' \
        % (process, decayType, Q)
        f.write(text)
        text = '[Mother]\nZ=%d\nA=%d\nRadius=%f\nBeta2=%f\nBeta4=%f\nBeta6=%f\nSpinParity=%d\n\n' \
        % (Zm, A, mRad, beta2m, beta4m, beta6m, mJpi)
        f.write(text)
        text = '[Daughter]\nZ=%d\nA=%d\nRadius=%f\nBeta2=%f\nBeta4=%f\nBeta6=%f\nSpinParity=%d\nExcitationEnergy=%f\n' \
        % (Zd, A, dRad, beta2d, beta4d, beta6d, dJpi, dE)
        f.write(text)
    return name

class BSG_UI(QtGui.QMainWindow):
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self, None)
        self.ui = Ui_MainWindow()
        
        self.ui.setupUi(self)
        
        self.ui.cb_process.addItems(("B-", "B+", "EC"))
        self.ui.cb_type.addItems(("Fermi", "Gamow-Teller", "Mixed"))
        self.ui.dsb_MixingRatio.setEnabled(False)
        self.ui.cb_type.currentIndexChanged[int].connect(self.enableMixingRatio)
        
        self.ui.b_save_ini.clicked.connect(self.writeIniFile)
        
    def enableMixingRatio(self):
        if self.ui.cb_type.currentText() == "Mixed":
            self.ui.dsb_MixingRatio.setEnabled(True)
        else:
            self.ui.dsb_MixingRatio.setEnabled(False)
        
    def saveProfile(self):
        pass
        
    def writeIniFile(self):
        Zm = self.ui.sb_ZM.value()
        Zd = self.ui.sb_ZD.value()
        Am = self.ui.sb_AM.value()
        Ad = self.ui.sb_AD.value()
        
        process = self.ui.cb_process.currentText()

        if Am != Ad:
            QtGui.QErrorMessage(self).showMessage("Mass numbers do not agree.")
            return
        if Zd != (Zm + 1 if process == "B-" else Zm -1):
            QtGui.QErrorMessage(self).showMessage("Z values do not agree with process")
            return
            
        Q = self.ui.dsb_Q.value()
        decayType = self.ui.cb_type.currentText()
        mJ = self.ui.dsb_JM.value()
        dJ = self.ui.dsb_JM.value()
        
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save file", "", "")[0]
    def writeConfigFile(self):
        pass
        
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    mw = BSG_UI()
    mw.show()
    app.exec_()
    
    sys.exit()