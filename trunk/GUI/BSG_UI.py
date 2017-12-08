# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:39 2014

@author: leendert
"""
import sys
from shell import Shell, CommandError

from PySide import QtCore, QtGui, QtUiTools

from MainWindowGUI import Ui_MainWindow

import ConfigParser

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
beta6d, mRad, dRad, mJpi, dJpi, phl=None, dE=0.0, mE=0.0, name=None, prefix='', **kwargs):
    
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
    config.set('Mother', 'Z', Zm)
    config.set('Mother', 'A', A)
    config.set('Mother', 'Radius', mRad)
    config.set('Mother', 'beta2', beta2m)
    config.set('Mother', 'beta4', beta4m)
    config.set('Mother', 'beta6', beta6m)
    config.set('Mother', 'SpinParity', mJpi)
    config.set('Mother', 'ExcitationEnergy', mE)
    config.set('Daughter', 'Z', Zd)
    config.set('Daughter', 'A', A)
    config.set('Daughter', 'Radius', dRad)
    config.set('Daughter', 'beta2', beta2d)
    config.set('Daughter', 'beta4', beta4d)
    config.set('Daughter', 'beta6', beta6d)
    config.set('Daughter', 'SpinParity', dJpi)
    config.set('Daughter', 'ExcitationEnergy', dE)
    with open(name, 'wb') as configFile:
        config.write(configFile)
    
def writeConfigFile(name, directory, prefix, method, margin, forceSpin, reverseGhallagher,\
overrideCoupling, potential, V0n, V0p, Xn, Xp, V0Sn, V0Sp, gA, gAeff, gP, gM):
    config = ConfigParser.RawConfigParser(allow_no_value=True)
    config.optionxform = str
    config.add_section('General')
    config.add_section('Computational')
    config.add_section('Constants')
    config.set('General', 'Folder', directory)
    config.set('General', 'Prefix', prefix)
    config.set('Computational', 'Method', method)
    config.set('Computational', 'EnergyMargin', margin)
    config.set('Computational', 'ForceSpin', forceSpin)
    config.set('Computational', 'ReversedGhallagher', reverseGhallagher)
    config.set('Computational', 'OverrideSPCoupling', overrideCoupling)
    config.set('Computational', 'Potential', potential)
    config.set('Computational', 'Vneutron', V0n)
    config.set('Computational', 'Vproton', V0p)
    config.set('Computational', 'Xneutron', Xn)
    config.set('Computational', 'Xproton', Xp)
    config.set('Computational', 'V0Sneutron', V0Sn)
    config.set('Computational', 'V0Sproton', V0Sp)
    config.set('Constants', 'gA', gA)
    config.set('Constants', 'gAeff', gAeff)
    config.set('Constants', 'gP', gP)
    config.set('Constants', 'gM', gM)
    with open(name, 'wb') as configFile:
        config.write(configFile)

class BSG_UI(QtGui.QMainWindow):
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self, None)
        self.ui = Ui_MainWindow()
        
        self.ui.setupUi(self)
        
        self.sh = Shell()
        
        self.checkBoxes =  {self.ui.cb_phasespace: "phase", self.ui.cb_fermi: "fermi", self.ui.cb_radiative: "radiative",\
        self.ui.cb_efs: "l0", self.ui.cb_ens: "U", self.ui.cb_efs_deformation: "deformation",\
        self.ui.cb_C: "C", self.ui.cb_isovector: "isovector", self.ui.cb_coupled: "connect", self.ui.cb_C_deformation: "deformation",\
        self.ui.cb_relativistic: "relativistic", self.ui.cb_recoil: "recoil", self.ui.cb_Q: "Q", \
        self.ui.cb_screening: "screening", self.ui.cb_exchange: "exchange"}
        
        self.ui.b_enable_all.clicked.connect(self.enableAll)
        self.ui.b_disable_all.clicked.connect(self.disableAll)
        
        self.ui.cb_shape.addItems(("Fermi", "Mod. Gauss."))
        self.ui.cb_C_shape.addItems(("UCS", "Mod. Gauss."))
        #self.ui.rb_dW0.toggled.connect(self.enableDeltaW0)
        #self.ui.rb_overlap.toggled.connect(self.enableDeltaW0)
        
        self.ui.cb_process.addItems(("B-", "B+", "EC"))
        self.ui.cb_type.addItems(("Fermi", "Gamow-Teller", "Mixed"))
        self.ui.cb_PiM.addItems(("+", "-"))
        self.ui.cb_PiD.addItems(("+", "-"))
        self.ui.dsb_MixingRatio.setEnabled(False)
        self.ui.cb_type.currentIndexChanged[int].connect(self.enableMixingRatio)
        
        self.ui.cb_nmeMethod.addItem("ESP")
        self.ui.cb_potential.addItems(("SHO", "WS", "DWS"))
        
        self.ui.b_saveConfigFile.clicked.connect(self.writeConfigFile)
        self.ui.b_loadConfigFile.clicked.connect(self.loadConfigFile)
        
        self.ui.b_save_ini.clicked.connect(self.writeIniFile)
        self.ui.b_load_ini.clicked.connect(self.loadIniFile)
        self.ui.b_runBSG.clicked.connect(self.runBSG)
        self.ui.b_changeIni.clicked.connect(self.changeBSGIniFile)
        self.ui.b_changeConfig.clicked.connect(self.changeBSGConfigFile)
        self.ui.b_changeProfile.clicked.connect(self.changeBSGProfileFile)
        self.ui.b_bsgExchange.clicked.connect(self.changeBSGExchange)
        
        self.log('Initialized...')
        
    def log(self, message):
        self.ui.txt_lastActions.insertPlainText(message + '\n')
        self.ui.txt_lastActions.moveCursor(QtGui.QTextCursor.End)
        
    def runBSG(self):
        execName = self.ui.l_bsgExec.text()
        iniName = self.ui.l_bsgIni.text()
        configName = self.ui.l_bsgConfig.text()
        outputName = self.ui.le_outputName.text()
        exchangeName = self.ui.l_bsgExchange.text()
        if iniName == '(blank)':
            QtGui.QErrorMessage(self).showMessage("Choose an ini file.")
            return
        if configName == '(blank)':
            QtGui.QErrorMessage(self).showMessage("Choose a config file.")
            return
        if exchangeName == '(blank)':
            QtGui.QErrorMessage(self).showMessage("Choose the Exchange Data file.")
            return
        command = "./{0} -i {1} -o {2} -c {3} -e {4}".format(execName, iniName, outputName, configName, exchangeName)
        print(command)
        try:
            self.sh.run(command)
            if self.sh.output(raw=True) != '':
                QtGui.QErrorMessage(self).showMessage(self.sh.output(raw=True))
                return
            spectrum = np.genfromtxt(outputName + '.raw')
            self.ui.gv_plotSpectrum.plot(x=spectrum[:, 1], y=spectrum[:,2])
            self.log("Spectrum calculation OK")
        except CommandError as e:
            QtGui.QErrorMessage(self).showMessage(('CommandError({0}): {1}'.format(e.errno, e.strerror)))
        
    def changeBSGExec(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose BSG exec")[0]
        if filename == '':
            return
        self.ui.l_bsgExec.setText(filename)
        self.ui.l_bsgExec.setToolTip(filename)
        
    def changeBSGIniFile(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose .ini file")[0]
        if filename == '':
            return
        self.loadIniFile(filename)
        self.ui.l_bsgIni.setText(filename)
        self.ui.l_bsgIni.setToolTip(filename)
        
    def changeBSGConfigFile(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose config file")[0]
        if filename == '':
            return
        self.loadConfigFile(filename)
        self.ui.l_bsgConfig.setText(filename)
        self.ui.l_bsgConfig.setToolTip(filename)
        
    def changeBSGProfileFile(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose profile file")[0]
        if filename == '':
            return
        self.loadProfileFile(filename)
        self.ui.l_bsgProfile.setText(filename)
        self.ui.l_bsgProfile.setToolTip(filename)
        
    def changeBSGExchange(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose Exchange data file")[0]
        if filename == '':
            return
        self.ui.l_bsgExchange.setText(filename)
        self.ui.l_bsgExchange.setToolTip(filename)
        
    def enableMixingRatio(self):
        if self.ui.cb_type.currentText() == "Mixed":
            self.ui.dsb_MixingRatio.setEnabled(True)
        else:
            self.ui.dsb_MixingRatio.setEnabled(False)
#    def enableDeltaW0(self):
#        if self.ui.rb_dW0.isChecked():
#            self.ui.dsb_dW0.setEnabled(True)
#        else:
#            self.ui.dsb_dW0.setEnabled(False)
            
    def enableAll(self):
        for key in self.checkBoxes:
            if not key.isChecked():
                key.toggle()
    def disableAll(self):
        for key in self.checkBoxes:
            if key.isChecked():
                key.toggle()
        
    def saveProfile(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save profile")[0]
        if filename == '':
            return
        
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
        
        mRad = self.ui.dsb_RM.value()
        dRad = self.ui.dsb_RD.value()
        
        beta2m = self.ui.dsb_Beta2M.value()
        beta4m = self.ui.dsb_Beta4M.value()
        beta6m = self.ui.dsb_Beta6M.value()
        beta2d = self.ui.dsb_Beta2D.value()
        beta4d = self.ui.dsb_Beta4D.value()
        beta6d = self.ui.dsb_Beta6D.value()
        
        mJpi = int(mJ*2)*(1 if self.ui.cb_PiM.currentText() == "+" else -1)
        dJpi = int(dJ*2)*(1 if self.ui.cb_PiD.currentText() == "+" else -1)
        
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save .ini file")[0]
        if filename == '':
            return
        
        writeIniFile(Zm, Zd, Am, Q, process, decayType, beta2m, beta4m, beta6m, beta2d, beta4d, beta6d, mRad, dRad, mJpi, dJpi, name=filename)
        self.log('Ini file written to %s.' % filename)
    def loadIniFile(self, filename = None):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self, "Choose .ini file")[0]
        if filename == '':
            return
        try:
            config = ConfigParser.RawConfigParser(allow_no_value=True)
            config.read(filename)
            
            self.ui.l_iniFilename.setText(filename.split('/')[-1])
            self.ui.cb_process.setCurrentIndex(self.ui.cb_process.findText(config.get('Transition', 'Process')))
            self.ui.cb_type.setCurrentIndex(self.ui.cb_type.findText(config.get('Transition', 'Type')))
            self.ui.dsb_MixingRatio.setValue(config.getfloat('Transition', 'MixingRatio'))
            self.ui.dsb_Q.setValue(config.getfloat('Transition', 'QValue'))
            self.ui.sb_ZM.setValue(config.getint('Mother', 'Z'))
            self.ui.sb_AM.setValue(config.getint('Mother', 'A'))
            self.ui.dsb_RM.setValue(config.getfloat('Mother', 'Radius'))
            self.ui.dsb_JM.setValue(abs(config.getint('Mother', 'SpinParity'))/2.)
            self.ui.cb_PiM.setCurrentIndex(0 if config.getint('Mother', 'SpinParity') > 0 else 1)
            self.ui.dsb_Beta2M.setValue(config.getfloat('Mother', 'Beta2'))
            self.ui.dsb_Beta4M.setValue(config.getfloat('Mother', 'Beta6'))
            self.ui.dsb_Beta6M.setValue(config.getfloat('Mother', 'Beta6'))
            self.ui.sb_ZD.setValue(config.getint('Daughter', 'Z'))
            self.ui.sb_AD.setValue(config.getint('Daughter', 'A'))
            self.ui.dsb_RD.setValue(config.getfloat('Daughter', 'Radius'))
            self.ui.dsb_JD.setValue(abs(config.getint('Daughter', 'SpinParity'))/2.)
            self.ui.cb_PiD.setCurrentIndex(0 if config.getint('Daughter', 'SpinParity') > 0 else 1)
            self.ui.dsb_Beta2D.setValue(config.getfloat('Daughter', 'Beta2'))
            self.ui.dsb_Beta4D.setValue(config.getfloat('Daughter', 'Beta6'))
            self.ui.dsb_Beta6D.setValue(config.getfloat('Daughter', 'Beta6'))
            self.log("Loaded Ini file: %s." % filename)
        except:
            QtGui.QErrorMessage(self).showMessage("Failed to load %s file" % filename)
    def writeConfigFile(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save config file")[0]
        if filename == '':
            return
        
        writeConfigFile(filename, self.ui.b_chooseConfigFolder.text(), \
        self.ui.le_configPrefix.text(), self.ui.cb_nmeMethod.currentText(),\
        self.ui.dsb_energyMargin.value(), self.ui.cb_forceSpin.isChecked(), self.ui.cb_reverseGhallagher.isChecked(),\
        self.ui.cb_overrideCoupling.isChecked(), self.ui.cb_potential.currentText(),\
        self.ui.dsb_V0n.value(), self.ui.dsb_V0p.value(), self.ui.dsb_Xn.value(),\
        self.ui.dsb_Xp.value(), self.ui.dsb_V0Sn.value(), self.ui.dsb_V0Sp.value(),\
        self.ui.dsb_gA.value(), self.ui.dsb_gAeff.value(), self.ui.dsb_gP.value(),\
        self.ui.dsb_gM.value())
        self.log("Config file written to %s." % filename)
        
    def loadConfigFile(self, filename = None):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self, "Choose config file")[0]
        if filename == '':
            return
        try:
            config = ConfigParser.RawConfigParser()
            config.read(filename)
            
            self.ui.b_chooseConfigFolder.setText(config.get('General', 'Folder'))
            self.ui.le_configPrefix.setText(config.get('General', 'Prefix'))
            self.ui.cb_nmeMethod.setCurrentIndex(self.ui.cb_nmeMethod.findText(config.get('Computational', 'Method')))
            self.ui.dsb_energyMargin.setValue(config.getfloat('Computational', 'EnergyMargin'))
            self.ui.cb_potential.setCurrentIndex(self.ui.cb_potential.findText(config.get('Computational', 'Potential')))
            self.ui.dsb_V0n.setValue(config.getfloat('Computational', 'Vneutron'))
            self.ui.dsb_V0p.setValue(config.getfloat('Computational', 'Vproton'))
            self.ui.dsb_Xn.setValue(config.getfloat('Computational', 'Xneutron'))
            self.ui.dsb_Xp.setValue(config.getfloat('Computational', 'Xproton'))
            self.ui.dsb_V0Sn.setValue(config.getfloat('Computational', 'V0Sneutron'))
            self.ui.dsb_V0Sp.setValue(config.getfloat('Computational', 'V0Sproton'))
            self.ui.dsb_gA.setValue(config.getfloat('Constants', 'gA'))
            self.ui.dsb_gAeff.setValue(config.getfloat('Constants', 'gAeff'))
            self.ui.dsb_gP.setValue(config.getfloat('Constants', 'gP'))
            self.ui.dsb_gM.setValue(config.getfloat('Constants', 'gM'))
            
            if (self.ui.cb_overrideCoupling.isChecked() and not config.getboolean('Computational', 'OverrideSPCoupling')) or (not self.ui.cb_overrideCoupling.isChecked() and config.getboolean('Computational', 'OverrideSPCoupling')):
                self.ui.cb_overrideCoupling.toggle()
            if (self.ui.cb_reverseGhallagher.isChecked() and not config.getboolean('Computational', 'ReversedGhallagher')) or (not self.ui.cb_reverseGhallagher.isChecked() and config.getboolean('Computational', 'ReversedGhallagher')):
                self.ui.cb_reverseGhallagher.toggle()
            if (self.ui.cb_forceSpin.isChecked() and not config.getboolean('Computational', 'ForceSpin')) or (not self.ui.cb_forceSpin.isChecked() and config.getboolean('Computational', 'ForceSpin')):
                self.ui.cb_forceSpin.toggle()
            self.log("Loaded config file: %s." % filename)
        except:
            QtGui.QErrorMessage(self).showMessage("Failed to load %s file" % filename)
        
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    mw = BSG_UI()
    mw.show()
    app.exec_()
    
    sys.exit()