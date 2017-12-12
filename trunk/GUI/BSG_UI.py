# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:39 2014

@author: leendert
"""
import sys
from shell import Shell, CommandError

from PySide import QtGui

from MainWindowGUI import Ui_MainWindow

import ConfigParser

import numpy as np

import os

import Utility as ut
        
def setCheckBoxState(cb, b):
    if (cb.isChecked() and not b) or (not cb.isChecked() and b):
        cb.toggle()

class BSG_UI(QtGui.QMainWindow):
    
    def __init__(self):
        QtGui.QMainWindow.__init__(self, None)
        self.ui = Ui_MainWindow()
        
        self.ui.setupUi(self)
        
        self.iniName = ''
        self.configName = ''
        self.execPath = ''
        self.exchangePath = ''
        
        self.spectrumCheckBoxes =  {self.ui.cb_phasespace: "Phasespace", self.ui.cb_fermi: "Fermi", self.ui.cb_radiative: "Radiative",\
        self.ui.cb_efs: "ESFiniteSize", self.ui.cb_ens: "ESShape", self.ui.cb_efs_deformation: "ESDeformation",\
        self.ui.cb_C: "C", self.ui.cb_isovector: "Isovector", self.ui.cb_coupled: "Connect", self.ui.cb_C_deformation: "CDeformation",\
        self.ui.cb_relativistic: "Relativistic", self.ui.cb_recoil: "Recoil", self.ui.cb_Q: "CoulombRecoil", \
        self.ui.cb_screening: "Screening", self.ui.cb_exchange: "Exchange", self.ui.cb_mismatch: "AtomicMismatch"}
        self.spectrumNME = {self.ui.dsb_bAc: "WeakMagnetism", self.ui.dsb_dAc: "InducedTensor",\
        self.ui.dsb_lambda: "Lambda"}
        self.spectrumDSB = {self.ui.dsb_begin: "Begin", self.ui.dsb_end: "End",\
        self.ui.dsb_step: "Step"}
        self.spectrumComboBoxes = {self.ui.cb_shape: "Shape"}
        self.constants = {self.ui.dsb_gA: "gA", self.ui.dsb_gAeff: "gAeff", self.ui.dsb_gP: "gP",\
        self.ui.dsb_gM: "gM"}
        self.computational = {self.ui.cb_nmeMethod: "Method", self.ui.dsb_energyMargin: "EnergyMargin",\
        self.ui.cb_forceSpin: "ForceSpin", self.ui.cb_reverseGhallagher: "ReversedGhallagher",\
        self.ui.cb_overrideCoupling: "OverrideSPCoupling", self.ui.cb_potential: "Potential",\
        self.ui.dsb_V0n: "Vneutron", self.ui.dsb_V0p: "Vproton", self.ui.dsb_V0Sn: "V0Sneutron",\
        self.ui.dsb_V0Sp: "V0Sproton", self.ui.dsb_Xn: "Xneutron", self.ui.dsb_Xp: "Xproton"}
        
        
        self.ui.gv_plotSpectrum.setLabels(bottom='Energy [keV]', left='Amplitude')
        self.ui.gv_plotSpectrum.showGrid(True, True, 0.5)
        
        self.sh = Shell()
        
        self.ui.b_enable_all.clicked.connect(self.enableAll)
        self.ui.b_disable_all.clicked.connect(self.disableAll)
        
        self.ui.cb_shape.addItems(("Fermi", "Mod. Gauss."))
        self.ui.cb_C_shape.addItems(("UCS", "Mod. Gauss."))
        
        self.ui.cb_process.addItems(("B-", "B+", "EC"))
        self.ui.cb_type.addItems(("Fermi", "Gamow-Teller", "Mixed"))
        self.ui.cb_PiM.addItems(("+", "-"))
        self.ui.cb_PiD.addItems(("+", "-"))
        self.ui.dsb_MixingRatio.setEnabled(False)
        self.ui.cb_type.currentIndexChanged[int].connect(self.enableMixingRatio)
        
        self.ui.cb_nmeMethod.addItem("ESP")
        self.ui.cb_potential.addItems(("SHO", "WS", "DWS"))
        
        self.ui.a_new.triggered.connect(self.newDecay)
        
        self.ui.a_loadIni.triggered.connect(self.loadIniFile)
        self.ui.a_loadConfig.triggered.connect(self.loadConfigFile)
        self.ui.a_saveConfig.triggered.connect(self.writeConfigFile)
        self.ui.a_saveIni.triggered.connect(self.writeIniFile)
        self.ui.a_bsgExec.triggered.connect(self.changeBSGExec)
        self.ui.a_exchangeData.triggered.connect(self.changeBSGExchange)
        
        self.ui.a_exit.triggered.connect(self.close)
        
        self.ui.a_clearPlots.triggered.connect(self.clearPlots)
        
        self.deformationFile = None
        self.radiiFile = None
        self.ensdfFolder = None
        self.ui.a_loadDeformation.triggered.connect(self.loadDeformation)
        self.ui.a_loadRadii.triggered.connect(self.loadRadii)
        self.ui.a_findTransition.triggered.connect(self.loadESNDFBranches)
        
        self.ui.a_about.triggered.connect(self.about)
        self.ui.a_feedback.triggered.connect(self.submitFeedback)
        
        self.ui.b_runBSG.clicked.connect(self.runBSG)
        
        self.plotColors = ('r','g','b','w','y')
        self.currentPlotIndex = 0
        
        self.findDefaults()
        
        self.log('Initialized...')
        
        self.statusBar().showMessage("Ready!")
        
    def log(self, message):
        self.ui.txt_lastActions.insertPlainText(message + '\n')
        self.ui.txt_lastActions.moveCursor(QtGui.QTextCursor.End)
        
    def status(self, message):
        self.statusBar().showMessage(message)
        
    def about(self):
        QtGui.QMessageBox.information(self, "About", "Beta spectrum Generator GUI v 0.1")
        
    def submitFeedback(self):
        QtGui.QMessageBox.information(self, "Send feedback", "Not yet implemented. Send emails to leendert.hayen@kuleuven.be")
        
    def findDefaults(self):
        self.log("Looking for defaults in current folder...")
        files = [f for f in os.listdir('.') if os.path.isfile(f)]
        for f in files:
            fullPath = os.path.join(os.getcwd(), f)
            if f == 'BSG':
                self.log("Found BSG Executable")
                self.execPath = fullPath
            elif f == 'config.txt':
                self.log("Found config.txt")
                self.loadConfigFile(fullPath)
            elif f == 'ExchangeData.dat':
                self.log("Found ExchangeData.dat")
                self.exchangePath = fullPath
                
    def newDecay(self):
        isotope, ok = QtGui.QInputDialog.getText(self, 'Isotope', 'Enter the isotope name (e.g. 238U)')
        if not ok:
            return
        import re
        A = int(re.findall(r'\d+', isotope)[0])
        name = isotope.replace(str(A), '').capitalize()
        Z = ut.atoms.index(name)+1
        
        self.ui.sb_AM.setValue(A)
        self.ui.sb_ZM.setValue(Z)
        
        button = QtGui.QMessageBox.question(self, 'ENSDF search', 'Search for transitions in ENSDF database?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        
        if button == QtGui.QMessageBox.Yes:
            self.loadESNDFBranches(Z, A)
        else:
            process = QtGui.QInputDialog.getItem(self, "Choose beta process", "Beta process", ('B-', 'B+', 'EC'), editable=False)[0]
            self.ui.cb_process.setCurrentIndex(self.ui.cb_process.findText(process))
            self.ui.sb_AD.setValue(A)
            self.ui.sb_ZD.setValue(Z+1 if process == 'B-' else Z-1)

        button = QtGui.QMessageBox.question(self, 'Deformation', 'Enter deformation info from the Moeller database?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        
        if button == QtGui.QMessageBox.Yes:
            self.loadDeformation()
            
        button = QtGui.QMessageBox.question(self, "Charge radii", 'Load charge radii from database?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
        
        if button == QtGui.QMessageBox.Yes:
            self.loadRadii()
            
                
    def clearPlots(self):
        self.ui.gv_plotSpectrum.clear()
        
    def loadESNDFBranches(self, Z, A):
        print('Entered ensdfBranches')
        if not self.ensdfFolder:
            self.ensdfFolder = QtGui.QFileDialog.getExistingDirectory(self, "Choose the ENSDF directory", ".")
            if self.ensdfFolder == '':
                self.ensdfFolder = None
                return
        fullPath = self.ensdfFolder + '/ensdf.{:03d}'.format(A)
        
        bc = ut.buildBetaBranches(fullPath, Z, A)
        items = []
        for i in bc.branches:
            s = '{}'
            items.append(s)
        
        
    def loadDeformation(self):
        if not self.deformationFile:
            self.deformationFile = QtGui.QFileDialog.getOpenFileName(self, "Choose the deformation file")[0]
        if self.deformationFile == '':
            return
        mDef = ut.loadDeformation(self.ui.sb_ZM.value(), self.ui.sb_AM.value(), self.deformationFile)
        dDef = ut.loadDeformation(self.ui.sb_ZD.value(), self.ui.sb_AD.value(), self.deformationFile)
        
        self.ui.dsb_Beta2M.setValue(mDef[0])
        self.ui.dsb_Beta4M.setValue(mDef[1])
        self.ui.dsb_Beta6M.setValue(mDef[2])
        self.ui.dsb_Beta2D.setValue(dDef[0])
        self.ui.dsb_Beta4D.setValue(dDef[1])
        self.ui.dsb_Beta6D.setValue(dDef[2])
        
    def loadRadii(self):
        if not self.radiiFile:
            self.radiiFile = QtGui.QFileDialog.getOpenFileName(self, "Choose the radii file")[0]
        if self.radiiFile == '':
            return
        mRad = ut.loadChargeRadius(self.ui.sb_ZM.value(), self.ui.sb_AM.value(), self.radiiFile)
        dRad = ut.loadChargeRadius(self.ui.sb_ZD.value(), self.ui.sb_AD.value(), self.radiiFile)
        
        self.ui.dsb_RM.setValue(mRad)
        self.ui.dsb_RD.setValue(dRad)
        
    def runBSG(self):
        self.status("Calculating...")
        outputName = self.ui.le_outputName.text()
        if self.execPath == '':
            QtGui.QErrorMessage(self).showMessage("Set the path for the generator executable.")
        if self.iniName == '':
            QtGui.QErrorMessage(self).showMessage("Choose an ini file.")
            return
        if self.configName == '':
            QtGui.QErrorMessage(self).showMessage("Choose a config file.")
            return
        if self.exchangePath == '':
            QtGui.QErrorMessage(self).showMessage("Choose the Exchange Data file.")
            return
        command = "{0} -i {1} -o {2} -c {3} -e {4}".format(self.execPath, self.iniName, outputName, self.configName, self.exchangePath)
        for key in self.spectrumCheckBoxes:
            command += " --Spectrum.{0}={1}".format(self.spectrumCheckBoxes[key], key.isChecked())
        for key in self.spectrumComboBoxes:
            command += " --Spectrum.{0}={1}".format(self.spectrumComboBoxes[key], key.currentText())
        for key in self.spectrumDSB:
            command += " --Spectrum.{0}={1}".format(self.spectrumDSB[key], key.value())
        if self.ui.cb_enforceNME.isChecked():
            for key in self.spectrumNME:
                command += " --Spectrum.{0}={1}".format(self.spectrumNME[key], key.value())
        print(command)
        try:
            self.sh.run(command)
            if self.sh.output(raw=True) != '':
                QtGui.QErrorMessage(self).showMessage(self.sh.output(raw=True))
                return
            spectrum = np.genfromtxt(outputName + '.raw')
            self.ui.gv_plotSpectrum.plot(x=spectrum[:, 1], y=spectrum[:,2], pen=self.plotColors[self.currentPlotIndex%len(self.plotColors)])
            self.currentPlotIndex += 1
            self.log("Spectrum calculation OK")
        except CommandError as e:
            QtGui.QErrorMessage(self).showMessage(('CommandError({0}): {1}'.format(e.errno, e.strerror)))
        self.status("Ready!")
        
    def changeBSGExec(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose BSG exec")[0]
        if filename == '':
            return
        self.execPath = filename
        
    def changeBSGExchange(self):
        filename = QtGui.QFileDialog.getOpenFileName(self, "Choose Exchange data file")[0]
        if filename == '':
            return
        self.exchangePath = filename
        
    def enableMixingRatio(self):
        if self.ui.cb_type.currentText() == "Mixed":
            self.ui.dsb_MixingRatio.setEnabled(True)
        else:
            self.ui.dsb_MixingRatio.setEnabled(False)
            
    def enableAll(self):
        for key in self.spectrumCheckBoxes:
            if not key.isChecked():
                key.toggle()
                
    def disableAll(self):
        for key in self.spectrumCheckBoxes:
            if key.isChecked():
                key.toggle()
        
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
        
        ut.writeIniFile(Zm, Zd, Am, Q, process, decayType, beta2m, beta4m, beta6m, beta2d, beta4d, beta6d, mRad, dRad, mJpi, dJpi, name=filename)
        self.log('Ini file written to %s.' % filename)
    def loadIniFile(self, filename = None):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self, "Choose .ini file")[0]
        if filename == '':
            return
        try:
            config = ConfigParser.RawConfigParser(allow_no_value=True)
            config.read(filename)
            
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
            self.ui.l_iniFilename.setText(filename.split('/')[-1])
            self.ui.l_iniFilename.setToolTip(filename)
            
            self.iniName = filename
            self.log("Loaded Ini file: %s." % filename)
        except:
            QtGui.QErrorMessage(self).showMessage("Failed to load %s file" % filename)
    def writeConfigFile(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save config file")[0]
        if filename == '':
            return
        try:
            ut.writeConfigFile(filename, self.ui.b_chooseConfigFolder.text(), self.computational,\
            self.constants, self.spectrumCheckBoxes, self.spectrumDSB, self.spectrumNME,\
            self.ui.cb_enforceNME.isChecked(), self.spectrumComboBoxes)
            self.log("Config file written to %s." % filename)
        except:
            QtGui.QErrorMessage(self).showMessage("Writing ini file failed.")
        
    def loadConfigFile(self, filename = None):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self, "Choose config file")[0]
        if filename == '':
            return
        #try:
        config = ConfigParser.RawConfigParser()
        config.read(filename)
        
        self.ui.b_chooseConfigFolder.setText(config.get('General', 'Folder'))
        for key in self.computational:
            try:
                key.setValue(config.getfloat('Computational', self.computational[key]))
            except:
                try:
                    setCheckBoxState(key, config.getboolean('Computational', self.computational[key]))
                except:
                    try:
                        key.setCurrentIndex(key.findText(config.get('Computational', self.computational[key])))
                    except:
                        pass
        
        for key in self.spectrumCheckBoxes:
            try:
                setCheckBoxState(key, config.getboolean('Spectrum', self.spectrumCheckBoxes[key]))
            except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                break
        for key in self.spectrumDSB:
            try:
                key.setValue(config.getfloat('Spectrum', self.spectrumDSB[key]))
            except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                break
        for key in self.spectrumNME:
            try:
                key.setValue(config.getfloat('Spectrum', self.spectrumNME[key]))
                setCheckBoxState(self.ui.cb_enforceNME, True)
            except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                break
        for key in self.spectrumComboBoxes:
            try:
                key.setCurrentIndex(key.findText(config.get('Spectrum', self.spectrumComboBoxes[key])))
            except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                break
            
        self.configName = filename
        self.log("Loaded config file: %s." % filename)
        #except:
            #QtGui.QErrorMessage(self).showMessage("Failed to load %s file" % filename)
        
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    mw = BSG_UI()
    mw.show()
    app.exec_()
    
    sys.exit()