#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:39 2014

@author: leendert
"""
import sys
from shell import Shell, CommandError

import qdarkstyle
import os

from PySide import QtGui

from MainWindowGUI import Ui_MainWindow

import ConfigParser

import numpy as np

import utils.utilities as ut

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
        self.robtdFile = ''

        self.unsavedTransitionChanges = list()

        self.spectrumCheckBoxes =  {self.ui.cb_phasespace: "Phasespace", self.ui.cb_fermi: "Fermi", self.ui.cb_radiative: "Radiative",\
        self.ui.cb_efs: "ESFiniteSize", self.ui.cb_ens: "ESShape", self.ui.cb_efs_deformation: "ESDeformation",\
        self.ui.cb_C: "C", self.ui.cb_isovector: "Isovector", self.ui.cb_coupled: "Connect", self.ui.cb_C_deformation: "CDeformation",\
        self.ui.cb_relativistic: "Relativistic", self.ui.cb_recoil: "Recoil", self.ui.cb_Q: "CoulombRecoil", \
        self.ui.cb_screening: "Screening", self.ui.cb_exchange: "Exchange", self.ui.cb_mismatch: "AtomicMismatch"}

        self.spectrumNME = {self.ui.dsb_bAc: "WeakMagnetism", self.ui.dsb_dAc: "InducedTensor",\
        self.ui.dsb_lambda: "Lambda"}

        self.spectrumDSB = {self.ui.dsb_begin: "Begin", self.ui.dsb_end: "End",\
        self.ui.dsb_stepSize: "StepSize", self.ui.sb_steps: "Steps"}

        self.spectrumComboBoxes = {self.ui.cb_shape: "Shape"}

        self.constants = {self.ui.dsb_gA: "gA", self.ui.dsb_gAeff: "gAeff", self.ui.dsb_gP: "gP",\
        self.ui.dsb_gM: "gM"}

        self.computationalComboBoxes = {self.ui.cb_nmeMethod: "Method", self.ui.cb_potential: "Potential"}

        self.computationalCheckBoxes = {self.ui.cb_forceSpin: "ForceSpin",
        self.ui.cb_reverseGhallagher: "ReversedGhallagher", self.ui.cb_overrideCoupling: "OverrideSPCoupling"}

        self.computationalDSB = {self.ui.dsb_energyMargin: "EnergyMargin",
        self.ui.dsb_V0n: "Vneutron", self.ui.dsb_V0p: "Vproton", self.ui.dsb_V0Sn: "V0Sneutron",
        self.ui.dsb_V0Sp: "V0Sproton", self.ui.dsb_Xn: "Xneutron", self.ui.dsb_Xp: "Xproton"}


        self.ui.gv_plotSpectrum.setLabels(bottom='Energy [keV]', left='Amplitude [a.u.]')
        self.ui.gv_plotSpectrum.showGrid(True, True, 0.5)

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

        self.ui.cb_nmeMethod.addItems(("ROBTD", "ESP"))
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

        self.ui.b_chooseROBTDFile.clicked.connect(self.setROBTDFile)
        self.ui.b_chooseConfigFile.clicked.connect(self.loadConfigFile)

        self.ui.dsb_bAc.valueChanged.connect(self.enforceMatrixElements)
        self.ui.dsb_dAc.valueChanged.connect(self.enforceMatrixElements)
        self.ui.dsb_lambda.valueChanged.connect(self.enforceMatrixElements)

        self.ui.cb_nmeMethod.currentIndexChanged[unicode].connect(self.setESPVisibility)

        self.ui.b_resetSSO.clicked.connect(self.resetSpectrumShapeOptions)
        self.ui.b_resetNSO.clicked.connect(self.resetNuclearStructureOptions)

        self.ui.rb_stepSize.toggled.connect(self.ui.dsb_stepSize.setEnabled)
        self.ui.rb_steps.toggled.connect(self.ui.sb_steps.setEnabled)

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

    def resetSpectrumShapeOptions(self):
        self.loadConfigFile(self.configName, nuclearStructure=False)
        self.log('Resetting spectrum shape options back to values from %s.' % self.configName)

    def resetNuclearStructureOptions(self):
        self.loadConfigFile(self.configName, spectrumShape=False)
        self.log('Resetting nuclear structure options back to values from %s.' % self.configName)

    def setESPVisibility(self, value):
        if value == 'ESP':
            self.ui.gb_ESP.setEnabled(True)
        else:
            self.ui.gb_ESP.setEnabled(False)

    def addUnsavedChanges(self, component):
        self.unsavedTransitionChanges.append(component)

    def clearUnsavedChanges(self):
        self.unsavedTransitionChanges = list()

    def askSaveChanges(self):
        if len(self.unsavedTransitionChanges):
            s = 'There are unsaved changes in the transition information:\n\n'
            for i in self.unsavedTransitionChanges:
                s += i + '\n'
            s += '\nThese will not be taken into account. Save info?'
            button = QtGui.QMessageBox.question(self, 'Unsaved changes',
            s, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if button == QtGui.QMessageBox.Yes:
                self.writeIniFile()

    def checkUnsavedTransitionChanges(self):
        self.clearUnsavedChanges()
        if self.iniName == '':
            self.addUnsavedChanges("No information currently saved.")
        else:
            config = ConfigParser.RawConfigParser(allow_no_value=True)
            config.read(self.iniName)

            if config.get('Transition', 'Process') != self.ui.cb_process.currentText():
                self.addUnsavedChanges('Process (Prev. %s)' % config.get('Transition', 'Process'))
            if config.get('Transition', 'Type') != self.ui.cb_type.currentText():
                self.addUnsavedChanges('Type (Prev. %s)' % config.get('Transition', 'Type'))
            if config.getfloat('Transition', 'MixingRatio') != self.ui.dsb_MixingRatio.value():
                self.addUnsavedChanges('Mixing Ratio (Prev. %s)' % config.get('Transition', 'MixingRatio'))
            if config.getfloat('Transition', 'QValue') != self.ui.dsb_Q.value():
                self.addUnsavedChanges('Q Value (Prev. %s)' % config.get('Transition', 'QValue'))
            try:
                if config.getfloat('Transition', 'PartialHalflife') != self.ui.dsb_halflife.value():
                    self.addUnsavedChanges('Partial Halflife (Prev. %s)' % config.get('Transition', 'PartialHalflife'))
            except ConfigParser.NoOptionError:
                pass
            try:
                if config.getfloat('Transition', 'LogFt') != self.ui.dsb_logft.value():
                    self.addUnsavedChanges('Log Ft (Prev. %s)' % config.get('Transition', 'LogFt'))
            except ConfigParser.NoOptionError:
                pass

            if config.getint('Mother', 'Z') != self.ui.sb_ZM.value():
                self.addUnsavedChanges('Mother Z (Prev. %d)' % config.getint('Mother', 'Z'))
            if config.getint('Mother', 'A') != self.ui.sb_AM.value():
                self.addUnsavedChanges('Mother A (Prev. %d)' % config.getint('Mother', 'A'))
            if config.getfloat('Mother', 'Radius') != self.ui.dsb_RM.value():
                self.addUnsavedChanges('Mother Radius (Prev. %.2f)' % config.getfloat('Mother', 'Radius'))
            if abs(config.getint('Mother', 'SpinParity')/2.) != self.ui.dsb_JM.value():
                self.addUnsavedChanges('Mother Spin (Prev. %.1f)' % abs(config.getint('Mother', 'SpinParity')/2.))
            if None:
                pass
            if config.getfloat('Mother', 'Beta2') != self.ui.dsb_Beta2M.value():
                self.addUnsavedChanges('Mother beta2 (Prev. %.3f)' % config.getfloat('Mother', 'Beta2'))
            if config.getfloat('Mother', 'Beta4') != self.ui.dsb_Beta4M.value():
                self.addUnsavedChanges('Mother beta4 (Prev. %.3f)' % config.getfloat('Mother', 'Beta4'))
            if config.getfloat('Mother', 'Beta6') != self.ui.dsb_Beta6M.value():
                self.addUnsavedChanges('Mother beta6 (Prev. %.3f)' % config.getfloat('Mother', 'Beta6'))

            if config.getint('Daughter', 'Z') != self.ui.sb_ZD.value():
                self.addUnsavedChanges('Daughter Z (Prev. %d)' % config.getint('Daughter', 'Z'))
            if config.getint('Daughter', 'A') != self.ui.sb_AD.value():
                self.addUnsavedChanges('Daughter A (Prev. %d)' % config.getint('Daughter', 'A'))
            if config.getfloat('Daughter', 'Radius') != self.ui.dsb_RD.value():
                self.addUnsavedChanges('Daughter Radius (Prev. %.2f)' % config.getfloat('Daughter', 'Radius'))
            if abs(config.getint('Daughter', 'SpinParity')/2.) != self.ui.dsb_JD.value():
                self.addUnsavedChanges('Daughter Spin (Prev. %.1f)' % abs(config.getint('Daughter', 'SpinParity')/2.))
            if None:
                pass
            if config.getfloat('Daughter', 'Beta2') != self.ui.dsb_Beta2D.value():
                self.addUnsavedChanges('Daughter beta2 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta2'))
            if config.getfloat('Daughter', 'Beta4') != self.ui.dsb_Beta4D.value():
                self.addUnsavedChanges('Daughter beta4 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta4'))
            if config.getfloat('Daughter', 'Beta6') != self.ui.dsb_Beta6D.value():
                self.addUnsavedChanges('Daughter beta6 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta6'))

        self.askSaveChanges()

    def about(self):
        QtGui.QMessageBox.information(self, "About", "Beta spectrum Generator GUI v 0.2")

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
            elif f == 'FRDM2012.dat':
                self.log('Found FRDM2012.dat')
                self.deformationFile = fullPath
            elif f == 'ChargeRadii.dat':
                self.log('Found ChargeRadii.dat')
                self.radiiFile = fullPath

    def enforceMatrixElements(self):
        if not self.ui.cb_enforceNME.isChecked():
            self.ui.cb_enforceNME.toggle()

    def newDecay(self):
        self.checkUnsavedTransitionChanges()

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

        button = QtGui.QMessageBox.question(self, "Save as", 'Save transition in ini file?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if button == QtGui.QMessageBox.Yes:
            self.writeIniFile()

    def clearPlots(self):
        self.ui.gv_plotSpectrum.clear()
        self.currentPlotIndex = 0

    def loadESNDFBranches(self, Z, A):
        print('Entered ensdfBranches')
        if not self.ensdfFolder:
            self.ensdfFolder = QtGui.QFileDialog.getExistingDirectory(self, "Choose the ENSDF directory", ".")
            if self.ensdfFolder == '':
                self.ensdfFolder = None
                return
        fullPath = self.ensdfFolder + '/ensdf.{:03d}'.format(A)

        bbs = ut.getENSDFBetaBranches(str(fullPath), Z, A, 0)
        items = []
        for i in bbs:
            s = '{}: {}{} ([{}] {}) -> {}{} ([{}] {}); Endpoint: {} keV, t1/2: {:.3f} s'.format(i.process, A, ut.atoms[Z-1], i.motherLevel.Jpi, i.motherLevel.E,\
            A, ut.atoms[Z if i.process == 'B-' else Z-2], i.daughterLevel.Jpi, i.daughterLevel.E, i.E, i.partialHalflife)

            items.append(s)

        choice  = QtGui.QInputDialog.getItem(self, "Transition", "Choose the beta transition:", items, editable=False)[0]
        branch = bbs[items.index(choice)]

        self.ui.sb_ZM.setValue(Z)
        self.ui.sb_AM.setValue(A)
        self.ui.sb_ZD.setValue(Z+1 if branch.process == 'B-' else Z-1)
        self.ui.sb_AD.setValue(A)
        self.ui.cb_process.setCurrentIndex(self.ui.cb_process.findText(branch.process))
        self.ui.cb_type.setCurrentIndex(1)
        self.ui.dsb_motherEn.setValue(branch.motherLevel.E)
        self.ui.dsb_daughterEn.setValue(branch.daughterLevel.E)
        self.ui.dsb_Q.setValue(branch.E-branch.motherLevel.E+branch.daughterLevel.E)
        mJpi = branch.motherLevel.Jpi
        self.ui.dsb_JM.setValue(float(mJpi.split('/2')[0])/2. if '/2' in mJpi else float(mJpi[:-1]))
        dJpi = branch.daughterLevel.Jpi
        self.ui.dsb_JD.setValue(float(dJpi.split('/2')[0])/2. if '/2' in dJpi else float(dJpi[:-1]))
        self.ui.cb_PiM.setCurrentIndex(self.ui.cb_PiM.findText(mJpi[-1]))
        self.ui.cb_PiD.setCurrentIndex(self.ui.cb_PiD.findText(dJpi[-1]))
        self.ui.dsb_halflife.setValue(branch.partialHalflife)
        self.ui.dsb_logft.setValue(branch.logft)

        self.setTransitionLabel()

    def loadDeformation(self):
        if not self.deformationFile:
            self.deformationFile = QtGui.QFileDialog.getOpenFileName(self, "Choose the deformation file")[0]
        if self.deformationFile == '':
            return
        mDef = ut.loadDeformation(self.ui.sb_ZM.value(), self.ui.sb_AM.value(), self.deformationFile)
        dDef = ut.loadDeformation(self.ui.sb_ZD.value(), self.ui.sb_AD.value(), self.deformationFile)

        self.ui.dsb_Beta2M.setValue(mDef[0])
        self.ui.dsb_Beta4M.setValue(mDef[2])
        self.ui.dsb_Beta6M.setValue(mDef[3])
        self.ui.dsb_Beta2D.setValue(dDef[0])
        self.ui.dsb_Beta4D.setValue(dDef[2])
        self.ui.dsb_Beta6D.setValue(dDef[3])

    def loadRadii(self):
        if not self.radiiFile:
            self.radiiFile = QtGui.QFileDialog.getOpenFileName(self, "Choose the radii file")[0]
        if self.radiiFile == '':
            return
        mRad = ut.loadChargeRadius(self.ui.sb_ZM.value(), self.ui.sb_AM.value(), self.radiiFile)
        dRad = ut.loadChargeRadius(self.ui.sb_ZD.value(), self.ui.sb_AD.value(), self.radiiFile)

        self.ui.dsb_RM.setValue(mRad)
        self.ui.dsb_RD.setValue(dRad)

    def setTransitionLabel(self):
        Am = self.ui.sb_AM.value()
        Ad = self.ui.sb_AD.value()
        Zm = self.ui.sb_ZM.value()
        Zd = self.ui.sb_ZD.value()
        JpiM = str(int(self.ui.dsb_JM.value())) if self.ui.dsb_JM.value() % 1 == 0. else str(int(2*self.ui.dsb_JM.value())) + '/2'
        JpiM += self.ui.cb_PiM.currentText()
        JpiD = str(int(self.ui.dsb_JD.value())) if self.ui.dsb_JD.value() % 1 == 0. else str(int(2*self.ui.dsb_JD.value())) + '/2'
        JpiD += self.ui.cb_PiD.currentText()
        s = '{}{} ([{}] {}) -> {}{} ([{}] {})'.format(Am, ut.atoms[Zm-1],
        JpiM, self.ui.dsb_motherEn.value(), Ad, ut.atoms[Zd-1],
        JpiD, self.ui.dsb_daughterEn.value())
        self.ui.l_transitionName.setText(s)
        self.ui.l_transitionName.setToolTip(s)
        self.ui.l_transitionName.setStatusTip(s)

    def runBSG(self):
        self.checkUnsavedTransitionChanges()

        outputName = self.ui.le_outputName.text()
        if self.execPath == '':
            QtGui.QErrorMessage(self).showMessage("Set the path for the generator executable.")
        if self.iniName == '':
            QtGui.QErrorMessage(self).showMessage("No ini file found.")
            return
        if self.ui.rb_configFile.isChecked() and self.configName == '':
            QtGui.QErrorMessage(self).showMessage("No configuration file specified.")
            return
        if self.ui.cb_exchange.isChecked() and self.exchangePath == '':
            QtGui.QErrorMessage(self).showMessage("Choose the Exchange Data file, or turn off atomic exchange.")
            return

        self.log("Performing BSG Calculation...")
        self.status("Calculating...")
        command = "{0} -i {1} -o {2} -c {3}".format(self.execPath, self.iniName, outputName, self.configName)
        if self.ui.cb_exchange.isChecked():
            command += " -e {}".format(self.exchangePath)

        if not self.ui.rb_configFile.isChecked():
            for key in self.spectrumCheckBoxes:
                command += " --Spectrum.{0}={1}".format(self.spectrumCheckBoxes[key], key.isChecked())
            for key in self.spectrumComboBoxes:
                command += " --Spectrum.{0}={1}".format(self.spectrumComboBoxes[key], key.currentText())
            for key in self.spectrumDSB:
                if key.isEnabled():
                    command += " --Spectrum.{0}={1}".format(self.spectrumDSB[key], key.value())
            if self.ui.cb_enforceNME.isChecked():
                for key in self.spectrumNME:
                    command += " --Spectrum.{0}={1}".format(self.spectrumNME[key], key.value())
            for key in self.computationalComboBoxes:
                command += ' --Computational.{0}={1}'.format(self.computationalDSB[key], key.currentText())
            for key in self.computationalCheckBoxes:
                command += ' -- Computational.{0}={1}'.format(self.computationalComboBoxes[key], key.isChecked())
            for key in self.computationalDSB:
                command += ' --Computational.{0}={1}'.format(self.computationalDSB[key], key.value())

        print("Executing command: %s" % command)
        try:
            sh = Shell()
            sh.run(command)
            if sh.output(raw=True) != '':
                QtGui.QErrorMessage(self).showMessage(sh.output(raw=True))
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

    def setROBTDFile(self, filename=None):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self, "Choose ROBTD data file")[0]
        if filename == '':
            return
        self.robtdFile = filename
        self.ui.b_chooseROBTDFile.setText(filename.split('/')[-1])
        self.ui.b_chooseROBTDFile.setToolTip(filename)
        self.ui.b_chooseROBTDFile.setStatusTip(filename)

    def setCurrentIniFile(self, filename):
        self.iniName = filename
        self.ui.l_iniFilename.setText(filename.split('/')[-1])
        self.ui.l_iniFilename.setToolTip(filename)
        self.ui.l_iniFilename.setStatusTip(filename)

    def setCurrentConfigFile(self, filename):
        self.configName = filename
        self.ui.b_chooseConfigFile.setText(filename.split('/')[-1])
        self.ui.b_chooseConfigFile.setToolTip(filename)
        self.ui.b_chooseConfigFile.setStatusTip(filename)

        self.ui.b_resetSSO.setText("Reset to %s" % filename.split('/')[-1])
        self.ui.b_resetNSO.setText("Reset to %s" % filename.split('/')[-1])

        if not self.ui.rb_configFile.isChecked():
            self.ui.rb_configFile.toggle()

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
        dJ = self.ui.dsb_JD.value()

        mRad = self.ui.dsb_RM.value()
        dRad = self.ui.dsb_RD.value()

        beta2m = self.ui.dsb_Beta2M.value()
        beta4m = self.ui.dsb_Beta4M.value()
        beta6m = self.ui.dsb_Beta6M.value()
        beta2d = self.ui.dsb_Beta2D.value()
        beta4d = self.ui.dsb_Beta4D.value()
        beta6d = self.ui.dsb_Beta6D.value()

        mE = self.ui.dsb_motherEn.value()
        dE = self.ui.dsb_daughterEn.value()

        mJpi = int(mJ*2)*(1 if self.ui.cb_PiM.currentText() == "+" else -1)
        dJpi = int(dJ*2)*(1 if self.ui.cb_PiD.currentText() == "+" else -1)

        phl = self.ui.dsb_halflife.value()
        phl = phl if phl > 0 else None

        logft = self.ui.dsb_logft.value()
        logft = logft if logft > 0 else None

        robtdFile = self.robtdFile if not self.robtdFile == '' else None

        filename = QtGui.QFileDialog.getSaveFileName(self, "Save .ini file")[0]
        if filename == '':
            return

        ut.writeIniFile(Zm, Zd, Am, Q, process, decayType, beta2m, beta4m, beta6m,
        beta2d, beta4d, beta6d, mRad, dRad, mJpi, dJpi, mE=mE, dE=dE, name=filename,
        phl=phl, logft=logft, robtdFile=robtdFile)
        self.log('Ini file written to %s.' % filename)
        self.setCurrentIniFile(filename)
        return filename

    def loadIniFile(self, filename = None):
        self.checkUnsavedTransitionChanges()

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
            try:
                self.ui.dsb_halflife.setValue(config.getfloat('Transition', 'PartialHalflife'))
            except ConfigParser.NoOptionError:
                pass
            try:
                self.ui.dsb_logft.setValue(config.getfloat('Transition', 'LogFt'))
            except ConfigParser.NoOptionError:
                pass
            try:
                self.setROBTDFile(config.get('Transition', 'ROBTDFile'))
            except ConfigParser.NoOptionError:
                pass
            self.setCurrentIniFile(filename)
            self.setTransitionLabel()
            self.log("Loaded Ini file: %s." % filename)
        except:
            QtGui.QErrorMessage(self).showMessage("Failed to load %s file" % filename)


    def writeConfigFile(self):
        filename = QtGui.QFileDialog.getSaveFileName(self, "Save config file")[0]
        if filename == '':
            return
        try:
            ut.writeConfigFile(filename, self.computationalCheckBoxes, self.computationalComboBoxes,
            self.computationalDSB, self.constants, self.spectrumCheckBoxes, self.spectrumDSB,
            self.spectrumNME, self.ui.cb_enforceNME.isChecked(), self.spectrumComboBoxes)
            self.setCurrentConfigFile(filename)
            self.log("Config file written to %s." % filename)
        except:
            QtGui.QErrorMessage(self).showMessage("Writing ini file failed.")

    def loadConfigFile(self, filename = None, nuclearStructure=True, spectrumShape=True):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self, "Choose config file")[0]
        if filename == '':
            return

        config = ConfigParser.RawConfigParser()
        config.read(filename)

        if nuclearStructure:
            for key in self.computationalDSB:
                try:
                    key.setValue(config.getfloat('Computational', self.computationalDSB[key]))
                except:
                    continue
            for key in self.computationalCheckBoxes:
                setCheckBoxState(key, False)
                try:
                    setCheckBoxState(key, config.getboolean('Computational', self.computationalCheckBoxes[key]))
                except:
                    continue
            for key in self.computationalComboBoxes:
                key.setCurrentIndex(0)
                try:
                    key.setCurrentIndex(key.findText(config.get('Computational', self.computationalComboBoxes[key])))
                except:
                    continue
            for key in self.constants:
                try:
                    key.setValue(config.getfloat('Constants', self.constants[key]))
                except:
                    continue

        if spectrumShape:
            for key in self.spectrumCheckBoxes:
                setCheckBoxState(key, True)
                try:
                    setCheckBoxState(key, config.getboolean('Spectrum', self.spectrumCheckBoxes[key]))
                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                    continue
            for key in self.spectrumDSB:
                try:
                    key.setValue(config.getfloat('Spectrum', self.spectrumDSB[key]))
                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                    continue
            for key in self.spectrumNME:
                key.setValue(0.)
                setCheckBoxState(self.ui.cb_enforceNME, False)
                try:
                    key.setValue(config.getfloat('Spectrum', self.spectrumNME[key]))
                    setCheckBoxState(self.ui.cb_enforceNME, True)
                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                    continue
            for key in self.spectrumComboBoxes:
                key.setCurrentIndex(0)
                try:
                    key.setCurrentIndex(key.findText(config.get('Spectrum', self.spectrumComboBoxes[key])))
                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
                    continue

        self.setCurrentConfigFile(filename)
        self.log("Loaded config file: %s." % filename)


if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    mw = BSG_UI()

    # setup stylesheet
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyside())
    mw.show()
    app.exec_()

    sys.exit()
