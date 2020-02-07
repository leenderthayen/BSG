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

from PyQt4 import QtCore, QtGui
# from PySide2.QtWidgets import QApplication, QMainWindow
# from PySide2 import QtGui

from ui.MainWindowGUI import Ui_MainWindow

import configparser

import numpy as np

import utils.utilities as ut

from managers_gui.input_mgr import InputManager
from managers_gui.config_mgr import ConfigManager
from managers_gui.database_mgr import DatabaseManager
from managers_gui.help_mgr import HelpManager

def setCheckBoxState(cb, b):
    if (cb.isChecked() and not b) or (not cb.isChecked() and b):
        cb.toggle()

class BSG_UI(QtGui.QMainWindow):

    def __init__(self):
        QtGui.QMainWindow.__init__(self, None)
        self.ui = Ui_MainWindow()

        self.ui.setupUi(self)

        self.inputMgr = InputManager(self)
        self.configMgr = ConfigManager(self)
        self.databaseMgr = DatabaseManager(self)
        self.helpMgr = HelpManager(self)

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

        self.ui.a_exit.triggered.connect(self.close)

        self.ui.a_clearPlots.triggered.connect(self.clearPlots)

        self.ui.dsb_bAc.valueChanged.connect(self.enforceMatrixElements)
        self.ui.dsb_dAc.valueChanged.connect(self.enforceMatrixElements)
        self.ui.dsb_lambda.valueChanged.connect(self.enforceMatrixElements)

        self.ui.cb_nmeMethod.currentIndexChanged[str].connect(self.setESPVisibility)

        self.ui.rb_stepSize.toggled.connect(self.ui.dsb_stepSize.setEnabled)
        self.ui.rb_steps.toggled.connect(self.ui.sb_steps.setEnabled)

        self.ui.a_new.triggered.connect(self.newDecay)
        self.ui.b_runBSG.clicked.connect(self.runBSG)


        self.plotColors = ('r','g','b','w','y')
        self.currentPlotIndex = 0

        self.log('Initialized...')

        self.statusBar().showMessage("Ready!")

    def log(self, message):
        self.ui.txt_lastActions.insertPlainText(message + '\n')
        self.ui.txt_lastActions.moveCursor(QtGui.QTextCursor.End)

    def status(self, message):
        self.statusBar().showMessage(message)

    def setESPVisibility(self, value):
        if value == 'ESP':
            self.ui.gb_ESP.setEnabled(True)
        else:
            self.ui.gb_ESP.setEnabled(False)

    def enforceMatrixElements(self):
        if not self.ui.cb_enforceNME.isChecked():
            self.ui.cb_enforceNME.toggle()

    def clearPlots(self):
        self.ui.gv_plotSpectrum.clear()
        self.currentPlotIndex = 0

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

    def setTransitionLabel(self):
        Am = self.ui.sb_AM.value()
        Ad = self.ui.sb_AD.value()
        Zm = self.ui.sb_ZM.value()
        Zd = self.ui.sb_ZD.value()
        JpiM = str(int(self.ui.dsb_JM.value())) if self.ui.dsb_JM.value() % 1 == 0. else str(int(2*self.ui.dsb_JM.value())) + '/2'
        JpiM += self.ui.cb_PiM.currentText()
        JpiD = str(int(self.ui.dsb_JD.value())) if self.ui.dsb_JD.value() % 1 == 0. else str(int(2*self.ui.dsb_JD.value())) + '/2'
        JpiD += self.ui.cb_PiD.currentText()
        s = '{}{} ([{}] {}) -> {}{} ([{}] {})'.format(Am, ut.atoms[Zm],
        JpiM, self.ui.dsb_motherEn.value(), Ad, ut.atoms[Zd],
        JpiD, self.ui.dsb_daughterEn.value())
        self.ui.l_transitionName.setText(s)
        self.ui.l_transitionName.setToolTip(s)
        self.ui.l_transitionName.setStatusTip(s)

    def newDecay(self):
        self.inputMgr.checkUnsavedTransitionChanges()

        isotope, ok = QtGui.QInputDialog.getText(self, 'Isotope', 'Enter the isotope name (e.g. 238U)')
        if not ok:
            return
        import re
        A = int(re.findall(r'\d+', isotope)[0])
        name = isotope.replace(str(A), '').capitalize()
        Z = ut.atoms.index(name)

        self.ui.sb_AM.setValue(A)
        self.ui.sb_ZM.setValue(Z)

        button = QtGui.QMessageBox.question(self, 'ENSDF search', 'Search for transitions in ENSDF database?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if button == QtGui.QMessageBox.Yes:
            success = self.databaseMgr.loadESNDFBranches(Z, A)
            if not success:
                return
        else:
            process = QtGui.QInputDialog.getItem(self, "Choose beta process", "Beta process", ('B-', 'B+', 'EC'), editable=False)[0]
            self.ui.cb_process.setCurrentIndex(self.ui.cb_process.findText(process))
            self.ui.sb_AD.setValue(A)
            self.ui.sb_ZD.setValue(Z+1 if process == 'B-' else Z-1)

        button = QtGui.QMessageBox.question(self, 'Deformation', 'Enter deformation info from the Moeller database?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if button == QtGui.QMessageBox.Yes:
            self.databaseMgr.loadDeformation()

        button = QtGui.QMessageBox.question(self, "Charge radii", 'Load charge radii from database?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if button == QtGui.QMessageBox.Yes:
            self.databaseMgr.loadRadii()

        button = QtGui.QMessageBox.question(self, "Save as", 'Save transition in ini file?', QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)

        if button == QtGui.QMessageBox.Yes:
            self.inputMgr.writeIniFile()

    def runBSG(self):
        self.inputMgr.checkUnsavedTransitionChanges()

        outputName = self.ui.le_outputName.text()
        if self.databaseMgr.execPath == '':
            QtGui.QErrorMessage(self).showMessage("Set the path for the generator executable.")
        if self.inputMgr.iniName == '':
            QtGui.QErrorMessage(self).showMessage("No ini file found.")
            return
        if self.ui.rb_configFile.isChecked() and self.configMgr.configName == '':
            QtGui.QErrorMessage(self).showMessage("No configuration file specified.")
            return
        if self.ui.cb_exchange.isChecked() and self.databaseMgr.exchangePath == '':
            QtGui.QErrorMessage(self).showMessage("Choose the Exchange Data file, or turn off atomic exchange.")
            return

        self.log("Performing BSG Calculation...")
        self.status("Calculating...")
        command = "{0} -i {1} -o {2}".format(self.databaseMgr.execPath, self.inputMgr.iniName, outputName)
        if self.ui.cb_exchange.isChecked():
            command += " -e {}".format(self.databaseMgr.exchangePath)

        if self.ui.rb_configFile.isChecked():
            command += " -c {}".format(self.configMgr.configName)
        else:
            for key in self.configMgr.spectrumCheckBoxes:
                command += " --Spectrum.{0}={1}".format(self.configMgr.spectrumCheckBoxes[key], key.isChecked())
            for key in self.configMgr.spectrumComboBoxes:
                command += " --Spectrum.{0}={1}".format(self.configMgr.spectrumComboBoxes[key], key.currentText())
            for key in self.configMgr.spectrumDSB:
                if key.isEnabled():
                    command += " --Spectrum.{0}={1}".format(self.configMgr.spectrumDSB[key], key.value())
            if self.ui.cb_enforceNME.isChecked():
                for key in self.configMgr.spectrumNME:
                    command += " --Spectrum.{0}={1}".format(self.configMgr.spectrumNME[key], key.value())
            for key in self.configMgr.computationalComboBoxes:
                command += ' --Computational.{0}={1}'.format(self.configMgr.computationalComboBoxes[key], key.currentText())
            for key in self.configMgr.computationalCheckBoxes:
                command += ' -- Computational.{0}={1}'.format(self.configMgr.computationalCheckBoxes[key], key.isChecked())
            for key in self.configMgr.computationalDSB:
                command += ' --Computational.{0}={1}'.format(self.configMgr.computationalDSB[key], key.value())

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

if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    mw = BSG_UI()

    # setup stylesheet
    app.setStyleSheet(qdarkstyle.load_stylesheet_pyqt())
    mw.show()
    app.exec_()

    sys.exit()
