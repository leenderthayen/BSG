#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""
from PyQt4 import QtCore, QtGui
# from PySide2.QtWidgets import QApplication, QMainWindow
# from PySide2 import QtGui
#from bsg_gui.ui.MainWindowGUI import Ui_MainWindow
import configparser
import utils.utilities as ut

class ConfigManager:
    def __init__(self,bsg_ui):
        self.bsg_ui = bsg_ui
        self.configName = ''
        self.bsg_ui.ui.a_loadConfig.triggered.connect(self.loadConfigFile)
        self.bsg_ui.ui.a_saveConfig.triggered.connect(self.writeConfigFile)
        self.bsg_ui.ui.b_chooseConfigFile.clicked.connect(self.loadConfigFile)
        self.bsg_ui.ui.b_resetSSO.clicked.connect(self.resetSpectrumShapeOptions)
        self.bsg_ui.ui.b_resetNSO.clicked.connect(self.resetNuclearStructureOptions)


        self.spectrumCheckBoxes =  {"Phasespace": (True, self.bsg_ui.ui.cb_phasespace), "Fermi": (True,self.bsg_ui.ui.cb_fermi), "Radiative": (True, self.bsg_ui.ui.cb_radiative),\
        "ESFiniteSize": (True, self.bsg_ui.ui.cb_efs), "ESShape": (True, self.bsg_ui.ui.cb_ens), "ESDeformation": (True, self.bsg_ui.ui.cb_efs_deformation),\
        "NSShape": (True, self.bsg_ui.ui.cb_C), "Isovector": (True, self.bsg_ui.ui.cb_isovector), "Connect": (True, self.bsg_ui.ui.cb_coupled), "CDeformation": (True, self.bsg_ui.ui.cb_C_deformation),\
        "Relativistic": (True, self.bsg_ui.ui.cb_relativistic), "Recoil": (True, self.bsg_ui.ui.cb_recoil), "CoulombRecoil": (True, self.bsg_ui.ui.cb_Q), \
        "Screening": (True, self.bsg_ui.ui.cb_screening), "Exchange": (True, self.bsg_ui.ui.cb_exchange), "AtomicMismatch": (True, self.bsg_ui.ui.cb_mismatch), "EnableNrOfBins": (False, self.bsg_ui.ui.rb_steps), "EnforceNME": (False,self.bsg_ui.ui.cb_enforceNME)}

        self.spectrumNME = {"WeakMagnetism": (0.0, self.bsg_ui.ui.dsb_bAc), "InducedTensor": (0.0, self.bsg_ui.ui.dsb_dAc),\
        "Lambda": (0.0, self.bsg_ui.ui.dsb_lambda)}

        self.spectrumDSB = {"Begin": (0.0, self.bsg_ui.ui.dsb_begin), "End": (0.0, self.bsg_ui.ui.dsb_end),\
        "StepSize": (0.10, self.bsg_ui.ui.dsb_stepSize), "Steps": (0.0, self.bsg_ui.ui.sb_steps)}

        self.spectrumComboBoxes = {"ESShapeForm": ("Fermi", self.bsg_ui.ui.cb_shape), "NSShapeForm": ("UCS", self.bsg_ui.ui.cb_C_shape)}

        self.spectrum = self.spectrumCheckBoxes
        self.spectrum.update(self.spectrumComboBoxes)
        self.spectrum.update(self.spectrumNME)
        self.spectrum.update(self.spectrumDSB)

        self.constants = {"gA": (1.2723, self.bsg_ui.ui.dsb_gA), "gAeff": (1.0, self.bsg_ui.ui.dsb_gAeff), "gP": (0.0, self.bsg_ui.ui.dsb_gP),\
        "gM": (4.706, self.bsg_ui.ui.dsb_gM)}

        self.computationalComboBoxes = {"Method": ("ROBTD", self.bsg_ui.ui.cb_nmeMethod), "Potential": ("SHO", self.bsg_ui.ui.cb_potential)}

        self.computationalCheckBoxes = {"ForceSpin": (True, self.bsg_ui.ui.cb_forceSpin),
        "ReversedGhallagher": (False, self.bsg_ui.ui.cb_reverseGhallagher), "OverrideSPCoupling": (False, self.bsg_ui.ui.cb_overrideCoupling)}

        self.computationalDSB = {"EnergyMargin": (1.0, self.bsg_ui.ui.dsb_energyMargin),
        "Vneutron": (49.60, self.bsg_ui.ui.dsb_V0n), "Vproton": (49.60, self.bsg_ui.ui.dsb_V0p), "V0Sneutron": (7.20, self.bsg_ui.ui.dsb_V0Sn), "V0Sproton": (7.20, self.bsg_ui.ui.dsb_V0Sp), "Xneutron": (0.86, self.bsg_ui.ui.dsb_Xn), "Xproton": (0.86, self.bsg_ui.ui.dsb_Xp)}

        self.computational = self.computationalCheckBoxes
        self.computational.update(self.computationalComboBoxes)
        self.computational.update(self.computationalDSB)

        self.general = {"Folder": '.'}

        self.config_settings = {"General": self.general, "Spectrum": self.spectrum, "Constants": self.constants, "Computational": self.computational}


    def loadConfigFile(self, filename = None, nuclearStructure=True, spectrumShape=True):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose config file")[0]
        if filename == '':
            return

        config = configparser.ConfigParser()
        config.read(filename)

        for conf_key in self.config_settings:
            for key in self.config_settings[conf_key]:
                try:
                    t = self.config_settings[conf_key][key]
                    if isinstance(t[0],bool):
                        t[0] = config.getboolean(conf_key, key)
                        setCheckBoxState(t[1], False)
                        setCheckBoxState(t[1], t[0])
                    elif isinstance(t[0],int):
                        t[0] = config.getint(conf_key, key)
                        t[1].setValue(t[0])
                    elif isinstance(t[0],float):
                        t[0] = config.getfloat(conf_key, key)
                        t[1].setValue(t[0])
                    else:
                        t[0] = config.get(conf_key, key)
                        t[1].setCurrentIndex(0)
                        t[1].setCurrentIndex(t[1].findText(t[0]))
                except (configparser.NoSectionError, configparser.NoOptionError):
                    continue

#        if nuclearStructure:
#            for key in self.computationalDSB:
#                try:
#                    key.setValue(config.getfloat('Computational', self.computationalDSB[key]))
#                except:
#                    continue
#            for key in self.computationalCheckBoxes:
#                setCheckBoxState(key, False)
#                try:
#                    setCheckBoxState(key, config.getboolean('Computational', self.computationalCheckBoxes[key]))
#                except:
#                    continue
#            for key in self.computationalComboBoxes:
#                key.setCurrentIndex(0)
#                try:
#                    key.setCurrentIndex(key.findText(config.get('Computational', self.computationalComboBoxes[key])))
#                except:
#                    continue
#            for key in self.constants:
#                try:
#                    key.setValue(config.getfloat('Constants', self.constants[key]))
#                except:
#                    continue
#
#        if spectrumShape:
#            for key in self.spectrumCheckBoxes:
#                setCheckBoxState(key, True)
#                try:
#                    setCheckBoxState(key, config.getboolean('Spectrum', self.spectrumCheckBoxes[key]))
#                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
#                    continue
#            for key in self.spectrumDSB:
#                try:
#                    key.setValue(config.getfloat('Spectrum', self.spectrumDSB[key]))
#                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
#                    continue
#            for key in self.spectrumNME:
#                key.setValue(0.)
#                setCheckBoxState(self.ui.cb_enforceNME, False)
#                try:
#                    key.setValue(config.getfloat('Spectrum', self.spectrumNME[key]))
#                    setCheckBoxState(self.ui.cb_enforceNME, True)
#                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
#                    continue
#            for key in self.spectrumComboBoxes:
#                key.setCurrentIndex(0)
#                try:
#                    key.setCurrentIndex(key.findText(config.get('Spectrum', self.spectrumComboBoxes[key])))
#                except (ConfigParser.NoSectionError, ConfigParser.NoOptionError):
#                    continue

        self.setCurrentConfigFile(filename)
        self.bsg_ui.log("Loaded config file: %s." % filename)

    def writeConfigFile(self):
        self.updateConfigData()
        filename = QtGui.QFileDialog.getSaveFileName(self.bsg_ui, "Save config file")
        if filename == '':
            return
        try:
#            ut.writeConfigFile(filename, self.computationalCheckBoxes, self.computationalComboBoxes,
#            self.computationalDSB, self.constants, self.spectrumCheckBoxes, self.spectrumDSB,
#            self.spectrumNME, self.ui.cb_enforceNME.isChecked(), self.spectrumComboBoxes)
            ut.writeConfigFile(filename, self.config_settings)
            self.setCurrentConfigFile(filename)
            self.bsg_ui.log("Config file written to %s." % filename)
        except:
            QtGui.QErrorMessage(self.bsg_ui).showMessage("Writing config file failed.")

    def updateConfigData(self):
        for conf_key in self.config_settings:
            for key in self.config_settings[conf_key]:
                t = self.config_settings[conf_key][key]
                if isinstance(t[0],bool):
                    t[0] = t[1].isChecked()
                elif isinstance(t[0],int):
                    t[0] = t[1].value()
                elif isinstance(t[0],float):
                    t[0] = t[1].value()
                else:
                    t[0] = t[1].currentText()

    def setCurrentConfigFile(self, filename):
        self.configName = filename
        self.bsg_ui.ui.b_chooseConfigFile.setText(filename.split('/')[-1])
        self.bsg_ui.ui.b_chooseConfigFile.setToolTip(filename)
        self.bsg_ui.ui.b_chooseConfigFile.setStatusTip(filename)

        self.bsg_ui.ui.b_resetSSO.setText("Reset to %s" % filename.split('/')[-1])
        self.bsg_ui.ui.b_resetNSO.setText("Reset to %s" % filename.split('/')[-1])

        if not self.bsg_ui.ui.rb_configFile.isChecked():
            self.bsg_ui.ui.rb_configFile.toggle()

    def resetSpectrumShapeOptions(self):
        self.loadConfigFile(self.configName, nuclearStructure=False)
        self.bsg_ui.log('Resetting spectrum shape options back to values from %s.' % self.configName)

    def resetNuclearStructureOptions(self):
        self.loadConfigFile(self.configName, spectrumShape=False)
        self.bsg_ui.log('Resetting nuclear structure options back to values from %s.' % self.configName)
