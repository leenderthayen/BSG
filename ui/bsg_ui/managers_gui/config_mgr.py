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

        self.couplingOptions = {"gV": (1.0, self.bsg_ui.ui.dsb_gAeff), "gA": (1.2723, self.bsg_ui.ui.dsb_gA),\
        "gM": (4.706, self.bsg_ui.ui.dsb_gM), "gP": (0.0, self.bsg_ui.ui.dsb_gP)}

        self.spectrumOptions = {"phasespace": (True, self.bsg_ui.ui.cb_phasespace), "fermi": (True,self.bsg_ui.ui.cb_fermi),\
        "esfinitesize": (True, self.bsg_ui.ui.cb_efs), "shapefactor": (True, self.bsg_ui.ui.cb_C), "isovector": (True, self.bsg_ui.ui.cb_isovector),\
        "relativistic": (True, self.bsg_ui.ui.cb_relativistic), "esdeformation": (True, self.bsg_ui.ui.cb_efs_deformation),\
        "esfermi": (True, self.bsg_ui.ui.cb_ens), "coulombrecoil": (True, self.bsg_ui.ui.cb_Q), "radiative": (True, self.bsg_ui.ui.cb_radiative),\
        "kinematicrecoil": (True, self.bsg_ui.ui.cb_recoil), "screening": (True, self.bsg_ui.ui.cb_screening), "exchange": (True, self.bsg_ui.ui.cb_exchange),\
        "atomicmismatch": (True, self.bsg_ui.ui.cb_mismatch)}

        self.calculationOptions = {"Begin": (0.0, self.bsg_ui.ui.dsb_begin), "End": (0.0, self.bsg_ui.ui.dsb_end),\
        "StepSize": (0.10, self.bsg_ui.ui.dsb_stepSize), "Steps": (0.0, self.bsg_ui.ui.sb_steps)}#missing: neutrino

        self.advancedOptions = {"Connect": (True, self.bsg_ui.ui.cb_coupled), "ESShape": ("Fermi", self.bsg_ui.ui.cb_shape),\
        "NSShape": ("UCS", self.bsg_ui.ui.cb_C_shape)}#missing: vold and vnew

        self.matrixElementsOptions = {"WeakMagnetism": (0.0, self.bsg_ui.ui.dsb_bAc), "InducedTensor": (0.0, self.bsg_ui.ui.dsb_dAc),\
        "Lambda": (0.0, self.bsg_ui.ui.dsb_lambda)}

        self.config_settings = {"Coupling": self.couplingOptions, "Spectrum": self.spectrumOptions, "Calculation": self.calculationOptions,\
        "Advanced": self.advancedOptions, "MatrixElements": self.matrixElementsOptions}

        ########################################
        #at the moment these are not used by BSG
        self.spectrumCheckBoxes =  {"CDeformation": (True, self.bsg_ui.ui.cb_C_deformation),\
        "EnableNrOfBins": (False, self.bsg_ui.ui.rb_steps), "EnforceNME": (False,self.bsg_ui.ui.cb_enforceNME)}

        self.computationalComboBoxes = {"Method": ("ROBTD", self.bsg_ui.ui.cb_nmeMethod), "Potential": ("SHO", self.bsg_ui.ui.cb_potential)}

        self.computationalCheckBoxes = {"ForceSpin": (True, self.bsg_ui.ui.cb_forceSpin),
        "ReversedGhallagher": (False, self.bsg_ui.ui.cb_reverseGhallagher), "OverrideSPCoupling": (False, self.bsg_ui.ui.cb_overrideCoupling)}

        self.computationalDSB = {"EnergyMargin": (1.0, self.bsg_ui.ui.dsb_energyMargin),
        "Vneutron": (49.60, self.bsg_ui.ui.dsb_V0n), "Vproton": (49.60, self.bsg_ui.ui.dsb_V0p), "V0Sneutron": (7.20, self.bsg_ui.ui.dsb_V0Sn), "V0Sproton": (7.20, self.bsg_ui.ui.dsb_V0Sp), "Xneutron": (0.86, self.bsg_ui.ui.dsb_Xn), "Xproton": (0.86, self.bsg_ui.ui.dsb_Xp)}
        ########################################


    def loadConfigFile(self, filename = None, nuclearStructure=True, spectrumShape=True):
        if not filename:
            filename = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose config file")
        if filename == '':
            return

        config = configparser.ConfigParser()
        config.read(filename)



        for conf_key in self.config_settings:
            for key in self.config_settings[conf_key]:
                try:
                    t = self.config_settings[conf_key][key]
                    if isinstance(t[0],bool):
                        newt0 = config.getboolean(conf_key, key)
                        t[1].setChecked(False)
                        t[1].setChecked(newt0)
                    elif isinstance(t[0],int):
                        newt0 = config.getint(conf_key, key)
                        t[1].setValue(newt0)
                    elif isinstance(t[0],float):
                        newt0 = config.getfloat(conf_key, key)
                        t[1].setValue(newt0)
                    else:
                        newt0 = config.get(conf_key, key)
                        t[1].setCurrentIndex(0)
                        t[1].setCurrentIndex(t[1].findText(newt0))
                    self.config_settings[conf_key][key] = (newt0,t[1])
                except (configparser.NoSectionError, configparser.NoOptionError):
                    continue

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
                    newt = t[1].isChecked()
                elif isinstance(t[0],int):
                    newt = t[1].value()
                elif isinstance(t[0],float):
                    newt = t[1].value()
                else:
                    newt = t[1].currentText()
                self.config_settings[conf_key][key] = (newt,t[1])

    def enableAll(self):
        for conf_key in self.config_settings:
            for key in self.config_settings[conf_key]:
                t = self.config_settings[conf_key][key]
                if isinstance(t[0],bool):
                    t[1].setChecked(True)

        self.updateConfigData()

    def disableAll(self):
        for conf_key in self.config_settings:
            for key in self.config_settings[conf_key]:
                t = self.config_settings[conf_key][key]
                if isinstance(t[0],bool):
                    t[1].setChecked(False)

        self.updateConfigData()


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
        print(self.configName)
        self.loadConfigFile(filename=self.configName, nuclearStructure=False)
        self.bsg_ui.log('Resetting spectrum shape options back to values from %s.' % self.configName)

    def resetNuclearStructureOptions(self):
        self.loadConfigFile(filename=self.configName, spectrumShape=False)
        self.bsg_ui.log('Resetting nuclear structure options back to values from %s.' % self.configName)
