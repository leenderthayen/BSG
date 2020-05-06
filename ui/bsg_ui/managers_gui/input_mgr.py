#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""
from PyQt5 import QtCore, QtGui
# from PySide2.QtWidgets import QApplication, QMainWindow
# from PySide2 import QtGui
#from bsg_gui.ui.MainWindowGUI import Ui_MainWindow
#import bsg_gui.utils.utilities as ut
import configparser
import utils.utilities as ut

class InputManager:
    def __init__(self,bsg_ui):
        self.bsg_ui = bsg_ui
        self.iniName = ''
        self.robtdFile = ''
        self.unsavedTransitionChanges = list()
        self.bsg_ui.ui.a_loadIni.triggered.connect(self.loadIniFile)
        self.bsg_ui.ui.a_saveIni.triggered.connect(self.writeIniFile)
        self.bsg_ui.ui.b_chooseROBTDFile.clicked.connect(self.setROBTDFile)

    def loadIniFile(self, filename = None):
        self.checkUnsavedTransitionChanges()

        if not filename:
            filename,ok = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose .ini file")
        if not ok:
            return
        try:
            config = configparser.ConfigParser()
            config.read(filename)

            self.bsg_ui.ui.cb_process.setCurrentIndex(self.bsg_ui.ui.cb_process.findText(config.get('Transition', 'Process')))
            self.bsg_ui.ui.cb_type.setCurrentIndex(self.bsg_ui.ui.cb_type.findText(config.get('Transition', 'Type')))
            self.bsg_ui.ui.dsb_MixingRatio.setValue(config.getfloat('Transition', 'MixingRatio'))
            self.bsg_ui.ui.dsb_Q.setValue(config.getfloat('Transition', 'QValue'))
            self.bsg_ui.ui.sb_ZM.setValue(config.getint('Mother', 'Z'))
            self.bsg_ui.ui.sb_AM.setValue(config.getint('Mother', 'A'))
            self.bsg_ui.ui.dsb_RM.setValue(config.getfloat('Mother', 'Radius'))
            self.bsg_ui.ui.dsb_JM.setValue(abs(config.getint('Mother', 'SpinParity'))/2.)
            self.bsg_ui.ui.cb_PiM.setCurrentIndex(0 if config.getint('Mother', 'SpinParity') > 0 else 1)
            self.bsg_ui.ui.dsb_Beta2M.setValue(config.getfloat('Mother', 'Beta2'))
            self.bsg_ui.ui.dsb_Beta4M.setValue(config.getfloat('Mother', 'Beta4'))
            self.bsg_ui.ui.dsb_Beta6M.setValue(config.getfloat('Mother', 'Beta6'))
            self.bsg_ui.ui.sb_ZD.setValue(config.getint('Daughter', 'Z'))
            self.bsg_ui.ui.sb_AD.setValue(config.getint('Daughter', 'A'))
            self.bsg_ui.ui.dsb_RD.setValue(config.getfloat('Daughter', 'Radius'))
            self.bsg_ui.ui.dsb_JD.setValue(abs(config.getint('Daughter', 'SpinParity'))/2.)
            self.bsg_ui.ui.cb_PiD.setCurrentIndex(0 if config.getint('Daughter', 'SpinParity') > 0 else 1)
            self.bsg_ui.ui.dsb_Beta2D.setValue(config.getfloat('Daughter', 'Beta2'))
            self.bsg_ui.ui.dsb_Beta4D.setValue(config.getfloat('Daughter', 'Beta6'))
            self.bsg_ui.ui.dsb_Beta6D.setValue(config.getfloat('Daughter', 'Beta6'))
            try:
                self.bsg_ui.ui.dsb_halflife.setValue(config.getfloat('Transition', 'PartialHalflife'))
            except configparser.NoOptionError:
                pass
            try:
                self.bsg_ui.ui.dsb_logft.setValue(config.getfloat('Transition', 'LogFt'))
            except configparser.NoOptionError:
                pass
            try:
                self.setROBTDFile(config.get('Transition', 'ROBTDFile'))
            except configparser.NoOptionError:
                pass
            self.setCurrentIniFile(filename)
            self.bsg_ui.setTransitionLabel()
            self.bsg_ui.log("Loaded Ini file: %s." % filename)
        except:
            QtGui.QErrorMessage(self.bsg_ui).showMessage("Failed to load %s file" % filename)

    def writeIniFile(self):
        Zm = self.bsg_ui.ui.sb_ZM.value()
        Zd = self.bsg_ui.ui.sb_ZD.value()
        Am = self.bsg_ui.ui.sb_AM.value()
        Ad = self.bsg_ui.ui.sb_AD.value()

        process = self.bsg_ui.ui.cb_process.currentText()

        if Am != Ad:
            QtGui.QErrorMessage(self.bsg_ui).showMessage("Mass numbers do not agree.")
            return
        if Zd != (Zm + 1 if process == "B-" else Zm -1):
            QtGui.QErrorMessage(self.bsg_ui).showMessage("Z values do not agree with process")
            return

        Q = self.bsg_ui.ui.dsb_Q.value()
        decayType = self.bsg_ui.ui.cb_type.currentText()
        mJ = self.bsg_ui.ui.dsb_JM.value()
        dJ = self.bsg_ui.ui.dsb_JD.value()

        mRad = self.bsg_ui.ui.dsb_RM.value()
        dRad = self.bsg_ui.ui.dsb_RD.value()

        beta2m = self.bsg_ui.ui.dsb_Beta2M.value()
        beta4m = self.bsg_ui.ui.dsb_Beta4M.value()
        beta6m = self.bsg_ui.ui.dsb_Beta6M.value()
        beta2d = self.bsg_ui.ui.dsb_Beta2D.value()
        beta4d = self.bsg_ui.ui.dsb_Beta4D.value()
        beta6d = self.bsg_ui.ui.dsb_Beta6D.value()

        mE = self.bsg_ui.ui.dsb_motherEn.value()
        dE = self.bsg_ui.ui.dsb_daughterEn.value()

        mJpi = int(mJ*2)*(1 if self.bsg_ui.ui.cb_PiM.currentText() == "+" else -1)
        dJpi = int(dJ*2)*(1 if self.bsg_ui.ui.cb_PiD.currentText() == "+" else -1)

        phl = self.bsg_ui.ui.dsb_halflife.value()
        phl = phl if phl > 0 else None

        logft = self.bsg_ui.ui.dsb_logft.value()
        logft = logft if logft > 0 else None

        robtdFile = self.robtdFile if not self.robtdFile == '' else None

        filename,_ = QtGui.QFileDialog.getSaveFileName(self.bsg_ui, "Save .ini file")
        if filename == '':
            return
        print(filename)
        ut.writeIniFile(Zm, Zd, Am, Q, process, decayType, beta2m, beta4m, beta6m,beta2d, beta4d, beta6d, mRad, dRad, mJpi, dJpi, mE=mE, dE=dE, name=filename,phl=phl, logft=logft, robtdFile=robtdFile)
        self.bsg_ui.log('Ini file written to %s.' % filename)
        self.setCurrentIniFile(filename)
        return filename

    def setCurrentIniFile(self, filename):
        self.iniName = filename
        self.bsg_ui.ui.l_iniFilename.setText(filename.split('/')[-1])
        self.bsg_ui.ui.l_iniFilename.setToolTip(filename)
        self.bsg_ui.ui.l_iniFilename.setStatusTip(filename)

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
            button = QtGui.QMessageBox.question(self.bsg_ui, 'Unsaved changes',
            s, QtGui.QMessageBox.Yes, QtGui.QMessageBox.No)
            if button == QtGui.QMessageBox.Yes:
                self.writeIniFile()

    def checkUnsavedTransitionChanges(self):
        self.clearUnsavedChanges()
        if self.iniName == '':
            self.addUnsavedChanges("No information currently saved.")
        else:
            config = configparser.ConfigParser()
            config.read(self.iniName)

            if config.get('Transition', 'Process') != self.bsg_ui.ui.cb_process.currentText():
                self.addUnsavedChanges('Process (Prev. %s)' % config.get('Transition', 'Process'))
            if config.get('Transition', 'Type') != self.bsg_ui.ui.cb_type.currentText():
                self.addUnsavedChanges('Type (Prev. %s)' % config.get('Transition', 'Type'))
            if config.getfloat('Transition', 'MixingRatio') != self.bsg_ui.ui.dsb_MixingRatio.value():
                self.addUnsavedChanges('Mixing Ratio (Prev. %s)' % config.get('Transition', 'MixingRatio'))
            if config.getfloat('Transition', 'QValue') != self.bsg_ui.ui.dsb_Q.value():
                self.addUnsavedChanges('Q Value (Prev. %s)' % config.get('Transition', 'QValue'))
            try:
                if config.getfloat('Transition', 'PartialHalflife') != self.bsg_ui.ui.dsb_halflife.value():
                    self.addUnsavedChanges('Partial Halflife (Prev. %s)' % config.get('Transition', 'PartialHalflife'))
            except configparser.NoOptionError:
                pass
            try:
                if config.getfloat('Transition', 'LogFt') != self.bsg_ui.ui.dsb_logft.value():
                    self.addUnsavedChanges('Log Ft (Prev. %s)' % config.get('Transition', 'LogFt'))
            except configparser.NoOptionError:
                pass

            if config.getint('Mother', 'Z') != self.bsg_ui.ui.sb_ZM.value():
                self.addUnsavedChanges('Mother Z (Prev. %d)' % config.getint('Mother', 'Z'))
            if config.getint('Mother', 'A') != self.bsg_ui.ui.sb_AM.value():
                self.addUnsavedChanges('Mother A (Prev. %d)' % config.getint('Mother', 'A'))
            if config.getfloat('Mother', 'Radius') != self.bsg_ui.ui.dsb_RM.value():
                self.addUnsavedChanges('Mother Radius (Prev. %.2f)' % config.getfloat('Mother', 'Radius'))
            if abs(config.getint('Mother', 'SpinParity')/2.) != self.bsg_ui.ui.dsb_JM.value():
                self.addUnsavedChanges('Mother Spin (Prev. %.1f)' % abs(config.getint('Mother', 'SpinParity')/2.))
            if None:
                pass
            if config.getfloat('Mother', 'Beta2') != self.bsg_ui.ui.dsb_Beta2M.value():
                self.addUnsavedChanges('Mother beta2 (Prev. %.3f)' % config.getfloat('Mother', 'Beta2'))
            if config.getfloat('Mother', 'Beta4') != self.bsg_ui.ui.dsb_Beta4M.value():
                self.addUnsavedChanges('Mother beta4 (Prev. %.3f)' % config.getfloat('Mother', 'Beta4'))
            if config.getfloat('Mother', 'Beta6') != self.bsg_ui.ui.dsb_Beta6M.value():
                self.addUnsavedChanges('Mother beta6 (Prev. %.3f)' % config.getfloat('Mother', 'Beta6'))

            if config.getint('Daughter', 'Z') != self.bsg_ui.ui.sb_ZD.value():
                self.addUnsavedChanges('Daughter Z (Prev. %d)' % config.getint('Daughter', 'Z'))
            if config.getint('Daughter', 'A') != self.bsg_ui.ui.sb_AD.value():
                self.addUnsavedChanges('Daughter A (Prev. %d)' % config.getint('Daughter', 'A'))
            if config.getfloat('Daughter', 'Radius') != self.bsg_ui.ui.dsb_RD.value():
                self.addUnsavedChanges('Daughter Radius (Prev. %.2f)' % config.getfloat('Daughter', 'Radius'))
            if abs(config.getint('Daughter', 'SpinParity')/2.) != self.bsg_ui.ui.dsb_JD.value():
                self.addUnsavedChanges('Daughter Spin (Prev. %.1f)' % abs(config.getint('Daughter', 'SpinParity')/2.))
            if None:
                pass
            if config.getfloat('Daughter', 'Beta2') != self.bsg_ui.ui.dsb_Beta2D.value():
                self.addUnsavedChanges('Daughter beta2 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta2'))
            if config.getfloat('Daughter', 'Beta4') != self.bsg_ui.ui.dsb_Beta4D.value():
                self.addUnsavedChanges('Daughter beta4 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta4'))
            if config.getfloat('Daughter', 'Beta6') != self.bsg_ui.ui.dsb_Beta6D.value():
                self.addUnsavedChanges('Daughter beta6 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta6'))

        self.askSaveChanges()


    def setROBTDFile(self, filename=None):
        if not filename:
            filename,ok = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose ROBTD data file")
            if not ok:
                return
        self.robtdFile = filename
        self.bsg_ui.ui.b_chooseROBTDFile.setText(filename.split('/')[-1])
        self.bsg_ui.ui.b_chooseROBTDFile.setToolTip(filename)
        self.bsg_ui.ui.b_chooseROBTDFile.setStatusTip(filename)
