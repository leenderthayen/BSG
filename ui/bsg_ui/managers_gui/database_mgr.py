#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""
from PyQt4 import QtCore, QtGui
# from PySide2.QtWidgets import QApplication, QMainWindow
# from PySide2 import QtGui
import os
import utils.utilities as ut
#from bsg_gui.ui.MainWindowGUI import Ui_MainWindow

class DatabaseManager:
    def __init__(self,bsg_ui):
        self.bsg_ui = bsg_ui
        self.deformationFile = None
        self.radiiFile = None
        self.ensdfFolder = None
        self.execPath = ''
        self.exchangePath = ''
        self.bsg_ui.ui.a_loadDeformation.triggered.connect(self.loadDeformation)
        self.bsg_ui.ui.a_loadRadii.triggered.connect(self.loadRadii)
        self.bsg_ui.ui.a_findTransition.triggered.connect(self.loadESNDFBranches)
        self.bsg_ui.ui.a_bsgExec.triggered.connect(self.changeBSGExec)
        self.bsg_ui.ui.a_exchangeData.triggered.connect(self.changeBSGExchange)
        self.findDefaults()

    def findDefaults(self):
        self.bsg_ui.log("Looking for defaults in environment variables...")
        bsg_exec = os.environ.get('BSGPATH')
        if bsg_exec:
            self.execPath = bsg_exec
            self.bsg_ui.log("Found bsg_exec")
        exchangePath = os.environ.get('BSGEXCHANGEPATH')
        if exchangePath:
            self.exchangePath = exchangePath
            self.bsg_ui.log("Found Exchange data")
        ensdf = os.environ.get('ENSDFDIR')
        if ensdf:
            self.ensdfFolder = ensdf
            self.bsg_ui.log("Found ENSDF database")
        frdm = os.environ.get('FRDMPATH')
        if frdm:
            self.deformationFile = frdm
            self.bsg_ui.log("Found FRDM")
        chargeRadii = os.environ.get('CHARGERADIIPATH')
        if chargeRadii:
            self.radiiFile = chargeRadii
            self.bsg_ui.log("Found Charge Radii")


    def loadESNDFBranches(self, Z, A):
        if not self.ensdfFolder:
            self.ensdfFolder = QtGui.QFileDialog.getExistingDirectory(self.bsg_ui, "Choose the ENSDF directory", ".")
            if self.ensdfFolder == '':
                self.ensdfFolder = None
                return
        fullPath = self.ensdfFolder + '/ensdf.{:03d}'.format(A)
        print(fullPath)
        bbs = ut.getENSDFBetaBranches(str(fullPath), Z, A, 0)
        print(bbs)
        items = []
        for i in bbs:
            s = '{}: {}{} ([{}] {}) -> {}{} ([{}] {}); Endpoint: {} keV, t1/2: {:.3f} s'.format(i.process, A, ut.atoms[Z], i.motherLevel.Jpi, i.motherLevel.E,\
            A, ut.atoms[Z+1 if i.process == 'B-' else Z-1], i.daughterLevel.Jpi, i.daughterLevel.E, i.E, i.partialHalflife)

            items.append(s)

        if len(items):
            choice, ok  = QtGui.QInputDialog.getItem(self.bsg_ui, "Transition", "Choose the beta transition:", items, editable=False)
            if ok:
                branch = bbs[items.index(choice)]

                self.bsg_ui.ui.sb_ZM.setValue(Z)
                self.bsg_ui.ui.sb_AM.setValue(A)
                self.bsg_ui.ui.sb_ZD.setValue(Z+1 if branch.process == 'B-' else Z-1)
                self.bsg_ui.ui.sb_AD.setValue(A)
                self.bsg_ui.ui.cb_process.setCurrentIndex(self.bsg_ui.ui.cb_process.findText(branch.process))
                self.bsg_ui.ui.cb_type.setCurrentIndex(1)
                self.bsg_ui.ui.dsb_motherEn.setValue(branch.motherLevel.E)
                self.bsg_ui.ui.dsb_daughterEn.setValue(branch.daughterLevel.E)
                self.bsg_ui.ui.dsb_Q.setValue(branch.E-branch.motherLevel.E+branch.daughterLevel.E)
                mJpi = branch.motherLevel.Jpi
                self.bsg_ui.ui.dsb_JM.setValue(float(mJpi.split('/2')[0])/2. if '/2' in mJpi else float(mJpi[:-1]))
                dJpi = branch.daughterLevel.Jpi
                self.bsg_ui.ui.dsb_JD.setValue(float(dJpi.split('/2')[0])/2. if '/2' in dJpi else float(dJpi[:-1]))
                self.bsg_ui.ui.cb_PiM.setCurrentIndex(self.bsg_ui.ui.cb_PiM.findText(mJpi[-1]))
                self.bsg_ui.ui.cb_PiD.setCurrentIndex(self.bsg_ui.ui.cb_PiD.findText(dJpi[-1]))
                self.bsg_ui.ui.dsb_halflife.setValue(branch.partialHalflife)
                self.bsg_ui.ui.dsb_logft.setValue(branch.logft)

                self.bsg_ui.setTransitionLabel()
                return True
        else:
            QtGui.QErrorMessage(self.bsg_ui).showMessage("No transitions found.")
        return False

    def loadDeformation(self):
        if not self.deformationFile:
            self.deformationFile = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose the deformation file")[0]
        if self.deformationFile == '':
            return
        mDef = ut.loadDeformation(self.bsg_ui.ui.sb_ZM.value(), self.bsg_ui.ui.sb_AM.value(), self.deformationFile)
        dDef = ut.loadDeformation(self.bsg_ui.ui.sb_ZD.value(), self.bsg_ui.ui.sb_AD.value(), self.deformationFile)

        self.bsg_ui.ui.dsb_Beta2M.setValue(mDef[0])
        self.bsg_ui.ui.dsb_Beta4M.setValue(mDef[2])
        self.bsg_ui.ui.dsb_Beta6M.setValue(mDef[3])
        self.bsg_ui.ui.dsb_Beta2D.setValue(dDef[0])
        self.bsg_ui.ui.dsb_Beta4D.setValue(dDef[2])
        self.bsg_ui.ui.dsb_Beta6D.setValue(dDef[3])

    def loadRadii(self):
        if not self.radiiFile:
            self.radiiFile = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose the radii file")[0]
        if self.radiiFile == '':
            return
        mRad = ut.loadChargeRadius(self.bsg_ui.ui.sb_ZM.value(), self.bsg_ui.ui.sb_AM.value(), self.radiiFile)
        dRad = ut.loadChargeRadius(self.bsg_ui.ui.sb_ZD.value(), self.bsg_ui.ui.sb_AD.value(), self.radiiFile)

        self.bsg_ui.ui.dsb_RM.setValue(mRad)
        self.bsg_ui.ui.dsb_RD.setValue(dRad)

    def changeBSGExec(self):
        filename = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose BSG exec")[0]
        if filename == '':
            return
        self.execPath = filename

    def changeBSGExchange(self):
        filename = QtGui.QFileDialog.getOpenFileName(self.bsg_ui, "Choose Exchange data file")[0]
        if filename == '':
            return
        self.exchangePath = filename
