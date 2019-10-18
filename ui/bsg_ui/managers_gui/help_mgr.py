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

class HelpManager:
    def __init__(self,bsg_ui):
        self.bsg_ui = bsg_ui
        self.bsg_ui.ui.a_ENSDF.triggered.connect(self.visitENSDF)
        self.bsg_ui.ui.a_FRDM2012.triggered.connect(self.visitFRDM2012)
        self.bsg_ui.ui.a_ChargeRadii.triggered.connect(self.visitChargeRadii)
        self.bsg_ui.ui.a_about.triggered.connect(self.about)
        self.bsg_ui.ui.a_feedback.triggered.connect(self.submitFeedback)


    def visitENSDF(self):
        QtGui.QDesktopServices.openUrl('https://www.nndc.bnl.gov/ensarchivals/')

    def visitFRDM2012(self):
        QtGui.QDesktopServices.openUrl('https://www.sciencedirect.com/science/article/pii/S0092640X1600005X')

    def visitChargeRadii(self):
        QtGui.QDesktopServices.openUrl('https://journals.aps.org/prc/abstract/10.1103/PhysRevC.94.064315')

    def about(self):
        QtGui.QMessageBox.information(self, "About", "Beta spectrum Generator GUI v 1.0")

    def submitFeedback(self):
        QtGui.QMessageBox.information(self, "Send feedback", "Send emails to leendert.hayen@gmail.com")
