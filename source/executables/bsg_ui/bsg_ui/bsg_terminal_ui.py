#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:39 2014

@author: leendert
"""
import sys
from shell import Shell, CommandError
#
#import qdarkstyle
#import os
#
#from PySide2.QtWidgets import QApplication, QMainWindow
#from PySide2 import QtGui
#
#from ui.MainWindowGUI import Ui_MainWindow

import configparser
import os

import numpy as np

import utils.utilities as ut

from managers_tui.input_mgr import InputManager
from managers_tui.config_mgr import ConfigManager
from managers_tui.database_mgr import DatabaseManager
from managers_tui.help_mgr import HelpManager


class BSG_UI():

    def __init__(self):
        self.inputMgr = InputManager()
        self.configMgr = ConfigManager()
        self.databaseMgr = DatabaseManager()
        self.helpMgr = HelpManager()

        self.log('Initialized...')


    def log(self, message):
        print(message + '\n')

    def loadExistingDecay(self,filename=None):
        self.inputMgr.loadIniFile(filename)
    
    def newDecay(self,filename=None):
        self.inputMgr.checkUnsavedTransitionChanges()
        
        isotope = input('Enter the isotope name (e.g. 238U): ')
        import re
        A = int(re.findall(r'\d+', isotope)[0])
        name = isotope.replace(str(A), '').capitalize()
        try:
            Z = ut.atoms.index(name)
        except ValueError as e:
            self.log('ValueError({0}): {1}'.format(e.errno, e.strerror))
            self.log('Isotope not in the list')
        
        answer = input('Search for transitions in ENSDF database (yes/no)? ')
        
        if answer == 'yes':
            branch = self.databaseMgr.loadESNDFBranches(Z, A)
            if branch == None:
                self.log('Transition not found in the ENSDF database')
                return
        else:
            process = input('Choose beta process(B-, B+ or EC) ')
            while True:
                str_temp = input('What is the end point energy of the transition (in keV) (Default is 1000. keV)? ')
                if str_temp == '':
                    end_point = 1000.
                    break
                else:
                    try:
                        end_point = float(str_temp)
                        break
                    except ValueError as e:
                        print('Input is not a float, try again!')
            end_point = input('What is the end point energy of the transition (in keV)? ')
            while True:
                str_temp = input('What is the logft value of the transition?(Default is 6.) ')
                if str_temp == '':
                    logft = 6.
                    break
                else:
                    try:
                        logft = float(str_temp)
                        break
                    except ValueError as e:
                        print('Input is not a float, try again!')
            while True:
                str_temp = input('What is the half-life of the transition?(Default is 1. s) ')
                if str_temp == '':
                    halflife = 1.
                    break
                else:
                    try:
                        halflife = float(str_temp)
                        break
                    except ValueError as e:
                        print('Input is not a float, try again!')

            branch = ut.BetaBranch(process,end_point,0,100.,0.,logft,Level(0., 0., '0+'),Level(0., 0., '0+'),halflife)
        self.inputMgr.setCurrentTransition(A,Z,branch)
#            self.ui.cb_process.setCurrentIndex(self.ui.cb_process.findText(process))
#            self.ui.sb_AD.setValue(A)
#            self.ui.sb_ZD.setValue(Z+1 if process == 'B-' else Z-1)

        answer = input('Enter deformation info from the Moeller database? (yes/no) ')

        if answer == 'yes':
            mDef, dDef = self.databaseMgr.loadDeformation()
            self.inputMgr.setCurrentDeformations(mDef,dDef)
        
        answer = input('Load charge radii from database? (yes/no) ')

        if answer == 'yes':
            mRad, dRad = self.databaseMgr.loadRadii()
            self.inputMgr.setCurrentRadii(mRad,dRad)
        
        answer = input('Save transition in ini file? (yes/no) ')
            
        if answer == 'yes':
            self.inputMgr.writeIniFile(filename)

    def runBSG(self):
#        self.inputMgr.checkUnsavedTransitionChanges()

        outputName = input("Enter a name for the output files: ")
        if self.databaseMgr.execPath == '':
            print("Set the path for the generator executable.")
        if self.inputMgr.iniName == '':
            print("No ini file found.")
            return

        self.log("Performing BSG Calculation...")
        self.log("Calculating...")
        command = "{0} -i {1} -o {2}".format(self.databaseMgr.execPath, self.inputMgr.iniName, outputName)
        command += " -e {}".format(self.databaseMgr.exchangePath)

        command += " -c {}".format(self.configMgr.configName)

        print("Executing command: %s" % command)
        try:
            sh = Shell()
            sh.run(command)
            if sh.output(raw=True) != '':
                print(sh.output(raw=True))
            self.log("Spectrum calculation OK")
        except CommandError as e:
            self.log('CommandError({0}): {1}'.format(e.errno, e.strerror))
        self.log("Ready!")

def exit_error():
    print("Usage: python bsg_terminal_ui.py -i <.ini filename>(optional) -c <config filename>(optional)\n")
    sys.exit()

if __name__ == '__main__':
    expected_nArgs = [1,3,5]
    if len(sys.argv) not in expected_nArgs:
        exit_error()
    elif len(sys.argv) == 5:
        if sys.argv[1] != '-i' and sys.argv[1] != '-c':
            exit_error()
        elif sys.argv[1] == '-i' and sys.argv[3] != '-c':
            exit_error()
        elif sys.argv[1] == '-c' and sys.argv[3] != '-i':
            exit_error()
    elif len(sys.argv) == 3:
        if sys.argv[1] != '-i' and sys.argv[1] != '-c':
            exit_error()

    mw = BSG_UI()
    while True:
        arg_array = [None,None]
        if len(sys.argv) == 5:
            if sys.argv[1] == "-i" and sys.argv[3] == "-c" :
                arg_array[0] = sys.argv[2]
                arg_array[1] = sys.argv[4]
            else:
                arg_array[0] = sys.argv[4]
                arg_array[1] = sys.argv[2]
        elif len(sys.argv) == 3:
            if sys.argv[1] == "-i":
                arg_array[0] = sys.argv[2]
            else:
                arg_array[1] = sys.argv[2]

        answer = input("create new .ini file from databases (yes), or load from existing .ini file (no)? ")
        if answer == 'yes':
            mw.newDecay(arg_array[0])
        else:
            mw.loadExistingDecay(arg_array[0])
        answer = input("create new config file from scratch (yes), or load from existing config file (no)? ")
        if answer == 'yes':
            mw.configMgr.createNewConfigFile(arg_array[1])
        else:
            mw.configMgr.loadConfigFile(arg_array[1])

        mw.runBSG()

        answer = input("Would you like to Proceed (yes) or exit the program (no) ? ")
        if answer == 'no':
            break

    sys.exit()
