#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""

import configparser
import utils.utilities as ut

class InputManager:
    def __init__(self):
        self.iniName = ''
        self.unsavedTransitionChanges = list()
        self.robtdFile = ''
        self.process = ''
        self.type = 'Gamow-Teller'
        self.mixingRatio = 0.
        self.qValue = 0.
        self.mZ = 0
        self.mA = 0
        self.mRad = 0.
        self.mJpi = 0
        self.mBeta2 = 0.
        self.mBeta4 = 0.
        self.mBeta6 = 0.
        self.mE = 0.
        self.dZ = 0
        self.dA = 0
        self.dRad = 0.
        self.dJpi = 0
        self.dBeta2 = 0.
        self.dBeta4 = 0.
        self.dBeta6 = 0.
        self.dE = 0.
        self.partialHalfLife = 1.
        self.LogFt = 6.

    def loadIniFile(self, filename):
        if not filename:
            filename = input("Give the name of the .ini file you want to load: ")
            if filename == '':
                return
        self.iniName = filename
#        self.checkUnsavedTransitionChanges()
        try:
            config = configparser.ConfigParser()
            config.read(self.iniName)
            print("File exists")
            self.process = config.get('Transition', 'Process')
            self.type = config.get('Transition', 'Type')
            self.mixingRatio = config.getfloat('Transition', 'MixingRatio')
            self.qValue = config.getfloat('Transition', 'QValue')
            self.mZ = config.getint('Mother', 'Z')
            self.mA = config.getint('Mother', 'A')
            self.mRad = config.getfloat('Mother', 'Radius')
            self.mJpi = config.get('Mother', 'SpinParity')
            self.mBeta2 = config.getfloat('Mother', 'Beta2')
            self.mBeta4 = config.getfloat('Mother', 'Beta4')
            self.mBeta6 = config.getfloat('Mother', 'Beta6')
            self.mE = config.getfloat('Mother', 'ExcitationEnergy')
            self.dZ = config.getint('Daughter', 'Z')
            self.dA = config.getint('Daughter', 'A')
            self.dRad = config.getfloat('Daughter', 'Radius')
            self.dJpi = config.get('Daughter', 'SpinParity')
            self.dBeta2 = config.getfloat('Daughter', 'Beta2')
            self.dBeta4 = config.getfloat('Daughter', 'Beta4')
            self.dBeta6 = config.getfloat('Daughter', 'Beta6')
            self.dE = config.getfloat('Daughter', 'ExcitationEnergy')
            try:
                self.partialHalfLife = config.getfloat('Transition', 'PartialHalflife')
            except configparser.NoOptionError:
                pass
            try:
                self.logFt = config.getfloat('Transition', 'LogFt')
            except configparser.NoOptionError:
                pass
            try:
                self.setROBTDFile(config.get('Transition', 'ROBTDFile'))
            except configparser.NoOptionError:
                pass
            print("Loaded Ini file: %s." % self.iniName)
        except:
            print("Failed to load %s file" % self.iniName)

    def writeIniFile(self,filename):
        phl = self.partialHalfLife if self.partialHalfLife > 0 else None
        
        logft = self.logFt if self.logFt > 0 else None
        
        robtdFile = self.robtdFile if not self.robtdFile == '' else None
        if not filename:
            filename = input("Give the name of the .ini file: ")
            if filename == '':
                return
        self.iniName = filename

        ut.writeIniFile(str(self.mZ), str(self.dZ), str(self.mA), str(self.qValue), str(self.process), str(self.type), str(self.mBeta2), str(self.mBeta4), str(self.mBeta6),str(self.dBeta2), str(self.dBeta4), str(self.dBeta6), str(self.mRad), str(self.dRad), str(self.mJpi), str(self.dJpi), mE=str(self.mE), dE=str(self.dE), name=filename,phl=str(phl), logft=str(logft), robtdFile=robtdFile)
        print('Ini file written to %s.' % filename)

    def setCurrentTransition(self,A,Z,branch):
        self.mA = A
        self.mZ = Z
        self.process = branch.process
        self.dZ = Z+1 if branch.process == 'B-' else Z-1
        self.dA = A
        self.mE = branch.motherLevel.E
        self.dE = branch.daughterLevel.E
        self.qValue = branch.E-branch.motherLevel.E+branch.daughterLevel.E
        self.mJpi = int(branch.motherLevel.Jpi.split('/2')[0]) if '/2' in branch.motherLevel.Jpi else int(branch.motherLevel.Jpi[:-1])*2
        self.dJpi = int(branch.daughterLevel.Jpi.split('/2')[0]) if '/2' in branch.daughterLevel.Jpi else int(branch.daughterLevel.Jpi[:-1])*2
        self.partialHalfLife = branch.partialHalflife
        self.logFt = branch.logft
    
    def setCurrentDeformations(self,mDef,dDef):
        self.mBeta2 = mDef[0]
        self.mBeta4 = mDef[1]
        self.mBeta6 = mDef[2]
        self.dBeta2 = dDef[0]
        self.dBeta4 = dDef[1]
        self.dBeta6 = dDef[2]
    
    def setCurrentRadii(self,mRad,dRad):
        self.dRad = dRad
        self.mRad = mRad

    def addUnsavedChanges(self, component):
        self.unsavedTransitionChanges.append(component)

    def clearUnsavedChanges(self):
        self.unsavedTransitionChanges = list()

    def askSaveChanges(self):
        if len(self.unsavedTransitionChanges):
            s = 'There are unsaved changes in the transition information:\n\n'
            for i in self.unsavedTransitionChanges:
                s += i + '\n'
            s += '\nThese are not yet taken into account. Save info? (yes/no)'
            answer = input(s)
            if answer == 'yes':
                self.writeIniFile()

    def checkUnsavedTransitionChanges(self):
        self.clearUnsavedChanges()
        if self.iniName == '':
            self.addUnsavedChanges("No information currently saved.")
        else:
            try:
                config = configparser.ConfigParser()
                config.read(self.iniName)
            except Exception as e:
                print(e)
                print("{} not found. Exciting!".format(self.iniName))
                return
            if config.get('Transition', 'Process') != self.process:
                self.addUnsavedChanges('Process (Prev. %s)' % config.get('Transition', 'Process'))
            if config.get('Transition', 'Type') != self.type:
                self.addUnsavedChanges('Type (Prev. %s)' % config.get('Transition', 'Type'))
            if config.getfloat('Transition', 'MixingRatio') != self.mixingRatio:
                self.addUnsavedChanges('Mixing Ratio (Prev. %s)' % config.get('Transition', 'MixingRatio'))
            if config.getfloat('Transition', 'QValue') != self.qValue:
                self.addUnsavedChanges('Q Value (Prev. %s)' % config.get('Transition', 'QValue'))
            try:
                if config.getfloat('Transition', 'PartialHalflife') != self.partialHalfLife:
                    self.addUnsavedChanges('Partial Halflife (Prev. %s)' % config.get('Transition', 'PartialHalflife'))
            except configparser.NoOptionError:
                pass
            try:
                if config.getfloat('Transition', 'LogFt') != self.logFt:
                    self.addUnsavedChanges('Log Ft (Prev. %s)' % config.get('Transition', 'LogFt'))
            except configparser.NoOptionError:
                pass

            if config.getint('Mother', 'Z') != self.mZ:
                self.addUnsavedChanges('Mother Z (Prev. %d)' % config.getint('Mother', 'Z'))
            if config.getint('Mother', 'A') != self.mA:
                self.addUnsavedChanges('Mother A (Prev. %d)' % config.getint('Mother', 'A'))
            if config.getfloat('Mother', 'Radius') != self.mRad:
                self.addUnsavedChanges('Mother Radius (Prev. %.2f)' % config.getfloat('Mother', 'Radius'))
            if config.get('Mother', 'SpinParity') != self.mJpi:
                self.addUnsavedChanges('Mother Spin (Prev. %.1f)' % abs(config.getint('Mother', 'SpinParity')/2.))
            if None:
                pass
            if config.getfloat('Mother', 'Beta2') != self.mBeta2:
                self.addUnsavedChanges('Mother beta2 (Prev. %.3f)' % config.getfloat('Mother', 'Beta2'))
            if config.getfloat('Mother', 'Beta4') != self.mBeta4:
                self.addUnsavedChanges('Mother beta4 (Prev. %.3f)' % config.getfloat('Mother', 'Beta4'))
            if config.getfloat('Mother', 'Beta6') != self.mBeta6:
                self.addUnsavedChanges('Mother beta6 (Prev. %.3f)' % config.getfloat('Mother', 'Beta6'))
            
            if config.getint('Daughter', 'Z') != self.dZ:
                self.addUnsavedChanges('Daughter Z (Prev. %d)' % config.getint('Daughter', 'Z'))
            if config.getint('Daughter', 'A') != self.dA:
                self.addUnsavedChanges('Daughter A (Prev. %d)' % config.getint('Daughter', 'A'))
            if config.getfloat('Daughter', 'Radius') != self.dRad:
                self.addUnsavedChanges('Daughter Radius (Prev. %.2f)' % config.getfloat('Daughter', 'Radius'))
            if config.get('Daughter', 'SpinParity') != self.dJpi:
                self.addUnsavedChanges('Daughter Spin (Prev. %.1f)' % abs(config.getint('Daughter', 'SpinParity')/2.))
            if None:
                pass
            if config.getfloat('Daughter', 'Beta2') != self.dBeta2:
                self.addUnsavedChanges('Daughter beta2 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta2'))
            if config.getfloat('Daughter', 'Beta4') != self.dBeta4:
                self.addUnsavedChanges('Daughter beta4 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta4'))
            if config.getfloat('Daughter', 'Beta6') != self.dBeta6:
                self.addUnsavedChanges('Daughter beta6 (Prev. %.3f)' % config.getfloat('Daughter', 'Beta6'))

        self.askSaveChanges()

    def setROBTDFile(self, filename=''):
        self.robtdFile = filename
