#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""

import os

class DatabaseManager:
    def __init__(self):
        self.deformationFile = None
        self.radiiFile = None
        self.ensdfFolder = None
        self.execPath = ''
        self.exchangePath = ''
        self.findDefaults()

    def findDefaults(self):
        print("Looking for defaults in environment variables...")
        bsg_exec = os.environ.get('BSGPATH')
        if bsg_exec:
            self.execPath = bsg_exec
            print("Found bsg_exec")
        exchangePath = os.environ.get('BSGEXCHANGEPATH')
        if exchangePath:
            self.exchangePath = exchangePath
            print("Found Exchange data")
        ensdf = os.environ.get('ENSDFDIR')
        if ensdf:
            self.ensdfFolder = ensdf
            print("Found ENSDF database")
        frdm = os.environ.get('FRDMPATH')
        if frdm:
            self.deformationFile = frdm
            print("Found FRDM")
        chargeRadii = os.environ.get('CHARGERADIIPATH')
        if chargeRadii:
            self.radiiFile = chargeRadii
            print("Found Charge Radii")


    def loadESNDFBranches(self, Z, A):
        if not self.ensdfFolder:
            print("ENSDFDIR environment variable not set")
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
            print('Available beta transitions: ')
            for i,e in enumerate(items):
                print(i+1,': ',e)
            while True:
                try:
                    choice = int(input('Choose a transition (1- {})'.format(len(items))))
                    if choice>0 and choice <= len(items):
                        branch = bbs[choice]
                        return branch
                    else:
                        print("Wrong index! Choose a number between (1 and {})".format(len(items)))
                except ValueError as e:
                    print('Input is not a float, try again!')
        else:
            print("No transitions found.")
            return None

    def loadDeformation(self):
        if not self.deformationFile:
            print("No deformation file was set")
            return
        mDef = ut.loadDeformation(self.mZ, self.mA, self.deformationFile)
        dDef = ut.loadDeformation(self.dZ, self.dA, self.deformationFile)

        return np.delete(mDef,1),np.delete(dDef,1)
#        self.bsg_ui.ui.dsb_Beta2M.setValue(mDef[0])
#        self.bsg_ui.ui.dsb_Beta4M.setValue(mDef[2])
#        self.bsg_ui.ui.dsb_Beta6M.setValue(mDef[3])
#        self.bsg_ui.ui.dsb_Beta2D.setValue(dDef[0])
#        self.bsg_ui.ui.dsb_Beta4D.setValue(dDef[2])
#        self.bsg_ui.ui.dsb_Beta6D.setValue(dDef[3])

    def loadRadii(self):
        if not self.radiiFile:
            print("No radii file was set")
            return
        mRad = ut.loadChargeRadius(self.mZ, self.mA, self.radiiFile)
        dRad = ut.loadChargeRadius(self.dZ, self.dA, self.radiiFile)

        return mRad,dRad
#        self.bsg_ui.ui.dsb_RM.setValue(mRad)
#        self.bsg_ui.ui.dsb_RD.setValue(dRad)


