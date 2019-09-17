#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""

import configparser

class ConfigManager:
    def __init__(self,configName):
        self.configName = configName


        self.spectrumBoolOptions =  {"Phasespace": True, "Fermi": True, "Radiative": True,\
        "ESFiniteSize": True, "ESShape": True, "ESDeformation": True,\
        "NSShape": True, "Isovector": True, "Connect": False, "CDeformation": True,\
        "Relativistic": True, "Recoil": True, "CoulombRecoil": True, \
        "Screening": True, "Exchange": True, "AtomicMismatch": True, "EnableNrOfBins": False, "EnforceNME": False}
        
        self.spectrumStringOptions = {"ESShapeForm": "Fermi", "NSShapeForm": "UCS"}
        
        self.spectrumNME = {"WeakMagnetism": 0.0, "InducedTensor": 0.0, "Lambda": 0.0}
        
        self.spectrumDSB = {"Begin": 0.0, "End": 0.0,"StepSize": 0.10, "Steps": 0}
        
        self.spectrum = self.spectrumBoolOptions
        self.spectrum.update(self.spectrumStringOptions)
        self.spectrum.update(self.spectrumNME)
        self.spectrum.update(self.spectrumDSB)
        
        self.constants = {"gA": 1.2723, "gAeff": 1.00, "gP": 0.0,"gM": 4.706}
        
        self.computationalStringOptions = {"Method": "ROBTD", "Potential": "SHO"}
        
        self.computationalBoolOptions = {"ForceSpin": True,"ReversedGhallagher": False, "OverrideSPCoupling": False}
        
        self.computationalDSB = {"EnergyMargin": 1.0, "Vneutron": 49.60, "Vproton": 49.60, "V0Sneutron": 7.20,
            "V0Sproton": 7.20, "Xneutron": 0.86, "Xproton": 0.86}
        
        self.computational = self.computationalStringOptions
        self.computational.update(self.computationalBoolOptions)
        self.computational.update(self.computationalDSB)
        
        self.general = {"Folder": '.'}

        self.config_settings = {"General": self.general, "Spectrum": self.spectrum, "Constants": self.constants, "Computational": self.computational}


    def loadConfigFile(self):
        if self.configName == '':
            return

        config = configparser.ConfigParser()
        config.read(self.configName)

        
        for conf_key in self.config_settings:
            for key in self.config_settings[conf_key]:
                try:
                    if isinstance(self.config_settings[conf_key][key],bool):
                        self.config_settings[conf_key][key] = config.getboolean(conf_key, key)
                    elif isinstance(self.config_settings[conf_key][key],int):
                        self.config_settings[conf_key][key] = config.getint(conf_key, key)
                    elif isinstance(self.config_settings[conf_key][key],float):
                        self.config_settings[conf_key][key] = config.getfloat(conf_key, key)
                    else:
                        self.config_settings[conf_key][key] = config.get(conf_key, key)
                except (configparser.NoSectionError, configparser.NoOptionError):
                    continue

        print("Loaded config file: %s." % self.configName)

    def writeConfigFile(self):
        if self.configName == '':
            return
        try:
            ut.writeConfigFile(self.configName, self.config_settings)
            print("Config file written to %s." % self.configName)
        except:
            print("Writing config file failed.")

    def createNewConfigFile(self):
        print("Spectrum Boolean options: answer with yes or no to decide which spectrum options should be turned on")
        for key in self.spectrumBoolOptions:
            answer = input(self.spectrumBoolOptions[key]+" (yes/no): ")
            if answer == 'yes':
                self.spectrumBoolOptions[key] = True
            else:
                self.spectrumBoolOptions[key] = False

        print("Spectrum string options: answer with one of the options (numerical)")
        for key in self.spectrumStringOptions:
            if key == 'ESShapeForm':
                while True:
                    try:
                        answer = int(input("What shape should be used for the electrostatic potential: (1) Fermi, or (2) Modified Gauss ? "))
                        if answer < 1 or answer > 2:
                            print("choose between '1' or '2'")
                            continue
                        elif answer == 1:
                            self.spectrumStringOptions[key] = 'Fermi'
                            break
                        else:
                            self.spectrumStringOptions[key] = 'Mod. Gauss'
                            break
                    except ValueError as e:
                        print('Input is not a float, try again!')
            elif key == 'NSShapeForm':
                while True:
                    try:
                        answer = int(input("What shape should be used for the nuclear potential: (1) UCS, or (2) Modified Gauss ? "))
                        if answer < 1 or answer > 2:
                            print("choose between '1' or '2'")
                            continue
                        elif answer == 1:
                            self.spectrumStringOptions[key] = 'Fermi'
                            break
                        else:
                            self.spectrumStringOptions[key] = 'Mod. Gauss'
                            break
                    except ValueError as e:
                        print('Input is not a number, try again!')

        print("Spectrum nuclear matrix elements: answer with float")
        for key in self.spectrumNME:
            self.spectrumNME[key] = float(input(key+ " : "))

        print("Spectrum binning optione: starting and final energy binsize (in keV) and number of bins in case binsize is set to a negative value")
        for key in self.spectrumDSB:
            self.spectrumDSB[key] = float(input(key+ " : "))

        print("Set the constants")
        for key in self.constants:
            self.constants[key] = float(input(key + " : "))

        print("Boolean options for some computational methods, answer with yes or no")
        for key in self.computationalBoolOptions:
            answer = input(self.computationalBoolOptions[key]+" (yes/no): ")
            if answer == 'yes':
                self.computationalBoolOptions[key] = True
            else:
                self.computationalBoolOptions[key] = False

        print("Spectrum string options: answer with one of the options (numerical)")
        while True:
            try:
                answer = int(input("Which method should be used for calculation of the nuclear matrix elements (1) ESP, or (2) ROBTD ? "))
                if answer < 1 or answer > 2:
                    print("choose between '1' or '2'")
                    continue
                elif answer == 1:
                    self.computationalStringOptions['Method'] = 'ESP'
                    break
                else:
                    self.computationalStringOptions['Method'] = 'ROBTD'
                    break
            except ValueError as e:
                print('Input is not a number, try again!')
        if self.computationalStringOptions['Method'] == 'ESP':
            while True:
                try:
                    answer = int(input("What shape should be used for the nuclear potential in the ESP model: (1) SHO, (2) WS or (3) DWS ? "))
                    if answer < 1 or answer > 3:
                        print("choose between '1' or '2'")
                        continue
                    elif answer == 1:
                        self.computationalStringOptions['Potential'] = 'SHO'
                        break
                    elif answer == 2:
                        self.computationalStringOptions['Potential'] = 'WS'
                        break
                    else:
                        self.computationalStringOptions['Potential'] = 'DWS'
                        break
                except ValueError as e:
                    print('Input is not a float, try again!')

            print("Set the parameters for the ESP model (float)")
            for key in self.computationalDSB:
                self.computationalDSB[key] = float(input(key + " : "))

        self.writeConfigFile()
