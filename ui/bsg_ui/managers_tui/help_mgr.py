#! /usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 3 18:24 2019

@author: leendert
"""

class HelpManager:
    def __init__(self,):
        self.about()
        self.visitENSDF()
        self.visitFRDM2012()
        self.visitChargeRadii()
        self.submitFeedback()


    def visitENSDF(self):
        print('https://www.nndc.bnl.gov/ensarchivals/')
    
    def visitFRDM2012(self):
        print('https://www.sciencedirect.com/science/article/pii/S0092640X1600005X')
    
    def visitChargeRadii(self):
        print('https://journals.aps.org/prc/abstract/10.1103/PhysRevC.94.064315')

    def about(self):
        print("Beta spectrum Generator GUI v 1.0")
    
    def submitFeedback(self):
        print("Send emails to leendert.hayen@gmail.com")
