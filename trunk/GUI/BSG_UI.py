# -*- coding: utf-8 -*-
"""
Created on Fri Sep 19 14:08:39 2014

@author: leendert
"""
import sys

from PySide import QtCore, QtGui, QtUiTools

class BSG_UI(object):
    
    def __init__(self):
        loader = QtUiTools.QUiLoader(); 
        gfile = QtCore.QFile("MainWindow.ui")
        gfile.open(QtCore.QFile.ReadOnly)
        self.ui = loader.load(gfile)
        gfile.close()


        self.ui.show()
        
    
if __name__ == '__main__':
    app = QtGui.QApplication(sys.argv)
    gon = BSG_UI()
    app.exec_()
    
    sys.exit()