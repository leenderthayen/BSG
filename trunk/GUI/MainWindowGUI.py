# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'MainWindow.ui'
#
# Created: Thu Dec  7 16:04:17 2017
#      by: pyside-uic 0.2.15 running on PySide 1.2.4
#
# WARNING! All changes made in this file will be lost!

from PySide import QtCore, QtGui

class Ui_MainWindow(object):
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(912, 748)
        self.centralwidget = QtGui.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        self.labelStatus = QtGui.QLabel(self.centralwidget)
        self.labelStatus.setGeometry(QtCore.QRect(20, 730, 371, 16))
        self.labelStatus.setObjectName("labelStatus")
        self.scrollArea = QtGui.QScrollArea(self.centralwidget)
        self.scrollArea.setGeometry(QtCore.QRect(520, 500, 361, 131))
        self.scrollArea.setWidgetResizable(True)
        self.scrollArea.setObjectName("scrollArea")
        self.scrollAreaWidgetContents = QtGui.QWidget()
        self.scrollAreaWidgetContents.setGeometry(QtCore.QRect(0, 0, 359, 129))
        self.scrollAreaWidgetContents.setObjectName("scrollAreaWidgetContents")
        self.lastActions = QtGui.QPlainTextEdit(self.scrollAreaWidgetContents)
        self.lastActions.setGeometry(QtCore.QRect(0, 0, 361, 131))
        self.lastActions.setObjectName("lastActions")
        self.scrollArea.setWidget(self.scrollAreaWidgetContents)
        self.optionsTabWidget = QtGui.QTabWidget(self.centralwidget)
        self.optionsTabWidget.setGeometry(QtCore.QRect(510, 10, 381, 451))
        self.optionsTabWidget.setObjectName("optionsTabWidget")
        self.iniTab = QtGui.QWidget()
        self.iniTab.setObjectName("iniTab")
        self.groupBox = QtGui.QGroupBox(self.iniTab)
        self.groupBox.setGeometry(QtCore.QRect(180, 0, 191, 151))
        self.groupBox.setObjectName("groupBox")
        self.label = QtGui.QLabel(self.groupBox)
        self.label.setGeometry(QtCore.QRect(0, 30, 67, 21))
        self.label.setObjectName("label")
        self.label_2 = QtGui.QLabel(self.groupBox)
        self.label_2.setGeometry(QtCore.QRect(0, 60, 67, 21))
        self.label_2.setObjectName("label_2")
        self.label_3 = QtGui.QLabel(self.groupBox)
        self.label_3.setGeometry(QtCore.QRect(0, 90, 101, 21))
        self.label_3.setObjectName("label_3")
        self.label_4 = QtGui.QLabel(self.groupBox)
        self.label_4.setGeometry(QtCore.QRect(0, 120, 67, 21))
        self.label_4.setObjectName("label_4")
        self.cb_type = QtGui.QComboBox(self.groupBox)
        self.cb_type.setGeometry(QtCore.QRect(110, 60, 81, 27))
        self.cb_type.setObjectName("cb_type")
        self.cb_process = QtGui.QComboBox(self.groupBox)
        self.cb_process.setGeometry(QtCore.QRect(110, 30, 81, 27))
        self.cb_process.setObjectName("cb_process")
        self.dsb_MixingRatio = QtGui.QDoubleSpinBox(self.groupBox)
        self.dsb_MixingRatio.setGeometry(QtCore.QRect(100, 90, 91, 27))
        self.dsb_MixingRatio.setObjectName("dsb_MixingRatio")
        self.dsb_Q = QtGui.QDoubleSpinBox(self.groupBox)
        self.dsb_Q.setGeometry(QtCore.QRect(100, 120, 91, 27))
        self.dsb_Q.setObjectName("dsb_Q")
        self.group_M = QtGui.QGroupBox(self.iniTab)
        self.group_M.setGeometry(QtCore.QRect(10, 150, 171, 261))
        self.group_M.setObjectName("group_M")
        self.lab_ZM = QtGui.QLabel(self.group_M)
        self.lab_ZM.setGeometry(QtCore.QRect(0, 20, 61, 31))
        self.lab_ZM.setObjectName("lab_ZM")
        self.lab_AM = QtGui.QLabel(self.group_M)
        self.lab_AM.setGeometry(QtCore.QRect(0, 50, 67, 31))
        self.lab_AM.setObjectName("lab_AM")
        self.lab_RM = QtGui.QLabel(self.group_M)
        self.lab_RM.setGeometry(QtCore.QRect(0, 80, 67, 31))
        self.lab_RM.setObjectName("lab_RM")
        self.lab_JM = QtGui.QLabel(self.group_M)
        self.lab_JM.setGeometry(QtCore.QRect(0, 110, 67, 31))
        self.lab_JM.setObjectName("lab_JM")
        self.lab_Beta2M = QtGui.QLabel(self.group_M)
        self.lab_Beta2M.setGeometry(QtCore.QRect(0, 170, 67, 31))
        self.lab_Beta2M.setObjectName("lab_Beta2M")
        self.lab_Beta4M = QtGui.QLabel(self.group_M)
        self.lab_Beta4M.setGeometry(QtCore.QRect(0, 200, 67, 31))
        self.lab_Beta4M.setObjectName("lab_Beta4M")
        self.lab_Beta6M = QtGui.QLabel(self.group_M)
        self.lab_Beta6M.setGeometry(QtCore.QRect(0, 226, 67, 31))
        self.lab_Beta6M.setObjectName("lab_Beta6M")
        self.lab_PiM = QtGui.QLabel(self.group_M)
        self.lab_PiM.setGeometry(QtCore.QRect(0, 140, 67, 31))
        self.lab_PiM.setObjectName("lab_PiM")
        self.sb_ZM = QtGui.QSpinBox(self.group_M)
        self.sb_ZM.setGeometry(QtCore.QRect(110, 20, 48, 27))
        self.sb_ZM.setObjectName("sb_ZM")
        self.sb_AM = QtGui.QSpinBox(self.group_M)
        self.sb_AM.setGeometry(QtCore.QRect(110, 50, 48, 27))
        self.sb_AM.setObjectName("sb_AM")
        self.dsb_RM = QtGui.QDoubleSpinBox(self.group_M)
        self.dsb_RM.setGeometry(QtCore.QRect(80, 80, 81, 27))
        self.dsb_RM.setSingleStep(0.1)
        self.dsb_RM.setObjectName("dsb_RM")
        self.dsb_JM = QtGui.QDoubleSpinBox(self.group_M)
        self.dsb_JM.setGeometry(QtCore.QRect(90, 110, 71, 27))
        self.dsb_JM.setDecimals(1)
        self.dsb_JM.setSingleStep(0.5)
        self.dsb_JM.setObjectName("dsb_JM")
        self.cb_PiM = QtGui.QComboBox(self.group_M)
        self.cb_PiM.setGeometry(QtCore.QRect(100, 140, 61, 27))
        self.cb_PiM.setObjectName("cb_PiM")
        self.dsb_Beta2M = QtGui.QDoubleSpinBox(self.group_M)
        self.dsb_Beta2M.setGeometry(QtCore.QRect(90, 170, 69, 27))
        self.dsb_Beta2M.setDecimals(3)
        self.dsb_Beta2M.setMinimum(-1.0)
        self.dsb_Beta2M.setMaximum(1.0)
        self.dsb_Beta2M.setSingleStep(0.1)
        self.dsb_Beta2M.setObjectName("dsb_Beta2M")
        self.dsb_Beta4M = QtGui.QDoubleSpinBox(self.group_M)
        self.dsb_Beta4M.setGeometry(QtCore.QRect(90, 200, 69, 27))
        self.dsb_Beta4M.setDecimals(3)
        self.dsb_Beta4M.setMinimum(-1.0)
        self.dsb_Beta4M.setMaximum(1.0)
        self.dsb_Beta4M.setSingleStep(0.1)
        self.dsb_Beta4M.setObjectName("dsb_Beta4M")
        self.dsb_Beta6M = QtGui.QDoubleSpinBox(self.group_M)
        self.dsb_Beta6M.setGeometry(QtCore.QRect(90, 230, 69, 27))
        self.dsb_Beta6M.setDecimals(3)
        self.dsb_Beta6M.setMinimum(-1.0)
        self.dsb_Beta6M.setMaximum(1.0)
        self.dsb_Beta6M.setSingleStep(0.1)
        self.dsb_Beta6M.setObjectName("dsb_Beta6M")
        self.groupBox_3 = QtGui.QGroupBox(self.iniTab)
        self.groupBox_3.setGeometry(QtCore.QRect(10, 0, 161, 141))
        self.groupBox_3.setObjectName("groupBox_3")
        self.label_11 = QtGui.QLabel(self.groupBox_3)
        self.label_11.setGeometry(QtCore.QRect(0, 20, 91, 21))
        self.label_11.setObjectName("label_11")
        self.l_iniFilename = QtGui.QLabel(self.groupBox_3)
        self.l_iniFilename.setGeometry(QtCore.QRect(0, 50, 161, 17))
        self.l_iniFilename.setObjectName("l_iniFilename")
        self.b_load_ini = QtGui.QPushButton(self.groupBox_3)
        self.b_load_ini.setGeometry(QtCore.QRect(10, 80, 131, 27))
        self.b_load_ini.setObjectName("b_load_ini")
        self.b_save_ini = QtGui.QPushButton(self.groupBox_3)
        self.b_save_ini.setGeometry(QtCore.QRect(10, 110, 131, 27))
        self.b_save_ini.setObjectName("b_save_ini")
        self.group_M_2 = QtGui.QGroupBox(self.iniTab)
        self.group_M_2.setGeometry(QtCore.QRect(180, 150, 191, 261))
        self.group_M_2.setObjectName("group_M_2")
        self.lab_ZD = QtGui.QLabel(self.group_M_2)
        self.lab_ZD.setGeometry(QtCore.QRect(0, 20, 61, 31))
        self.lab_ZD.setObjectName("lab_ZD")
        self.lab_AD = QtGui.QLabel(self.group_M_2)
        self.lab_AD.setGeometry(QtCore.QRect(0, 50, 67, 31))
        self.lab_AD.setObjectName("lab_AD")
        self.lab_RD = QtGui.QLabel(self.group_M_2)
        self.lab_RD.setGeometry(QtCore.QRect(0, 80, 67, 31))
        self.lab_RD.setObjectName("lab_RD")
        self.lab_JD = QtGui.QLabel(self.group_M_2)
        self.lab_JD.setGeometry(QtCore.QRect(0, 110, 67, 31))
        self.lab_JD.setObjectName("lab_JD")
        self.lab_Beta2D = QtGui.QLabel(self.group_M_2)
        self.lab_Beta2D.setGeometry(QtCore.QRect(0, 170, 67, 31))
        self.lab_Beta2D.setObjectName("lab_Beta2D")
        self.lab_Beta4D = QtGui.QLabel(self.group_M_2)
        self.lab_Beta4D.setGeometry(QtCore.QRect(0, 200, 67, 31))
        self.lab_Beta4D.setObjectName("lab_Beta4D")
        self.lab_Beta6D = QtGui.QLabel(self.group_M_2)
        self.lab_Beta6D.setGeometry(QtCore.QRect(0, 226, 67, 31))
        self.lab_Beta6D.setObjectName("lab_Beta6D")
        self.lab_PiD = QtGui.QLabel(self.group_M_2)
        self.lab_PiD.setGeometry(QtCore.QRect(0, 140, 67, 31))
        self.lab_PiD.setObjectName("lab_PiD")
        self.sb_ZD = QtGui.QSpinBox(self.group_M_2)
        self.sb_ZD.setGeometry(QtCore.QRect(140, 20, 48, 27))
        self.sb_ZD.setObjectName("sb_ZD")
        self.sb_AD = QtGui.QSpinBox(self.group_M_2)
        self.sb_AD.setGeometry(QtCore.QRect(140, 50, 48, 27))
        self.sb_AD.setObjectName("sb_AD")
        self.dsb_RD = QtGui.QDoubleSpinBox(self.group_M_2)
        self.dsb_RD.setGeometry(QtCore.QRect(110, 80, 81, 27))
        self.dsb_RD.setSingleStep(0.1)
        self.dsb_RD.setObjectName("dsb_RD")
        self.dsb_JD = QtGui.QDoubleSpinBox(self.group_M_2)
        self.dsb_JD.setGeometry(QtCore.QRect(120, 110, 71, 27))
        self.dsb_JD.setDecimals(1)
        self.dsb_JD.setSingleStep(0.5)
        self.dsb_JD.setObjectName("dsb_JD")
        self.cb_PiD = QtGui.QComboBox(self.group_M_2)
        self.cb_PiD.setGeometry(QtCore.QRect(130, 140, 61, 27))
        self.cb_PiD.setObjectName("cb_PiD")
        self.dsb_Beta2D = QtGui.QDoubleSpinBox(self.group_M_2)
        self.dsb_Beta2D.setGeometry(QtCore.QRect(120, 170, 69, 27))
        self.dsb_Beta2D.setDecimals(3)
        self.dsb_Beta2D.setMinimum(-1.0)
        self.dsb_Beta2D.setMaximum(1.0)
        self.dsb_Beta2D.setSingleStep(0.1)
        self.dsb_Beta2D.setObjectName("dsb_Beta2D")
        self.dsb_Beta4D = QtGui.QDoubleSpinBox(self.group_M_2)
        self.dsb_Beta4D.setGeometry(QtCore.QRect(120, 200, 69, 27))
        self.dsb_Beta4D.setDecimals(3)
        self.dsb_Beta4D.setMinimum(-1.0)
        self.dsb_Beta4D.setMaximum(1.0)
        self.dsb_Beta4D.setSingleStep(0.1)
        self.dsb_Beta4D.setObjectName("dsb_Beta4D")
        self.dsb_Beta6D = QtGui.QDoubleSpinBox(self.group_M_2)
        self.dsb_Beta6D.setGeometry(QtCore.QRect(120, 230, 69, 27))
        self.dsb_Beta6D.setDecimals(3)
        self.dsb_Beta6D.setMinimum(-1.0)
        self.dsb_Beta6D.setMaximum(1.0)
        self.dsb_Beta6D.setSingleStep(0.1)
        self.dsb_Beta6D.setObjectName("dsb_Beta6D")
        self.optionsTabWidget.addTab(self.iniTab, "")
        self.configTab = QtGui.QWidget()
        self.configTab.setObjectName("configTab")
        self.optionsTabWidget.addTab(self.configTab, "")
        self.line = QtGui.QFrame(self.centralwidget)
        self.line.setGeometry(QtCore.QRect(520, 470, 361, 20))
        self.line.setFrameShape(QtGui.QFrame.HLine)
        self.line.setFrameShadow(QtGui.QFrame.Sunken)
        self.line.setObjectName("line")
        self.lab_logo = QtGui.QLabel(self.centralwidget)
        self.lab_logo.setGeometry(QtCore.QRect(10, 10, 481, 81))
        self.lab_logo.setObjectName("lab_logo")
        self.mainTabWidget = QtGui.QTabWidget(self.centralwidget)
        self.mainTabWidget.setGeometry(QtCore.QRect(10, 120, 491, 611))
        self.mainTabWidget.setObjectName("mainTabWidget")
        self.bsgTab = QtGui.QWidget()
        self.bsgTab.setObjectName("bsgTab")
        self.toolBox = QtGui.QToolBox(self.bsgTab)
        self.toolBox.setGeometry(QtCore.QRect(10, 10, 461, 541))
        self.toolBox.setObjectName("toolBox")
        self.page = QtGui.QWidget()
        self.page.setGeometry(QtCore.QRect(0, 0, 461, 479))
        self.page.setObjectName("page")
        self.groupBox_2 = QtGui.QGroupBox(self.page)
        self.groupBox_2.setGeometry(QtCore.QRect(10, 0, 461, 491))
        self.groupBox_2.setTitle("")
        self.groupBox_2.setObjectName("groupBox_2")
        self.cb_phasespace = QtGui.QCheckBox(self.groupBox_2)
        self.cb_phasespace.setGeometry(QtCore.QRect(10, 30, 131, 22))
        self.cb_phasespace.setChecked(True)
        self.cb_phasespace.setObjectName("cb_phasespace")
        self.label_5 = QtGui.QLabel(self.groupBox_2)
        self.label_5.setGeometry(QtCore.QRect(10, 0, 67, 17))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setItalic(False)
        font.setUnderline(False)
        font.setBold(True)
        self.label_5.setFont(font)
        self.label_5.setObjectName("label_5")
        self.cb_fermi = QtGui.QCheckBox(self.groupBox_2)
        self.cb_fermi.setGeometry(QtCore.QRect(10, 60, 151, 22))
        self.cb_fermi.setChecked(True)
        self.cb_fermi.setObjectName("cb_fermi")
        self.label_6 = QtGui.QLabel(self.groupBox_2)
        self.label_6.setGeometry(QtCore.QRect(10, 120, 141, 17))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_6.setFont(font)
        self.label_6.setObjectName("label_6")
        self.cb_radiative = QtGui.QCheckBox(self.groupBox_2)
        self.cb_radiative.setGeometry(QtCore.QRect(10, 90, 181, 22))
        self.cb_radiative.setChecked(True)
        self.cb_radiative.setObjectName("cb_radiative")
        self.cb_efs = QtGui.QCheckBox(self.groupBox_2)
        self.cb_efs.setGeometry(QtCore.QRect(10, 150, 151, 22))
        self.cb_efs.setChecked(True)
        self.cb_efs.setObjectName("cb_efs")
        self.cb_ens = QtGui.QCheckBox(self.groupBox_2)
        self.cb_ens.setGeometry(QtCore.QRect(10, 180, 171, 22))
        self.cb_ens.setChecked(True)
        self.cb_ens.setObjectName("cb_ens")
        self.cb_efs_deformation = QtGui.QCheckBox(self.groupBox_2)
        self.cb_efs_deformation.setGeometry(QtCore.QRect(10, 240, 151, 31))
        self.cb_efs_deformation.setChecked(True)
        self.cb_efs_deformation.setObjectName("cb_efs_deformation")
        self.label_7 = QtGui.QLabel(self.groupBox_2)
        self.label_7.setGeometry(QtCore.QRect(230, 0, 111, 17))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_7.setFont(font)
        self.label_7.setObjectName("label_7")
        self.cb_recoil = QtGui.QCheckBox(self.groupBox_2)
        self.cb_recoil.setGeometry(QtCore.QRect(230, 30, 141, 22))
        self.cb_recoil.setChecked(True)
        self.cb_recoil.setObjectName("cb_recoil")
        self.cb_Q = QtGui.QCheckBox(self.groupBox_2)
        self.cb_Q.setGeometry(QtCore.QRect(230, 60, 141, 22))
        self.cb_Q.setChecked(True)
        self.cb_Q.setObjectName("cb_Q")
        self.label_8 = QtGui.QLabel(self.groupBox_2)
        self.label_8.setGeometry(QtCore.QRect(230, 90, 67, 17))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_8.setFont(font)
        self.label_8.setObjectName("label_8")
        self.cb_screening = QtGui.QCheckBox(self.groupBox_2)
        self.cb_screening.setGeometry(QtCore.QRect(230, 120, 171, 22))
        self.cb_screening.setChecked(True)
        self.cb_screening.setObjectName("cb_screening")
        self.cb_exchange = QtGui.QCheckBox(self.groupBox_2)
        self.cb_exchange.setGeometry(QtCore.QRect(230, 150, 171, 22))
        self.cb_exchange.setChecked(True)
        self.cb_exchange.setObjectName("cb_exchange")
        self.rb_overlap = QtGui.QRadioButton(self.groupBox_2)
        self.rb_overlap.setGeometry(QtCore.QRect(230, 180, 151, 22))
        self.rb_overlap.setChecked(True)
        self.rb_overlap.setObjectName("rb_overlap")
        self.rb_dW0 = QtGui.QRadioButton(self.groupBox_2)
        self.rb_dW0.setGeometry(QtCore.QRect(230, 210, 181, 22))
        self.rb_dW0.setChecked(False)
        self.rb_dW0.setObjectName("rb_dW0")
        self.dsb_dW0 = QtGui.QDoubleSpinBox(self.groupBox_2)
        self.dsb_dW0.setEnabled(False)
        self.dsb_dW0.setGeometry(QtCore.QRect(260, 240, 101, 27))
        self.dsb_dW0.setDecimals(3)
        self.dsb_dW0.setSingleStep(0.1)
        self.dsb_dW0.setObjectName("dsb_dW0")
        self.cb_shape = QtGui.QComboBox(self.groupBox_2)
        self.cb_shape.setGeometry(QtCore.QRect(30, 210, 121, 27))
        self.cb_shape.setObjectName("cb_shape")
        self.label_9 = QtGui.QLabel(self.groupBox_2)
        self.label_9.setGeometry(QtCore.QRect(10, 280, 141, 17))
        font = QtGui.QFont()
        font.setWeight(75)
        font.setBold(True)
        self.label_9.setFont(font)
        self.label_9.setObjectName("label_9")
        self.cb_C = QtGui.QCheckBox(self.groupBox_2)
        self.cb_C.setGeometry(QtCore.QRect(10, 300, 121, 31))
        self.cb_C.setChecked(True)
        self.cb_C.setObjectName("cb_C")
        self.cb_isovector = QtGui.QCheckBox(self.groupBox_2)
        self.cb_isovector.setGeometry(QtCore.QRect(10, 360, 181, 22))
        self.cb_isovector.setChecked(True)
        self.cb_isovector.setObjectName("cb_isovector")
        self.cb_coupled = QtGui.QCheckBox(self.groupBox_2)
        self.cb_coupled.setGeometry(QtCore.QRect(40, 380, 91, 31))
        self.cb_coupled.setObjectName("cb_coupled")
        self.cb_C_deformation = QtGui.QCheckBox(self.groupBox_2)
        self.cb_C_deformation.setGeometry(QtCore.QRect(10, 410, 121, 22))
        self.cb_C_deformation.setChecked(True)
        self.cb_C_deformation.setObjectName("cb_C_deformation")
        self.cb_C_shape = QtGui.QComboBox(self.groupBox_2)
        self.cb_C_shape.setGeometry(QtCore.QRect(30, 330, 121, 27))
        self.cb_C_shape.setObjectName("cb_C_shape")
        self.cb_relativistic = QtGui.QCheckBox(self.groupBox_2)
        self.cb_relativistic.setGeometry(QtCore.QRect(10, 440, 111, 22))
        self.cb_relativistic.setChecked(True)
        self.cb_relativistic.setObjectName("cb_relativistic")
        self.b_enable_all = QtGui.QPushButton(self.groupBox_2)
        self.b_enable_all.setGeometry(QtCore.QRect(240, 280, 201, 51))
        self.b_enable_all.setObjectName("b_enable_all")
        self.b_disable_all = QtGui.QPushButton(self.groupBox_2)
        self.b_disable_all.setGeometry(QtCore.QRect(240, 340, 201, 51))
        self.b_disable_all.setObjectName("b_disable_all")
        self.b_load_profile = QtGui.QPushButton(self.groupBox_2)
        self.b_load_profile.setGeometry(QtCore.QRect(240, 400, 99, 51))
        self.b_load_profile.setObjectName("b_load_profile")
        self.b_save_profile = QtGui.QPushButton(self.groupBox_2)
        self.b_save_profile.setGeometry(QtCore.QRect(340, 400, 99, 51))
        self.b_save_profile.setObjectName("b_save_profile")
        self.toolBox.addItem(self.page, "")
        self.page_2 = QtGui.QWidget()
        self.page_2.setGeometry(QtCore.QRect(0, 0, 461, 479))
        self.page_2.setObjectName("page_2")
        self.label_10 = QtGui.QLabel(self.page_2)
        self.label_10.setGeometry(QtCore.QRect(10, 0, 67, 17))
        self.label_10.setObjectName("label_10")
        self.toolBox.addItem(self.page_2, "")
        self.mainTabWidget.addTab(self.bsgTab, "")
        self.nmeTab = QtGui.QWidget()
        self.nmeTab.setObjectName("nmeTab")
        self.mainTabWidget.addTab(self.nmeTab, "")
        MainWindow.setCentralWidget(self.centralwidget)

        self.retranslateUi(MainWindow)
        self.optionsTabWidget.setCurrentIndex(0)
        self.mainTabWidget.setCurrentIndex(0)
        self.toolBox.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        MainWindow.setWindowTitle(QtGui.QApplication.translate("MainWindow", "Goniometer Controller", None, QtGui.QApplication.UnicodeUTF8))
        self.labelStatus.setText(QtGui.QApplication.translate("MainWindow", "Ready", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox.setTitle(QtGui.QApplication.translate("MainWindow", "Transition", None, QtGui.QApplication.UnicodeUTF8))
        self.label.setText(QtGui.QApplication.translate("MainWindow", "Process:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_2.setText(QtGui.QApplication.translate("MainWindow", "Type:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_3.setText(QtGui.QApplication.translate("MainWindow", "Mixing Ratio:", None, QtGui.QApplication.UnicodeUTF8))
        self.label_4.setText(QtGui.QApplication.translate("MainWindow", "Q value: ", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_Q.setSuffix(QtGui.QApplication.translate("MainWindow", " keV", None, QtGui.QApplication.UnicodeUTF8))
        self.group_M.setTitle(QtGui.QApplication.translate("MainWindow", "Mother", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_ZM.setText(QtGui.QApplication.translate("MainWindow", "Z:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_AM.setText(QtGui.QApplication.translate("MainWindow", "A:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_RM.setText(QtGui.QApplication.translate("MainWindow", "Radius:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_JM.setText(QtGui.QApplication.translate("MainWindow", "Spin:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_Beta2M.setText(QtGui.QApplication.translate("MainWindow", "Beta2:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_Beta4M.setText(QtGui.QApplication.translate("MainWindow", "Beta4:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_Beta6M.setText(QtGui.QApplication.translate("MainWindow", "Beta6:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_PiM.setText(QtGui.QApplication.translate("MainWindow", "Parity:", None, QtGui.QApplication.UnicodeUTF8))
        self.sb_ZM.setStatusTip(QtGui.QApplication.translate("MainWindow", "Atomic number", None, QtGui.QApplication.UnicodeUTF8))
        self.sb_AM.setStatusTip(QtGui.QApplication.translate("MainWindow", "Mass number", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_RM.setStatusTip(QtGui.QApplication.translate("MainWindow", "Charge radius", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_RM.setSuffix(QtGui.QApplication.translate("MainWindow", " fm", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_JM.setStatusTip(QtGui.QApplication.translate("MainWindow", "Spin in units of HBAR", None, QtGui.QApplication.UnicodeUTF8))
        self.groupBox_3.setTitle(QtGui.QApplication.translate("MainWindow", "File info", None, QtGui.QApplication.UnicodeUTF8))
        self.label_11.setText(QtGui.QApplication.translate("MainWindow", "Current file:", None, QtGui.QApplication.UnicodeUTF8))
        self.l_iniFilename.setText(QtGui.QApplication.translate("MainWindow", "Filename placeholder", None, QtGui.QApplication.UnicodeUTF8))
        self.b_load_ini.setText(QtGui.QApplication.translate("MainWindow", "Load .ini file", None, QtGui.QApplication.UnicodeUTF8))
        self.b_save_ini.setText(QtGui.QApplication.translate("MainWindow", "Save as...", None, QtGui.QApplication.UnicodeUTF8))
        self.group_M_2.setTitle(QtGui.QApplication.translate("MainWindow", "Daughter", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_ZD.setText(QtGui.QApplication.translate("MainWindow", "Z:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_AD.setText(QtGui.QApplication.translate("MainWindow", "A:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_RD.setText(QtGui.QApplication.translate("MainWindow", "Radius:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_JD.setText(QtGui.QApplication.translate("MainWindow", "Spin:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_Beta2D.setText(QtGui.QApplication.translate("MainWindow", "Beta2:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_Beta4D.setText(QtGui.QApplication.translate("MainWindow", "Beta4:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_Beta6D.setText(QtGui.QApplication.translate("MainWindow", "Beta6:", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_PiD.setText(QtGui.QApplication.translate("MainWindow", "Parity:", None, QtGui.QApplication.UnicodeUTF8))
        self.sb_ZD.setStatusTip(QtGui.QApplication.translate("MainWindow", "Atomic number", None, QtGui.QApplication.UnicodeUTF8))
        self.sb_AD.setStatusTip(QtGui.QApplication.translate("MainWindow", "Mass number", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_RD.setStatusTip(QtGui.QApplication.translate("MainWindow", "Charge radius", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_RD.setSuffix(QtGui.QApplication.translate("MainWindow", " fm", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_JD.setStatusTip(QtGui.QApplication.translate("MainWindow", "Spin in units of HBAR", None, QtGui.QApplication.UnicodeUTF8))
        self.optionsTabWidget.setTabText(self.optionsTabWidget.indexOf(self.iniTab), QtGui.QApplication.translate("MainWindow", "Ini file", None, QtGui.QApplication.UnicodeUTF8))
        self.optionsTabWidget.setTabText(self.optionsTabWidget.indexOf(self.configTab), QtGui.QApplication.translate("MainWindow", "Config file", None, QtGui.QApplication.UnicodeUTF8))
        self.lab_logo.setText(QtGui.QApplication.translate("MainWindow", "Logo", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_phasespace.setText(QtGui.QApplication.translate("MainWindow", "Phase Space", None, QtGui.QApplication.UnicodeUTF8))
        self.label_5.setText(QtGui.QApplication.translate("MainWindow", "Basic:", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_fermi.setText(QtGui.QApplication.translate("MainWindow", "Fermi Function", None, QtGui.QApplication.UnicodeUTF8))
        self.label_6.setText(QtGui.QApplication.translate("MainWindow", "Electrostatic:", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_radiative.setText(QtGui.QApplication.translate("MainWindow", "Radiative corrections", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_efs.setText(QtGui.QApplication.translate("MainWindow", "ES Finite Size", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_ens.setText(QtGui.QApplication.translate("MainWindow", "ES Nuclear Shape", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_efs_deformation.setText(QtGui.QApplication.translate("MainWindow", "ES Deformation", None, QtGui.QApplication.UnicodeUTF8))
        self.label_7.setText(QtGui.QApplication.translate("MainWindow", "Kinematic:", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_recoil.setText(QtGui.QApplication.translate("MainWindow", "Nuclear Recoil", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_Q.setText(QtGui.QApplication.translate("MainWindow", "Coulomb Recoil", None, QtGui.QApplication.UnicodeUTF8))
        self.label_8.setText(QtGui.QApplication.translate("MainWindow", "Atomic:", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_screening.setText(QtGui.QApplication.translate("MainWindow", "Screening", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_exchange.setText(QtGui.QApplication.translate("MainWindow", "Exchange", None, QtGui.QApplication.UnicodeUTF8))
        self.rb_overlap.setText(QtGui.QApplication.translate("MainWindow", "Atomic Overlap", None, QtGui.QApplication.UnicodeUTF8))
        self.rb_dW0.setText(QtGui.QApplication.translate("MainWindow", "Explicit $Delta W_0$", None, QtGui.QApplication.UnicodeUTF8))
        self.dsb_dW0.setSuffix(QtGui.QApplication.translate("MainWindow", " keV", None, QtGui.QApplication.UnicodeUTF8))
        self.label_9.setText(QtGui.QApplication.translate("MainWindow", "Nuclear Structure:", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_C.setText(QtGui.QApplication.translate("MainWindow", "Shape factor", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_isovector.setText(QtGui.QApplication.translate("MainWindow", "Isovector Correction", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_coupled.setToolTip(QtGui.QApplication.translate("MainWindow", "Use actual wave functions calculated using the options specified in the Config file options panel.", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_coupled.setText(QtGui.QApplication.translate("MainWindow", "Coupled", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_C_deformation.setText(QtGui.QApplication.translate("MainWindow", "Deformation", None, QtGui.QApplication.UnicodeUTF8))
        self.cb_relativistic.setText(QtGui.QApplication.translate("MainWindow", "Relativistic", None, QtGui.QApplication.UnicodeUTF8))
        self.b_enable_all.setText(QtGui.QApplication.translate("MainWindow", "Enable All", None, QtGui.QApplication.UnicodeUTF8))
        self.b_disable_all.setText(QtGui.QApplication.translate("MainWindow", "Disable All", None, QtGui.QApplication.UnicodeUTF8))
        self.b_load_profile.setText(QtGui.QApplication.translate("MainWindow", "Load\n"
"Profile", None, QtGui.QApplication.UnicodeUTF8))
        self.b_save_profile.setText(QtGui.QApplication.translate("MainWindow", "Save\n"
"Profile", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page), QtGui.QApplication.translate("MainWindow", "Shape Options", None, QtGui.QApplication.UnicodeUTF8))
        self.label_10.setText(QtGui.QApplication.translate("MainWindow", "TextLabel", None, QtGui.QApplication.UnicodeUTF8))
        self.toolBox.setItemText(self.toolBox.indexOf(self.page_2), QtGui.QApplication.translate("MainWindow", "File Options", None, QtGui.QApplication.UnicodeUTF8))
        self.mainTabWidget.setTabText(self.mainTabWidget.indexOf(self.bsgTab), QtGui.QApplication.translate("MainWindow", "Beta Spectrum Shape Generator", None, QtGui.QApplication.UnicodeUTF8))
        self.mainTabWidget.setTabText(self.mainTabWidget.indexOf(self.nmeTab), QtGui.QApplication.translate("MainWindow", "Nuclear Matrix Element Calculation", None, QtGui.QApplication.UnicodeUTF8))
