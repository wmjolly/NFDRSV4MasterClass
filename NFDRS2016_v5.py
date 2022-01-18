__author__ = 'mjolly'
## \mainpage United States National Fire Danger Rating System
## \details   Contains code required to calculate the United States National Fire Danger Rating System Components and Indices
## \author    W. Matt Jolly, USFS, Rocky Mountain Research Station, Fire Sciences Laboratory
## \contact mjolly@fs.fed.us
## \version   1.0
## \date 03 July 2019
## \copyright GNU Public License.
## \warning Only tested on Python 2.6.2
## \warning Be careful with class memory management in Python
## \remarks  Example usage:\n
# t = USNFDRSCalc(-30.0,"G",0.188503)\n
# t.Update(temp, rh, prcp, ws, yearday)\n
# print t.ERC

import math
import sys, os
from datetime import *
from DeadFuelMoisture import *
from LiveFuelMoisture import *
from NFDRSFuelModels import *


import matplotlib.pyplot as plt
import numpy
INC = 0	#GSI Function Types, Increasing (Tmin, Dayl)
DEC = 1 #GSI Function Types, Decreasing (VPD)
MM_2_IN = 0.0393701  # Conversion factor from mm to inches
KPH_2_MPH = 0.621371 # Conversion factor from KPH to MPH

## \fn fToC Converts Fahrenheit to Celcius
## \param f Temperature in deg F
## \return c Temperature in deg C
def fToc(f):
    return (f -32) * 5/9

## \fn cTof Converts Celcius to Fahrenheit
## \param c Temperature in deg C
## \return f Temperature in deg F
def cTof(c):
    return (c * 9/5) + 32

# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 21:41:05 2019

@author: mjolly
"""

## \class USNFDRSFuelModel
## \brief US NFDRS Fuel Model Class
class USNFDRSFuelModel:
    CTA = 0.046
    SG1 = 2000	## \var SG1 1hr Surface Area to Volume Ratio
    SG10 = 109	## 10hr Surface Area to Volume Ratio
    SG100 = 30	## 100hr Surface Area to Volume Ratio
    SG1000 = 8	## 1000hr Surface Area to Volume Ratio
    SGWOOD = 1	## Wood Surface Area to Volume Ratio
    SGHERB = 1	## Herbaceous Surface Area to Volume Ratio
    L1 = 1		## 1hr fuel loading
    L10 = 1		## 10hr fuel loading
    L100 = 1	## 100hr fuel loading
    L1000 = 1	## 1000hr fuel loading
    LWOOD = 1	## Woody fuel loading
    LHERB = 1	## Herbaceous fuel loading
    DEPTH = 1	## Fuel bed depth
    MXD = 15	## Dead fuel moisture of extinction
    DROUGHT = 0 ## Drought fuel loading
    HD = 8000
    SCM = 301
    WNDFC = 0.6
    FuelModel = ""
    FMCode = ""
    W1 = L1 * CTA
    W1P = 0.0
    W10 = L10 * CTA
    W100 = L100 * CTA
    W1000 = L1000 * CTA
    WWOOD = LWOOD * CTA
    WHERB = LHERB * CTA
    WHERBP = 0.0

    ## \fn __init__(self,FMCode)
    ## \param FMCode Fuel Model Code (A = Grass, F = Brush, G = Timber Understory, Anything else is set to the Slash model)
    def __init__(self,FMCode):
        self.FMCode = FMCode.upper()

        if self.FMCode == "A":	#Grass fuel model
            self.SG1 = 3000
            self.SG10 = 0
            self.SG100 = 0
            self.SG1000 = 0
            self.SGWOOD = 0
            self.SGHERB = 3000
            self.L1 = .2
            self.L10 = 0
            self.L100 = 0
            self.L1000 = 0
            self.LWOOD = 0
            self.LHERB = .3
            self.DEPTH = .8
            self.MXD = 15
            self.HD = 8000
            self.SCM = 301
            self.WNDFC = 0.6
        elif self.FMCode == "B":
            self.SG1 = 700
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1250
            self.SGHERB = 3000
            self.L1 = 3.5
            self.L10 = 4
            self.L100 = 0.5
            self.L1000 = 0
            self.LWOOD = 11.5
            self.LHERB = 0
            self.DEPTH = 4.5
            self.MXD = 15
            self.HD = 9500
            self.SCM = 58
            self.WNDFC = 0.5
        elif self.FMCode == "C":
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2500
            self.L1 = 0.4
            self.L10 = 1
            self.L100 = 0
            self.L1000 = 0
            self.LWOOD = 0.5
            self.LHERB = 0.8
            self.DEPTH = 0.75
            self.MXD = 20
            self.HD = 8000
            self.SCM = 32
            self.WNDFC = 0.5
        elif self.FMCode == "D":
            self.SG1 = 1250
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 1500
            self.L1 = 2
            self.L10 = 1
            self.L100 = 0
            self.L1000 = 0
            self.LWOOD = 3
            self.LHERB = 0.75
            self.DEPTH = 2
            self.MXD = 30
            self.HD = 9000
            self.SCM = 68
            self.WNDFC = 0.4
        elif self.FMCode == "E":
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 1.5
            self.L10 = 0.5
            self.L100 = 0.25
            self.L1000 = 0
            self.LWOOD = 0.5
            self.LHERB = 0.5
            self.DEPTH = 0.4
            self.MXD = 25
            self.HD = 8000
            self.SCM = 25
            self.WNDFC = 0.4
        elif self.FMCode == "F":	# Brush fuel model
            self.SG1 = 700
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 0
            self.SGWOOD = 1250
            self.SGHERB = 0
            self.L1 = 2.5
            self.L10 = 2
            self.L100 = 1.5
            self.L1000 = 0
            self.LWOOD = 9
            self.LHERB = 0
            self.DEPTH = 4.5
            self.MXD = 15
            self.HD = 9500
            self.SCM = 24
            self.WNDFC = 0.5
        elif self.FMCode == "G":	# Timber understory fuel model
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 2.5
            self.L10 = 2
            self.L100 = 5
            self.L1000 = 12
            self.LWOOD = .5
            self.LHERB = .5
            self.DEPTH = 1
            self.MXD = 25
            self.HD = 8000
            self.SCM = 30
            self.WNDFC = 0.4
        elif self.FMCode == "H":	# Timber understory fuel model
            self.L1 = 1.5
            self.L10 = 1
            self.L100 = 2
            self.L1000 = 2
            self.LHERB = .5
            self.LWOOD = .5
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 20
            self.DEPTH = 0.3
            self.WNDFC = 0.4
            self.SCM = 8
        elif self.FMCode == "I":	# Timber understory fuel model
            self.L1 = 12
            self.L10 = 12
            self.L100 = 10
            self.L1000 = 12
            self.LHERB = 0
            self.LWOOD = 0
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 25
            self.DEPTH = 2
            self.WNDFC = 0.5
            self.SCM = 65
        elif self.FMCode == "J":	# Timber understory fuel model
            self.L1 = 7
            self.L10 = 7
            self.L100 = 6
            self.L1000 = 5.5
            self.LHERB = 0
            self.LWOOD = 0
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 25
            self.DEPTH = 1.3
            self.WNDFC = 0.5
            self.SCM = 44
        elif self.FMCode == "K":	# Timber understory fuel model
            self.L1 = 2.5
            self.L10 = 2.5
            self.L100 = 2
            self.L1000 = 2.5
            self.LHERB = 0
            self.LWOOD = 0
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 25
            self.DEPTH = 0.6
            self.WNDFC = 0.5
            self.SCM = 23
        elif self.FMCode == "L":	# Timber understory fuel model
            self.L1 = 0.25
            self.L10 = 0
            self.L100 = 0
            self.L1000 = 0
            self.LHERB = 0.5
            self.LWOOD = 0
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 15
            self.DEPTH = 1
            self.WNDFC = 0.6
            self.SCM = 178
        elif self.FMCode == "N":	# Timber understory fuel model
            self.L1 = 1.5
            self.L10 = 1.5
            self.L100 = 0
            self.L1000 = 0
            self.LHERB = 0
            self.LWOOD = 2
            self.SG1 = 1600
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8700
            self.MXD = 25
            self.DEPTH = 3
            self.WNDFC = 0.6
            self.SCM = 167
        elif self.FMCode == "O":	# Timber understory fuel model
            self.L1 = 2
            self.L10 = 3
            self.L100 = 3
            self.L1000 = 2
            self.LHERB = 0
            self.LWOOD = 7
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 1500
            self.SGWOOD = 1500
            self.HD = 9000
            self.MXD = 30
            self.DEPTH = 4
            self.WNDFC = 0.5
            self.SCM = 99
        elif self.FMCode == "P":	# Timber understory fuel model
            self.L1 = 1
            self.L10 = 1
            self.L100 = 0.5
            self.L1000 = 2
            self.LHERB = 0.5
            self.LWOOD = 0.5
            self.SG1 = 1750
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 30
            self.DEPTH = 0.4
            self.WNDFC = 0.4
            self.SCM = 14
        elif self.FMCode == "Q":	# Timber understory fuel model
            self.L1 = 2
            self.L10 = 2.5
            self.L100 = 2
            self.L1000 = 1
            self.LHERB = 0.5
            self.LWOOD = 4
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 1500
            self.SGWOOD = 1200
            self.HD = 8000
            self.MXD = 25
            self.DEPTH = 3
            self.WNDFC = 0.4
            self.SCM = 59
        elif self.FMCode == "R":	# Timber understory fuel model
            self.L1 = 0.5
            self.L10 = 0.5
            self.L100 = 0.5
            self.L1000 = 0
            self.LHERB = 0.5
            self.LWOOD = 0.5
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 25
            self.DEPTH = 0.25
            self.WNDFC = 0.4
            self.SCM = 6
        elif self.FMCode == "S":	# Timber understory fuel model
            self.L1 = 0.5
            self.L10 = 0.5
            self.L100 = 0.5
            self.L1000 = 0.5
            self.LHERB = 0.5
            self.LWOOD = 0.5
            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 1500
            self.SGWOOD = 1200
            self.HD = 8000
            self.MXD = 25
            self.DEPTH = 0.4
            self.WNDFC = 0.6
            self.SCM = 17
        elif self.FMCode == "T":	# Timber understory fuel model
            self.L1 = 1
            self.L10 = 0.5
            self.L100 = 0
            self.L1000 = 0
            self.LHERB = 0.5
            self.LWOOD = 2.5
            self.SG1 = 2500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 15
            self.DEPTH = 1.25
            self.WNDFC = 0.6
            self.SCM = 96
        elif self.FMCode == "U":	# Western Pines
            self.L1 = 1.5
            self.L10 = 1.5
            self.L100 = 1
            self.L1000 = 0
            self.LHERB = 0.5
            self.LWOOD = 0.5
            self.SG1 = 1750
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGHERB = 2000
            self.SGWOOD = 1500
            self.HD = 8000
            self.MXD = 20
            self.DEPTH = 0.5
            self.WNDFC = 0.4
            self.SCM = 16
        
        elif self.FMCode == "V":  # Grass
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 0.1
            self.L10 = 0.0
            self.L100 = 0.0
            self.L1000 = 0.0
            self.LWOOD = 0.0
            self.LHERB = 1.0
            self.DEPTH = 1
            self.MXD = 15
            self.HD = 8000
            self.SCM = 108
            self.WNDFC = 0.6
            self.DROUGHT = 0

        elif self.FMCode == "W":  # Grass
            
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 0.5
            self.L10 = 0.5
            self.L100 = 0.0
            self.L1000 = 0.0
            self.LWOOD = 1.0
            self.LHERB = 0.6
            self.DEPTH = 1.5
            self.MXD = 15
            self.HD = 8000
            self.SCM = 62
            self.WNDFC = 0.4
            self.DROUGHT = 1
        elif self.FMCode == "X":  # Grass
            print("X")
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 4.5
            self.L10 = 2.45
            self.L100 = 0.0
            self.L1000 = 0.0
            self.LWOOD = 7.0
            self.LHERB = 1.55
            self.DEPTH = 4.4
            self.MXD = 25
            self.HD = 8000
            self.SCM = 104
            self.WNDFC = 0.4
            self.DROUGHT = 2.5


        elif self.FMCode == "Y":  # Timber understory fuel model
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 2.5
            self.L10 = 2.2
            self.L100 = 3.6
            self.L1000 = 10.16
            self.LWOOD = 0
            self.LHERB = 0.3
            self.DEPTH = 0.6
            self.MXD = 25
            self.HD = 8000
            self.SCM = 5
            self.WNDFC = 0.2
            self.DROUGHT = 5
 
        elif self.FMCode == "Z":  # Grass
            self.SG1 = 2000
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 1500
            self.SGHERB = 2000
            self.L1 = 4.5
            self.L10 = 4.25
            self.L100 = 4.0
            self.L1000 = 4.0
            self.LWOOD = 0.0
            self.LHERB = 0.0
            self.DEPTH = 1.5
            self.MXD = 25
            self.HD = 8000
            self.SCM = 19
            self.WNDFC = 0.4
            self.DROUGHT = 7
        else:	# Assign the Slash Fuel Model Parameters

            self.SG1 = 1500
            self.SG10 = 109
            self.SG100 = 30
            self.SG1000 = 8
            self.SGWOOD = 0
            self.SGHERB = 0
            self.L1 = 7
            self.L10 = 7
            self.L100 = 6
            self.L1000 = 5.5
            self.LWOOD = 0
            self.LHERB = 0
            self.DEPTH = 1.3
            self.MXD = 25
            self.HD = 8000
            self.SCM = 44
            self.WNDFC = 0.5

        if (self.SG1 <= 0):
            self.SG1 = 2000
        if (self.SG10 <= 0):
            self.SG10 = 109
        if (self.SG100 <=0 ):
            self.SG100 = 30
        if (self.SG1000 <= 0):
            self.SG1000 = 8
        if (self.SGWOOD <= 0):
            self.SGWOOD = 1
        if (self.SGHERB <= 0):
            self.SGHERB = 1
        self.W1 = self.L1 * self.CTA
        self.W10 = self.L10 * self.CTA
        self.W100 = self.L100 * self.CTA
        self.W1000 = self.L1000 * self.CTA
        self.WWOOD = self.LWOOD * self.CTA
        self.WHERB = self.LHERB * self.CTA
        
## \class USNFDRSCalc
## \brief US NFDRS Components and Indices Calculator Class
class USNFDRSCalc:

    MC1 = 5
    MC10 = 5
    MC100 = 20
    MC1000 = 15
    MCHERB = 100
    MCWOOD = 70
    
    Hist1000 = [15,15,15,15,15,15,15]
    HistBndryT =[15,15,15,15,15,15,15]
    KBDI = 100
    yKBDI = 100
    CummPrecip = 0.0
    fAvgPrecip = 13.36
    Y100 = 20
    Y1000 = 15
    YGSI = 0.0
    GSI = 0.0
    Lat = -30
    SlopeCls = 1
    SOW = 0
    NFDRSVersion = 78
    m_HerbMaxGSI = 1.0
    m_HerbGreenup = 0.5
    m_HerbMax = 250.0
    m_HerbMin = 30
    m_WoodyMaxGSI = 1.0
    m_WoodyGreenup = 0.5
    m_WoodyMax = 200.00
    m_WoodyMin = 60.0
    m_HerbSlope = 1.0
    m_HerbIntercept = 1.0
    m_WoodySlope = 1.0
    m_WoodyIntercept = 1.0
    MCHERB = m_HerbMin
    MCWOOD = m_WoodyMin
    MetSet = False
    beta = 0.188503 #Parameter relating precipitation amount to precipitation duration
    nConsectiveSnowDays = 0
    YesterdayJDay = 1
    num_updates = 0
    SnowCovered = False
    FuelModel = USNFDRSFuelModel("Y")
    OneHourFM = DeadFuelMoisture(0.64,"1hr")
    
    TenHourFM = DeadFuelMoisture(0.64,"10hr")
    HundredHourFM = DeadFuelMoisture(2.0,"100hr")
    ThousandHourFM = DeadFuelMoisture(3.81,"1000hr")
    GsiFM = LiveFuelMoisture(45,False,False)
    HerbFM = LiveFuelMoisture(45,True,False)
    WoodFM = LiveFuelMoisture(45,False,False)
    OneHourFM.setAdsorptionRate(0.462252733)
    TenHourFM.setAdsorptionRate(0.079548303)
    HundredHourFM.setAdsorptionRate(0.06)
    ThousandHourFM.setAdsorptionRate(0.06)
		
	#//these two lines were not done in the original constructor. 
	#/sets maximum stick moisture limit. Not really relevant for 100 and 1000 hour moistures.
    OneHourFM.setMaximumLocalMoisture(0.35)
    TenHourFM.setMaximumLocalMoisture(0.35)
    HundredHourFM.setMaximumLocalMoisture(0.35)
    ThousandHourFM.setMaximumLocalMoisture(0.35)
    
    #precalculate slopes and intercepts for GSI->LFM Conversion
    if(m_HerbGreenup == 1.0):
        m_HerbGreenup = 0.9999
    if(m_WoodyGreenup == 1.0):
        m_WoodyGreenup = 0.9999
    m_HerbSlope = (m_HerbMax - m_HerbMin) / (1.0 - m_HerbGreenup)
    m_HerbIntercept = m_HerbMax - m_HerbSlope
    m_WoodySlope = (m_WoodyMax - m_WoodyMin) / (1.0 - m_WoodyGreenup)
    m_WoodyIntercept = m_WoodyMax - m_WoodySlope


    # Calculated Indices
    SC = 0
    ERC = 0.0	## \var ERC Energy Release Component
    BI = 0.0
    def Init(self):
        MC1 = 5
        MC10 = 5
        MC100 = 20
        MC1000 = 15
        MCHERB = 100
        MCWOOD = 70
        
        Hist1000 = [15,15,15,15,15,15,15]
        HistBndryT =[15,15,15,15,15,15,15]
        KBDI = 100
        yKBDI = 100
        CummPrecip = 0.0
        fAvgPrecip = 13.36
        Y100 = 20
        Y1000 = 15
        YGSI = 0.0
        GSI = 0.0
        Lat = -30
        SlopeCls = 1
        SOW = 0
        NFDRSVersion = 78
        m_HerbMaxGSI = 1.0
        m_HerbGreenup = 0.5
        m_HerbMax = 250.0
        m_HerbMin = 30
        m_WoodyMaxGSI = 1.0
        m_WoodyGreenup = 0.5
        m_WoodyMax = 200.00
        m_WoodyMin = 60.0
        m_HerbSlope = 1.0
        m_HerbIntercept = 1.0
        m_WoodySlope = 1.0
        m_WoodyIntercept = 1.0
        MCHERB = m_HerbMin
        MCWOOD = m_WoodyMin
        MetSet = False
        beta = 0.188503 #Parameter relating precipitation amount to precipitation duration
        nConsectiveSnowDays = 0
        YesterdayJDay = 1
        num_updates = 0
        SnowCovered = False
        FuelModel = USNFDRSFuelModel("Y")
        OneHourFM = DeadFuelMoisture(0.64,"1hr")
        
        TenHourFM = DeadFuelMoisture(0.64,"10hr")
        HundredHourFM = DeadFuelMoisture(2.0,"100hr")
        ThousandHourFM = DeadFuelMoisture(3.81,"1000hr")
        GsiFM = LiveFuelMoisture(45,False,False)
        HerbFM = LiveFuelMoisture(45,True,False)
        WoodFM = LiveFuelMoisture(45,False,False)
        OneHourFM.setAdsorptionRate(0.462252733)
        TenHourFM.setAdsorptionRate(0.079548303)
        HundredHourFM.setAdsorptionRate(0.06)
        ThousandHourFM.setAdsorptionRate(0.06)
    		
    	#//these two lines were not done in the original constructor. 
    	#/sets maximum stick moisture limit. Not really relevant for 100 and 1000 hour moistures.
        OneHourFM.setMaximumLocalMoisture(0.35)
        TenHourFM.setMaximumLocalMoisture(0.35)
        HundredHourFM.setMaximumLocalMoisture(0.35)
        ThousandHourFM.setMaximumLocalMoisture(0.35)
        
        #precalculate slopes and intercepts for GSI->LFM Conversion
        if(m_HerbGreenup == 1.0):
            m_HerbGreenup = 0.9999
        if(m_WoodyGreenup == 1.0):
            m_WoodyGreenup = 0.9999
        m_HerbSlope = (m_HerbMax - m_HerbMin) / (1.0 - m_HerbGreenup)
        m_HerbIntercept = m_HerbMax - m_HerbSlope
        m_WoodySlope = (m_WoodyMax - m_WoodyMin) / (1.0 - m_WoodyGreenup)
        m_WoodyIntercept = m_WoodyMax - m_WoodySlope
    
    
        # Calculated Indices
        SC = 0
        ERC = 0.0	## \var ERC Energy Release Component
        BI = 0.0
    ## \fn  __init__(self,Lat,FMCode,beta)
    ## \brief Constructor for the USNFDRSCalc Class
    ## \param Lat Latitude (DD)
    ## \param FMCode Fuel Model Code (valid values are: A, F, G, all others values are set to J)
    ## \remark Example: t = USNFDRSCalc(-30.0,"G",0.188503)
    def __init__(self,Lat,FMCode,beta):
        
        self.Lat = Lat
        self.FuelModel = USNFDRSFuelModel(FMCode)
        self.beta = beta
        self.Init()

    ## \fn eqmc (self,Temp,RH)
    ## \brief Calculate the equilibrium moisture content
    ## \param Temp Temperature (deg F)
    def eqmc (self,Temp,RH):

        if(RH > 50):
            return  (21.0606 + 0.005565 * RH**2 - 0.00035 * RH * Temp - 0.483199 * RH)
        if(RH > 10) and (RH < 51):
            return (2.22749 + 0.160107 * RH - 0.014784 * Temp)
        else:
            return (0.03229 + 0.281073 * RH - 0.000578 * RH * Temp)

    ## \fn satvap (self,Temp)
    ## \brief Calculate the Saturation Vapor Pressure from the temperature passed
    ## \param Temp Temperature (deg C)
    def satvap (self,Temp):
        kTemp = Temp  + 273.16
        return exp(1.81 + (kTemp * 17.27 - 4717.31) / (kTemp - 35.86))

    ## \fn oneten (self, Temp, RH, SOW)
    ## \brief Calculate the one and ten hour moisture contents
    ## \param Temp Temperature (deg F)
    ## \param RH Relative Humidity (%)
    ## \param SOW State of the Weather (0-9)
    def oneten (self, Temp, RH, SOW):
        tfact = float(0)
        hfact = float(0)
        emc = float(0)

        # Determine the temperature and rh factors to adjust for fuel temperatures
        if SOW == 0:
            tfact = 25.0
            hfact = 0.75
        if SOW == 1:
            tfact = 19.0
            hfact = 0.83
        if SOW == 2:
            tfact = 12.0
            hfact = 0.92
        else:
            tfact = 5.0
            hfact = 1.0

        emc = self.eqmc(tfact + Temp, hfact * RH)
        self.MC1 = 1.03 * emc
        self.MC10 = 1.28 * emc
        self.SOW = SOW

    ## \fn CalcDaylightHours(self,Julian)
    ## \brief Calculate Daylight hours as a function of Latitude and Julian Date
    ## \param Julian Julian Year Day (1-366)
    ## Requires that Latitude already be set as part of the class intialization
    def CalcDaylightHours(self,Julian):
        phi = 0.0
        xfact = 0.0
        decl = 0.0
        tla = 0.0
        phi = tan(self.Lat * 0.01745) * -1.0
        xfact = (Julian - 80) * 0.01745
        decl = 23.5 * sin(xfact)
        decl = decl * 0.01745
        tla = phi * sin(decl)

        if (fabs(tla) < .01):
            tla = 0.01
        else:
            if (tla >= 1.0):
                tla = 0.99999999
            else:
                if (tla <= -1.0):
                    tla = -.9999999


        tla = atan(sqrt((1.0 - tla * tla))/tla)
        if (tla < 0.0):
            tla = tla + 3.141593
        return tla * 7.64

    ## \fn hundredthous (self, Temp,RH,MaxTemp,MaxRH,MinTemp, MinRH,Julian,PrecipDur)
    ## \brief Calculate the hundred and thousand hour moisture contents
    ## \param Temp Temperature (deg F)
    ## \param RH Relative Humidity (%)
    ## \param MaxTemp 24-hour Maximum Temperature (deg F)
    ## \param MaxRH 24-hour Maximum Relative Humidity (%)
    ## \param MinTemp 24-hour Minimum Temperature (deg F)
    ## \param MinRH 24-hour Minimum Relative Humidity (%)
    ## \param Julian Julian Year Day (1-366)
    ## \param 24-hour Precipitation Duration (hours)
    def hundredthous (self, MaxTemp,MaxRH,MinTemp, MinRH,Julian,PrecipDur):
        emcMin = 0.0
        emcMax = 0.0
        emcBar = 0.0
        bndryH = 0.0
        bndryT = 0.0
        bndryBar = 0.0
        Daylight = 0.0
        ambvp = 0.0

        Daylight = self.CalcDaylightHours( Julian)
        emcMin = self.eqmc (MaxTemp, MinRH)
        emcMax = self.eqmc (MinTemp, MaxRH)
        emcBar = (Daylight * emcMin + (24 - Daylight) * emcMax) / 24
        bndryH = ((24.0 - PrecipDur) * emcBar + PrecipDur * (.5 * PrecipDur + 41)) / 24
        self.MC100 = self.Y100 + .3156 * (bndryH - self.Y100);
        bndryT = ((24.0 - PrecipDur) * emcBar + PrecipDur * (2.7 * PrecipDur + 76)) / 24

        # Note this can be rewritten using Python push and pop
        for i in range (0,6):
            self.HistBndryT[i] = self.HistBndryT[i+1]
            bndryBar = bndryBar + self.HistBndryT[i]

        self.HistBndryT[6] = bndryT
        bndryBar = (bndryBar + bndryT) / 7
        self.MC1000 = self.Hist1000[0] + (bndryBar - self.Hist1000[0]) * .3068
        for i in range (0,6):
            self.Hist1000[i] = self.Hist1000[i+1]

        self.Hist1000[6] = self.MC1000
        self.Y100 = self.MC100
        self.Y1000 = self.MC1000
    ## \fn  GetInd(self,Val,ValMin, ValMax,Type)
    ## \brief Indicator functions for calculating the Growing Season Index
    ## \param Val Input Value
    ## \param ValMin Indicator Minimum Value
    ## \param ValMax Indicator Maximum Value
    ## \param Type INC (Increasing) DEC (Decreasing) Defined Globally
    def GetInd(self,Val,ValMin, ValMax,Type):

        if(Type == DEC):	# A Decreasing indicator function like the one used for VPD
            if( Val < ValMin):
                return 1
            elif(Val > ValMax):
                return 0
            else:
                return 1 - (Val - ValMin) / (ValMax - ValMin);

        elif(Type == INC): # An increasing indicator function like the one used for TMIN and DAYL
            if( Val < ValMin):
                return 0
            elif(Val > ValMax):
                return 1
            else:
                return (Val - ValMin) / (ValMax - ValMin)
        self.FM.W1P = self.FM.W1 + self.FM.WHERB * fctCur
        self.FM.WHERBP = self.FM.WHERB * (1 - fctCur)
    ## \fn live ( self, Temp, RH, Julian)
    ## \brief Calculate the live fuel moistures from the Growing Season Index
    ## \param Temp Temperature (deg C)
    ## \param RH Relative Humidity (%)
    ## \param Julian Julian Year Day (1-366)
    # Calculate the Growing Season Index and subsequent live fuel moistures
    def live ( self, Temp, RH, Julian):
        GSIINT = 21	#Running average period for GSI caculations
        # Caculate the Vapor Pressure Deficit indicator
        SatVP = self.satvap(Temp)
        VPD = SatVP - (RH /100) * SatVP
        indVPD = self.GetInd(VPD,1000,3000,DEC)

        # Calculate the daylength indicator
        Photo = self.CalcDaylightHours (Julian) * 3600
        indPhoto = self.GetInd(Photo,36000,39600,INC)

        #  Calculate the temperature indicator
        indTemp = self.GetInd(Temp,5,15,INC)
        iGSI = indTemp * indVPD * indPhoto

        # Calculate the weighted-average of the GSI (based on a 21 day smoothing period)
        self.GSI = ((GSIINT -1) * self.YGSI + iGSI) / GSIINT
        self.CalcWoodyFMFromGSI()
        self.CalcHerbFMFromGSI()

        # Store today's GSI for use in tomorrow's calculation
        self.YGSI = self.GSI

    ## \fn CalcHerbFMFromGSI(self)
    ## \brief Calculate Herbaceous Fuel Moisture from GSI
    def CalcHerbFMFromGSI(self):

        rescale = self.GSI / self.m_HerbMaxGSI
        ret = self.m_HerbMin
        rescale = min(1.0, rescale)
        rescale = max(0.0, rescale)
        if(rescale >= self.m_HerbGreenup):

            self.MCHERB = self.m_HerbSlope * rescale + self.m_HerbIntercept
        else:
            self.MCHERB = self.m_HerbMin

    ## \fn CalcWoodyFMFromGSI(self)
    ## \brief Calculate Woody Fuel Moisture from GSI
    def CalcWoodyFMFromGSI(self):

        rescale = self.GSI / self.m_WoodyMaxGSI
        rescale = min(1.0, rescale)
        rescale = max(0.0, rescale)
        if(rescale >= self.m_WoodyGreenup):

            self.MCWOOD = self.m_WoodySlope * rescale + self.m_WoodyIntercept
        else:
            self.MCWOOD = self.m_WoodyMin

    ## \fn CureHerb (self)
    ## \brief Calculate the percentage of load transfer due to herbaceous curing
    def CureHerb (self):
        fctCur = 0.0
        fctCur = 1.33 - .0111 * self.MCHERB
        if (fctCur < 0):
            fctCur = 0.0
        if (fctCur > 1):
            fctCur = 1.0
        self.FuelModel.W1P = self.FuelModel.W1 + self.FuelModel.WHERB * fctCur
        self.FuelModel.WHERBP = self.FuelModel.WHERB * (1 - fctCur)

   
    ## \fn SetMoistures(self,MC1,MC10,MC100,MC1000,MCWOOD,MCHERB)
    ## \brief Set the fuel moistures for a fire danger calculation
    ## \param MC1 1-hour Moisture Content (% dry wt)
    ## \param MC10 10-hour Moisture Content (% dry wt)
    ## \param MC100 100-hour Moisture Content (% dry wt)
    ## \param MC1000 1000-hour Moisture Content (% dry wt)
    ## \param MCWOOD Live Woody Moisture Content (% dry wt)
    ## \param MCHERB Live Herbaceous Moisture Content (% dry wt)
    def SetMoistures(self,MC1,MC10,MC100,MC1000,MCWOOD,MCHERB):
        self.MC1 = MC1
        self.MC10 = MC10
        self.MC100 = MC100
        self.MC1000 = MC1000
        self.MCHERB = MCHERB
        self.MCWOOD = MCWOOD
        
        
        
    
    def CalcKBDI (self,fPrecipAmt, iMaxTemp):
        net = 0
        idq = 0
        pptnet = 0.00
        xkbdi = 0.00
        xtemp = 0.00
        
        KBDI = int(self.yKBDI)
        if(fPrecipAmt == 0.0):
           self.CummPrecip = 0
        else:
            if(self.CummPrecip > 0.20):
                pptnet = fPrecipAmt
                self.CummPrecip = self.CummPrecip + fPrecipAmt
            else:
                self.CummPrecip = self.CummPrecip + fPrecipAmt
                if(self.CummPrecip > 0.20):
                    pptnet = self.CummPrecip - 0.20
        
        net = (100.0 * pptnet) + 0.0005
        
        net = KBDI - net
        
        if net > 0:
            KBDI = round(net,0)
        else:
            KBDI = 0
        if(iMaxTemp > 50):
            idq = (800.0 - KBDI) * (0.9679 * exp(0.0486 * iMaxTemp) - 8.299) * 0.001 / (1.0 + 10.88 * exp(-0.04409 * self.fAvgPrecip)) + 0.5
        self.KBDI = int(KBDI + idq)
        self.yKBDI = self.KBDI
     
    
    ## \fn CalcIndexes (self,WS)
    ## \brief Calculates the Energy Release Component, Spread Component and Burning Index
    ##
    ## \param WS Windspeed (KPH)
    def CalcIndexes (self,WS,KBDI,FuelTemperature):

        WS = ceil(WS)
        STD = .0555
        STL = .0555
        RHOD = 32
        RHOL = 32
        ETASD = 0.4173969
        ETASL = 0.4173969
        SG1 = self.FuelModel.SG1
        SG10 = self.FuelModel.SG10
        SG100 = self.FuelModel.SG100
        SG1000 = self.FuelModel.SG1000
        SGHERB = self.FuelModel.SGHERB
        SGWOOD = self.FuelModel.SGWOOD

        W1 = self.FuelModel.L1 * self.FuelModel.CTA
        W1P = 0.0
        W10 = self.FuelModel.L10 * self.FuelModel.CTA
        W100 = self.FuelModel.L100 * self.FuelModel.CTA
        W1000 = self.FuelModel.L100 * self.FuelModel.CTA
        WWOOD = self.FuelModel.LWOOD * self.FuelModel.CTA
        WHERB = self.FuelModel.LHERB    * self.FuelModel.CTA



        if (KBDI > 1000):
            WTOTD = W1 + W10 + W100
            WTOTL = WHERB + WWOOD
            WTOT = WTOTD + WTOTL
            PackingRatio = WTOT / self.FuelModel.DEPTH
            if (PackingRatio == 0):
                PackingRatio = 1.0
            WTOTD = WTOTD + W1000
            DroughtUnit = self.FuelModel.DROUGHT / 700.
            W1 = W1 + (W1 / WTOTD) * (KBDI - 100) * DroughtUnit
            W10 = W10 + (W10 / WTOTD) * (KBDI - 100) * DroughtUnit
            W100 = W100 + (W100 / WTOTD) * (KBDI - 100) * DroughtUnit
            W1000 = W1000 + (W1000 / WTOTD) * (KBDI - 100) * DroughtUnit
            WTOT = W1 + W10 + W100 + W1000 + WTOTL
            self.FuelModel.DEPTH = (WTOT - W1000) / PackingRatio
        # Transfer herbaceous load due to curing
        self.CureHerb()
        WTOTD = self.FuelModel.W1P + self.FuelModel.W10 + self.FuelModel.W100 + self.FuelModel.W1000
        WTOTL = self.FuelModel.WHERBP + self.FuelModel.WWOOD
        WTOT = WTOTD + WTOTL

        W1N = self.FuelModel.W1P * (1.0 - STD)
        W10N = self.FuelModel.W10 * (1.0 - STD)
        W100N = self.FuelModel.W100 * (1.0 - STD)
        WHERBN = self.FuelModel.WHERBP * (1.0 - STL)
        WWOODN = self.FuelModel.WWOOD * (1.0 - STL)
        WTOTLN = WTOTL * (1.0 - STL)
        RHOBED = (WTOT - self.FuelModel.W1000) / self.FuelModel.DEPTH
        RHOBAR = ((WTOTL * RHOL) + (WTOTD * RHOD)) / WTOT
        BETBAR = RHOBED / RHOBAR

        # If Total Live Loading is greater than zero, calculate the live moisture of extinction
        if (WTOTLN > 0):


            HN1 = W1N * exp(-138.0 / SG1)
            HN10 = W10N  * exp(-138.0 / SG10)
            HN100 = W100N * exp(-138.0 / SG100)

            if ((-500 / SGHERB) < -180.218):
                HNHERB = 0
            else:
                HNHERB = WHERBN * exp(-500.0 / SGHERB)

            if ((-500 / SGWOOD) < -180.218):
                HNWOOD = 0
            else:
                HNWOOD = WWOODN * exp(-500.0 / SGWOOD)

            if ((HNHERB + HNWOOD) == 0):
                WRAT = 0
            else:
                WRAT = (HN1 + HN10 + HN100) / (HNHERB + HNWOOD)

            MCLFE = ((self.MC1 * HN1) + (self.MC10 * HN10) + (self.MC100 * HN100)) / (HN1 + HN10 + HN100);
            MXL = (2.9 * WRAT * (1.0 - MCLFE / self.FuelModel.MXD) - 0.226) * 100

        else:
            MXL = 0


        if (MXL < self.FuelModel.MXD):
            MXL = self.FuelModel.MXD

        SA1 = (self.FuelModel.W1P / RHOD) * SG1
        SA10 = (self.FuelModel.W10 / RHOD) * SG10
        SA100 = (self.FuelModel.W100 / RHOD) * SG100
        SAHERB = (self.FuelModel.WHERBP / RHOL) * SGHERB
        SAWOOD = (self.FuelModel.WWOOD / RHOL) * SGWOOD
        SADEAD = SA1 + SA10 + SA100
        SALIVE = SAHERB + SAWOOD



        F1 = SA1 / SADEAD
        F10 = SA10 / SADEAD
        F100 = SA100 / SADEAD

        if (WTOTL <=0):
            FHERB = 0
            FWOOD = 0

        else:
            FHERB = SAHERB / SALIVE
            FWOOD = SAWOOD / SALIVE

        FDEAD = SADEAD / (SADEAD + SALIVE)
        FLIVE = SALIVE / (SADEAD + SALIVE)

        WDEADN = (F1 * W1N) + (F10 * W10N) + (F100 * W100N);

        if (SGWOOD > 1200 and SGHERB > 1200):
            WLIVEN = WTOTLN
        else:
            WLIVEN = (FWOOD * WWOODN) + (FHERB * WHERBN)


        SGBRD = (F1 * SG1) + (F10 * SG10) + (F100 * SG100)
        SGBRL = (FHERB * SGHERB) + (FWOOD * SGWOOD)
        SGBRT = (FDEAD * SGBRD) + (FLIVE * SGBRL)
        BETOP = 3.348 * pow(SGBRT, -0.8189)
        GMAMX = pow(SGBRT, 1.5) / (495.0 + 0.0594 * pow(SGBRT, 1.5))
        AD = 133 * pow(SGBRT, -0.7913)
        GMAOP = GMAMX * pow( (BETBAR/BETOP), AD) * exp(AD * (1.0 - (BETBAR / BETOP)))

        ZETA = exp((0.792 + 0.681 * pow(SGBRT, 0.5)) * (BETBAR + 0.1))
        ZETA = ZETA / (192.0 + 0.2595 * SGBRT)

        WTMCD = (F1 * self.MC1) + (F10 * self.MC10) + (F100 * self.MC100)
        WTMCL = (FHERB * self.MCHERB) + (FWOOD * self.MCWOOD)
        DEDRT = WTMCD / self.FuelModel.MXD
        LIVRT = WTMCL / MXL
        ETAMD = 1.0 - 2.59 * DEDRT + 5.11 * pow(DEDRT,2.0) - 3.52 * pow(DEDRT, 3.0)
        ETAML = 1.0 - 2.59 * LIVRT + 5.11 * pow(LIVRT,2.0) - 3.52 * pow(LIVRT, 3.0)
        if (ETAMD < 0):
            ETAMD = 0
        if (ETAMD > 1):
            ETAMD = 1
        if (ETAML < 0):
            ETAML = 0
        if (ETAML > 1):
            ETAML = 1

        B = 0.02526 * pow(SGBRT, 0.54)
        C = 7.47 * exp(-0.133 * pow(SGBRT,0.55))
        E = 0.715 * exp(-3.59 * pow(10.0, -4.0) * SGBRT)


        UFACT = C * pow(BETBAR / BETOP, -1 * E)
        IR = GMAOP * ((WDEADN * self.FuelModel.HD * ETASD * ETAMD) + (WLIVEN * self.FuelModel.HD * ETASL * ETAML))


        if (88.0 * WS * self.FuelModel.WNDFC > 0.9 * IR):

            PHIWND = UFACT * pow(0.9 * IR, B)
        else:
            PHIWND = UFACT * pow(WS * 88.0 * self.FuelModel.WNDFC, B)


        if(self.SlopeCls ==  1):
            PHISLP = 0.267 * pow(BETBAR, -0.3)
        elif(self.SlopeCls == 2):
            PHISLP = 0.533 * pow(BETBAR, -0.3)
        elif(self.SlopeCls == 3):
            PHISLP = 1.068 * pow(BETBAR, -0.3)
        elif(self.SlopeCls == 4):
            PHISLP = 2.134 * pow(BETBAR, -0.3)
        else:
            PHISLP = 4.273 * pow(BETBAR, -0.3)

        XF1 = F1 * exp(-138.0 / SG1) * (250.0 + 11.16 * self.MC1)
        XF10 = F10 * exp(-138.0 / SG10) * (250.0 + 11.16 * self.MC10)
        XF100 = F100 * exp(-138.0 / SG100) * (250.0 + 11.16 * self.MC100)
        XFHERB = FHERB * exp(-138.0 / SGHERB) * (250.0 + 11.16 * self.MCHERB)
        XFWOOD = FWOOD * exp(-138.0 / SGWOOD) * (250.0 + 11.16 * self.MCWOOD)
        HTSINK = RHOBED * (FDEAD * (XF1 + XF10 + XF100) + FLIVE * (XFHERB + XFWOOD))

        fROS = IR * ZETA * (1.0 + PHISLP + PHIWND) / HTSINK
        self.SC = ceil(fROS)
        F1E = self.FuelModel.W1P / WTOTD
        F10E = self.FuelModel.W10 / WTOTD
        F100E = self.FuelModel.W100 / WTOTD
        F1000E = self.FuelModel.W1000 / WTOTD
        if (WTOTL <=0):
            FHERBE = 0
            FWOODE = 0
        else:
            FHERBE = self.FuelModel.WHERBP / WTOTL
            FWOODE = self.FuelModel.WWOOD / WTOTL

        FDEADE = WTOTD / WTOT
        FLIVEE = WTOTL / WTOT

        WDEDNE = WTOTD * (1.0 - STD)
        WLIVNE = WTOTL * (1.0 - STL)

        SGBRDE = (F1E * SG1) + (F10E * SG10) + (F100E * SG100) + (F1000E * SG1000)
        SGBRLE = (FHERBE * SGHERB) + (FWOODE * SGWOOD)
        SGBRTE = (FDEADE * SGBRDE) + (FLIVEE * SGBRLE)
        BETOPE = 3.348 * pow(SGBRTE, -0.8189)
        GMAMXE = pow(SGBRTE, 1.5) / (495.0 + 0.0594 * pow(SGBRTE, 1.5))
        ADE = 133 * pow(SGBRTE, -0.7913)
        GMAOPE = GMAMXE * pow( (BETBAR/BETOPE), ADE) * exp(ADE * (1.0 - (BETBAR / BETOPE)))
        WTMCDE = (F1E * self.MC1) + (F10E * self.MC10) + (F100E * self.MC100) + (F1000E * self.MC1000)
        WTMCLE = (FHERBE * self.MCHERB) + (FWOODE * self.MCWOOD)
        DEDRTE = WTMCDE / self.FuelModel.MXD
        LIVRTE = WTMCLE / MXL
        ETAMDE = 1.0 - 2.0 * DEDRTE + 1.5 * pow(DEDRTE,2.0) - 0.5 * pow(DEDRTE, 3.0)
        ETAMLE = 1.0 - 2.0 * LIVRTE + 1.5 * pow(LIVRTE,2.0) - 0.5 * pow(LIVRTE, 3.0)
        if (ETAMDE < 0):
            ETAMDE = 0
        if (ETAMDE > 1):
            ETAMDE = 1
        if (ETAMLE < 0):
            ETAMLE = 0
        if (ETAMLE > 1):
            ETAMLE = 1

        IRE = (FDEADE * WDEDNE * self.FuelModel.HD * ETASD * ETAMDE)
        IRE = GMAOPE * (IRE + (FLIVEE * WLIVNE * self.FuelModel.HD * ETASL * ETAMLE))
        TAU = 384.0 / SGBRT

        self.ERC = 0.04 * IRE * TAU

        self.BI = 10.0 * .301 * (fROS * self.ERC)**0.46

        if ((self.SOW >= 5) and (self.SOW <= 7)):
            self.SC = 0
            self.BI = 0
        # Finally, calculate the Igntion Component
        TMPPRM = 0.0
        PNORM1 = 0.00232
        PNORM2 = 0.99767
        if (self.FuelModel.SCM <= 0):
            IC = 0
        TMPPRM = FuelTemperature
        QIGN = 144.5 - (0.266 * TMPPRM) - (0.00058 * TMPPRM * TMPPRM)   - (0.01 * TMPPRM * self.MC1) + 18.54 * (1.0 - exp(-0.151 * self.MC1)) + 6.4 * self.MC1
        #print(QIGN)
        if (QIGN >= 344.0):
            IC = 0
        
        CHI = (344.0 - QIGN) / 10.0
        #print(CHI)
        if (((CHI** 3.66) * 0.000923 / 50) <= PNORM1):
            IC = 0

        PI = (((CHI** 3.66) * 0.000923 / 50) - PNORM1) * 100.0 / PNORM2
        if (PI < 0):
            PI = 0
        if (PI > 100):
            PI = 100
        SCN = 100.0 * self.SC / self.FuelModel.SCM
        if (SCN > 100.0):
            SCN = 100.0
        PFI = (SCN** 0.5)
        IC = 0.10 * PI * PFI
        if (self.SC < 0.00001):
            self.IC = 0
            

    ## \fn PDurFunc(self,Precip)
    ## \brief Function to derive precipitation duration from precipitation amount
    ## \param pamnt Precipitation Amount (mm)
    ## \param bw Beta parameter for logistic function
    ## \return Returns precipitation amount in hours, bounded between 0 and 24
    def PDurFunc(self,Precip):

        Precip = Precip * MM_2_IN
        if(Precip == 0):

            return 0

        else:

            return ceil(24 * (1 - exp(-self.beta*Precip)))



    ## \fn Update(self, Temp, RH, Precip, WS, JDay)
    ## \brief Update the fire weather and recaculate the fire danger indices
    ## \param Temp Air Temperature (deg C)
    ## \param RH Relative Humidity (%)
    ## \param Precip Precipitation Amount (mm)
    ## \param WS Windspeed (KPH)
    ## \param JDay Julian Day (1-366)
    ## \return Updates the fuel moistures, ERC, SC and BI of the class
    ## \remark Example: t.Update(temp, rh, prcp, ws, yearday)
    def Update(self, Year, Month, Day, Hour, Julian, Temp, MinTemp, MaxTemp, RH, MinRH, PPTAmt, pcp24, SolarRad, WS, SnowDay, RegObsHr, JDay):
        
        self.MetSet = True
        #Temp = cTof(Temp)	# Convert temperature from C to F
        Tmax = Temp
        Tmin = Temp
        RHMax = RH
        RHMin = RH
        PDur = self.PDurFunc(PPTAmt) 	# Requires that the beta parameter be set during class initialization
        SOW = 0
        if(PPTAmt > 0):
            SOW = 6

        # Assume SI units for this function and make the necessary conversions
        #self.oneten(Temp,RH,SOW)
        #self.hundredthous(Tmax,RHMax, Tmin, RHMin,JDay,PDur)
        #self.live(fToc(Temp),RH,JDay)
        #self.CalcIndexes(WS)
       
        if(Hour == RegObsHr or self.num_updates == 0):
            if (SnowDay == "y"):
                self.SnowCovered = True
            else:
                self.SnowCovered = False
        
        temp = (Temp - 32.0) * 5.0 / 9.0
        rh = RH / 100.0
        sr = SolarRad
        pptamnt = PPTAmt * 2.54
        
        if (self.num_updates == 0):
            self.OneHourFM.setMoisture(0.2)
            self.TenHourFM.setMoisture(0.2)
            self.HundredHourFM.setMoisture(0.2)
            self.ThousandHourFM.setMoisture(0.2)
        MyMC1 = 0
        MyMC10 = 0
        MyMC100 = 0
        MyMC1000 = 0
        
        neltemp = temp
        nelrh = rh
        nelsr = sr
        nelppt = pptamnt
        #if (self.SnowCovered):
        #    neltemp = 0.
        #    nelrh = 0.999
        #    nelsr = 0.
        #    nelppt = 0.
        
        self.OneHourFM.updateAll(Year, Month, Day, Hour, 0, 0, neltemp,nelrh,nelsr,nelppt,0.02179999999)
        MC1 = MyMC1 = self.OneHourFM.medianRadialMoisture() * 100
        self.TenHourFM.updateAll(Year, Month, Day, Hour, 0, 0, neltemp, nelrh, nelsr,nelppt,0.02179999999)
        MC10 = MyMC10 = self.TenHourFM.medianRadialMoisture() * 100
        self.HundredHourFM.updateAll(Year, Month, Day, Hour, 0, 0, neltemp, nelrh, nelsr, nelppt, 0.02179999999)
        MC100 = MyMC100 = self.HundredHourFM.medianRadialMoisture() * 100
        self.ThousandHourFM.updateAll(Year, Month, Day, Hour, 0, 0, neltemp, nelrh, nelsr, nelppt, 0.02179999999)
        MC1000 = MyMC1000 = self.ThousandHourFM.medianRadialMoisture() * 100
        FuelTemperature = self.OneHourFM.surfaceTemperature()
        # Update live fuel moisture once per day
        if (Hour == RegObsHr):
            if (SnowDay):
                self.nConsectiveSnowDays =self.nConsectiveSnowDays + 1
            else:
                self.nConsectiveSnowDays = 0
            self.GsiFM.Update(Temp, MaxTemp, MinTemp, RH, MinRH, Julian)
            self.HerbFM.Update(Temp, MaxTemp, MinTemp, RH, MinRH, Julian)
            self.WoodFM.Update(Temp, MaxTemp, MinTemp, RH, MinRH, Julian)
            self.m_GSI = GsiFM.CalcRunningAvgGSI()
            if nConsectiveSnowDays >= SNOWDAYS_TRIGGER:
                SnowTrigger = True
                
            else:
                SnowTrigger = False
            MCHERB = HerbFM.GetMoisture(SnowTrigger)
            MCWOOD = WoodFM.GetMoisture(SnowTrigger)
        
        # Finishing wiring in KBDI calc here
        
        #self.KBDI = CalcKBDI(pcp24, MaxTemp)

        

		
       # SetFuelMoistures(self.MC1, self.MC10, self.MC100, self.MC1000, self.MCWOOD,self.MCHERB, FuelTemperature);
	#iCalcIndexes((int)WS, SlopeClass, &fSC, &fERC, &fBI, &fIC);
        self.num_updates = self.num_updates + 1
        YesterdayJDay = Julian





class WX:

    vtmax = []
    vtmin = []
    vrhmax = []
    vrhmin = []
    vobssr = []
    vprcpdur = []
    vobsdate = []
    vobstime = []
    vobsdailydate = []
    vtemp = []
    vrh = []
    vpamnt = []
    vsrad = []
    vwindsp = []
    vsnow = []
    # Read the data in FW13 format
    def read(self,infile):
        fin = open(infile)
        num_lines = sum(1 for line in open(infile))
        for line in fin.readlines():
            line = line.strip()
            if len(line) >= 74:
                statnum = line[3:9]
                obsdate = line[9:17]
                obstime = line[17:21]
                obstype = line[21:22]
                sow = line[22:23]
                temp = line[23:26]
                atmoist = line[26:29]
                winddir = line[29:32]
                windsp = line[32:35]
                tenhr = line[35:37]
                tmax = line[37:40]
                tmin = line[40:43]
                maxrh = line[43:46]
                minrh = line[46:49]
                prcpdur = line[49:51]
                prcpamnt = line[51:56]
                wetflag = line[56:57]
                herbgfact = line[57:59]
                shrubgfact = line[59:61]
                moisttype = line[61:62]
                meastype = line[62:63]
                seascode = line[63:64]
                solrad = line[64:68]
                winddirgust = line[68:71]
                windspgust = line[71:74]
                snowflag = line[74:75]
                #if (obstime == "1300"):
                if(temp != "" and atmoist != ""):
                    self.vtmax.append(float(ValidateText(tmax)))
                    self.vtmin.append(float(ValidateText(tmin)))
                    self.vrhmax.append(float(ValidateText(maxrh)))
                    self.vrhmin.append(float(ValidateText(minrh)))
                    self.vprcpdur.append(int(ValidateText(prcpdur)))
                    self.vobsdailydate.append(obsdate)
                    self.vobssr.append(float(ValidateText(solrad)))
                    # Store all the values here
                    self.vtemp.append(float(ValidateText(temp)))
                    self.vrh.append(float(ValidateText(atmoist)))
                    
                    self.vpamnt.append(float(ValidateText(prcpamnt)))
                    
                    self.vsrad.append(float(ValidateText(solrad)))
                    self.vwindsp.append(float(ValidateText(windsp)))
                    self.vobsdate.append(obsdate)
                    self.vobstime.append(obstime)
                    self.vsnow.append(snowflag)
           #print (statnum,obsdate,obstime,obstype,sow,temp,atmoist,winddir,windsp,  tenhr,tmax,tmin,maxrh,minrh,prcpdur,prcpamnt,wetflag,herbgfact,shrubgfact,moisttype,meastype,seascode,solrad,winddirgust,windspgust,snowflag)
        return


def ValidateText(txtin):

    if(len(txtin.strip()) > 0):
        return txtin
    else:
        return "-999.9"

def FToC(F):
    return ((F - 32) * 5/9)


def CalcDaylightHours( Julian,Lat):
    phi = 0.0
    xfact = 0.0
    decl = 0.0
    tla = 0.0
    phi = tan(Lat * 0.01745) * -1.0
    xfact = (Julian - 80) * 0.01745
    decl = 23.5 * sin(xfact)
    decl = decl * 0.01745
    tla = phi * sin(decl)

    if (fabs(tla) < .01):
        tla = 0.01
    else:
        if (tla >= 1.0):
            tla = 0.99999999
        else:
            if (tla <= -1.0):
                tla = -.9999999

    tla = atan(sqrt((1.0 - tla * tla)) / tla)
    if (tla < 0.0):
        tla = tla + 3.141593
    return tla * 7.64

def GetHourlyTemp(TMax,TMin, Dayl, H):

    nighthours = 24 - Dayl
    tavg = (TMax + TMin) / 2
    if (H < nighthours):
        return TMin

    else:
        DayHour = (H - nighthours) - 1
    if(DayHour < Dayl * 0.25):
        MyH = (DayHour) + 1
        return (MyH * (tavg - TMin) / (Dayl / 4)) + TMin
    if(DayHour < Dayl * 0.5):
        MyH = DayHour - Dayl / 4 + 1
        return ((MyH * (TMax - tavg)) / (Dayl / 4)) + tavg
    if(DayHour < (Dayl * 0.75)):
        MyH = DayHour - Dayl / 2 + 1
        return -((MyH * (TMax - tavg)) / (Dayl / 4)) + TMax
    if (DayHour >= (Dayl * 0.75)):
        MyH = DayHour - Dayl * 3 / 4 + 1
        return -((MyH * (tavg - TMin)) / (Dayl / 4)) + tavg


if __name__ == "__main__":
    print("Loaded")
    #fdr = USNFDRSCalc(45,"Y",0.188503)
    #fdr.SetMoistures(3,4,5,6,60,90)
    #fdr.CalcIndexes(10,100,80)
    #print(fdr.ERC)
    