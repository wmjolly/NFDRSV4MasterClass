# -*- coding: utf-8 -*-
"""
    \file LiveFuelMoisture.py
    \brief LiveFuelMoisture class definition and implementation.
    \author Copyright (C) 2015 by W. Matt Jolly, USFS, RMRS, Fire Sciences Laboratory
    \version 1.0.0 - Uses only standard Python

    \par Licensed under GNU GPL
    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    If you have not received a copy of the GNU General Public License
    along with this program; write to the Free Software Foundation, Inc.,
    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
    The GNU GPL is also available in full from
    http://www.opensource.org/licenses/gpl-license.php

Created on Mon Jun 17 11:13:25 2019

@author: mjolly
"""
from math import *
import matplotlib.pyplot as plt
import pandas as pd
from datetime import datetime


#define NOVALUE -9999.9
#define RADPERDAY 0.017214
#define RADPERDEG 0.01745329
#define MINDECL -0.4092797
#define SECPERRAD 13750.9871
#define DAYSOFF 10.25

    
class LiveFuelMoisture:
    m_Lat = 0.0
    m_IsHerb = False
    m_IsAnnual = False
    m_MaxGSI = 1.0
    m_GreenupThreshold=0.5
    m_MinLFMVal=30
    m_MaxLFMVal=250
    m_TminMin = -2.0
    m_TminMax = 5.0
    m_VPDMin = 900
    m_VPDMax = 4100
    m_DaylenMin=36000
    m_DaylenMax=39600
    iGSI = []
    m_UseVPDAvg = False
    m_LFIdaysAvg = 21
    hasGreenedUpThisYear = True
    canIncreaseHerb = True
    hasExceeded120ThisYear = False
    lastHerbFM = -1.0
    NOVALUE = -9999.9
    RADPERDAY = 0.017214
    RADPERDEG = 0.01745329
    MINDECL = -0.4092797
    SECPERRAD = 13750.9871
    DAYSOFF = 10.25
    
    def __init__(self,Lat,IsHerb=False,IsAnnual=False):    
        self.m_Lat = Lat
        self.m_IsHerb = IsHerb
        self.m_IsAnnual = IsAnnual
        
        if(self.m_IsHerb):
            self.SetLFMParameters(1.0,0.5,30,250)
        else:
            self.SetLFMParameters(1.0,0.5,60,250)
        self.SetLimits()
        self.SetMAPeriod()
        # Parameters to control green-up and curing of annual herbs
        self.m_IsAnnual = True # Overriden when the object is initialized
        self.ResetHerbState()
        self.m_UseVPDAvg = True
 
        if (len(self.iGSI) > 0):
            self.iGSI = []
            
    def SetLFMParameters(self,MaxGSI = 1.0,GreenupThreshold =0.5,MinLFMVal=30,MaxLFMVal=250):
        self.m_MaxGSI = MaxGSI
        self.m_GreenupThreshold = GreenupThreshold
        if(self.m_GreenupThreshold == 1.0):
            self.m_GreenupThreshold = 0.9999
        if(self.m_IsHerb):
            self.m_MaxLFMVal = MaxLFMVal
            self.m_MinLFMVal = MinLFMVal
            self.m_Slope = (self.m_MaxLFMVal - self.m_MinLFMVal) / (1.0 - self.m_GreenupThreshold)
            self.m_Intercept = self.m_MaxLFMVal - self.m_Slope
        else:
    
            self.m_MaxLFMVal = MaxLFMVal
            self.m_MinLFMVal = MinLFMVal
            self.m_Slope = (self.m_MaxLFMVal - self.m_MinLFMVal) / (1.0 - self.m_GreenupThreshold)
            self.m_Intercept = self.m_MaxLFMVal - self.m_Slope


    #/*! \brief Method to retrieve calculated live fuel moisture for a given instance of the class.
    #    \param[in] NONE.
    #    \return Calculated live fuel moisture (% dry wt)
     
    def GetMoisture(self,SnowDay):
    
         if (self.m_IsHerb):
             return self.CalcRunningAvgHerbFM(SnowDay)
         else:
              return self.CalcRunningAvgWoodyFM(SnowDay)
     

    def Update(self,TempF, MaxTempF, MinTempF, RH, MinRH, Jday):
        GSI = 0.0
        if not self.m_UseVPDAvg:
            GSI = self.CalcGSI(MinRH, MaxTempF, MinTempF, self.m_Lat, Jday)
        else:
            GSI = self.CalcGSI_VPDAvg(RH, TempF, MaxTempF, MinTempF, self.m_Lat, Jday)
	
        self.iGSI.append(GSI)

    #/*! \brief Update function for GSI-based fuel moistures based on daily weather data.
    #    Note: Latitude is set during initialization.
    #    \param[in] RH: Relative Humidity (1-100%).
    #    \param[in] TempF: Air temperature (F).
    #    \param[in] MinTempF: 24-hour minimum temperature (F)
    #    \param[in] Jday: Day of Year (1-366)
    #    \return NONE
    def UpdateSimple(self,RH, TempF, MinTempF,Jday):
        GSI = self.CalcGSI(RH, TempF, MinTempF, self.m_Lat, Jday)
        self.iGSI.push_back(GSI)



    def CalcGSI(self,minRH, maxTempF, minTempF, lat, doy):
        tMinInd = self.GetTminInd(minTempF)
        vpd = self.CalcVPD(max(minRH, 5.0), maxTempF)
        vpdInd = self.GetVPDInd(vpd)
        daylenInd = self.GetDaylInd(self.CalcDayl(lat, doy))
        GSI = tMinInd * vpdInd * daylenInd
        return GSI
    
    def CalcGSI_VPDAvg(self,RH, TempF, maxTempF, minTempF, lat, doy):
        tMinInd = self.GetTminInd(minTempF)
        tDew = self.CalcDPT(TempF, RH)
        vpd = self.CalcVPDavg(tDew, (maxTempF + minTempF) / 2)
        vpdInd = self.GetVPDInd(vpd)
        daylenInd = self.GetDaylInd(self.CalcDayl(lat, doy))
        GSI = tMinInd * vpdInd * daylenInd;
        return GSI

    def GetTminInd(self,Tmin):
        tmin =  (Tmin - 32.0) * 5.0 / 9.0; # Convert Tmin from Fahrenheit to celcuius
        if(self.m_TminMax == self.m_TminMin):
            return 0
        if( tmin < self.m_TminMin):
            return 0
        elif(tmin > self.m_TminMax):
            return 1
        else:
            return (tmin - self.m_TminMin) / (self.m_TminMax - self.m_TminMin)
    
    def GetVPDInd(self,VPD):
        if(self.m_VPDMax == self.m_VPDMin):
            return 0
        if( VPD < self.m_VPDMin):
            return 1
        elif(VPD > self.m_VPDMax):
            return 0
        else:
            return 1 - (VPD - self.m_VPDMin) / (self.m_VPDMax - self.m_VPDMin)

    def GetDaylInd(self,Dayl):
        if(self.m_DaylenMin == self.m_DaylenMax):
            return 0
        if( Dayl < self.m_DaylenMin):
            return 0
        elif(Dayl > self.m_DaylenMax):
            return 1
        else:
            return (Dayl - self.m_DaylenMin) / (self.m_DaylenMax - self.m_DaylenMin)



    #/*! \brief Set the environmental limits for the GSI calculations
    #    This function defines the environmental limits for the minimum temperature,
    #    vapor pressure deficit and daylength indicator functions used to calculate
    #    GSI.
    
    #    \param[in] TminMin: Lower limit for minimum temperature (C).
    #    \param[in] TminMax: Upper limit for minimum temperature (C).
    #    \param[in] VPDMin: Lower limit for vapor pressure deficit (Pa).
    #    \param[in] VPDMax: Upper limit for vapor pressure deficit (Pa).
    #    \param[in] DaylMin: Lower limit for daylength (minutes)
    #    \param[in] DaylMax: Upper limit for daylength (minutes)
    #    \return NONE
    # */
    def SetLimits(self,TminMin = -2.0,TminMax = 5.0,VPDMin = 900,VPDMax = 4100,DaylMin = 36000,DaylMax = 39600):
        self.m_TminMin = TminMin
        self.m_TminMax = TminMax
        self.m_VPDMin = VPDMin
        self.m_VPDMax = VPDMax
        self.m_DaylenMin = DaylMin
        self.m_DaylenMax = DaylMax
    def SetMAPeriod(self,MAPeriod=21):
        self.m_LFIdaysAvg = max(1, MAPeriod)
    def SetUseVPDAvg(self,Flag):
        self.m_UseVPDAvg = Flag
    def GetUseVPDAvg(self):
        return self.m_UseVPDAvg
    def GetIsAnnual(self):
        return self.m_IsAnnual
    def GetLFMParameters(self):
        return [self.m_MaxGSI,self.m_GreenupThrehold,self.m_MinLFMVal,self.m_MaxLFMVal]
    def GetMaxGSI(self):
        return self.m_MaxGSI
    def GetGreenupThreshold(self):
        return self.m_GreenupThreshold
    def GetMinLFMVal(self):
        return self.m_MinLFMVal    
    def GetMaxLFMVal(self):
        return self.m_MaxLFMVal
        
    def CalcRunningAvgGSI(self):
        numValid = 0
        days_start = 0
        days_end = 0
        gsi = 0.0

        if(len(self.iGSI) < self.m_LFIdaysAvg):
            days_start = 0
            days_end = len(self.iGSI)
        else:
            days_start = len(self.iGSI) - self.m_LFIdaysAvg
            days_end = len(self.iGSI)
        for i in range(days_start,days_end):
            if self.iGSI[i] >= 0:
                numValid = numValid + 1
                gsi = gsi + self.iGSI[i]
        
        if numValid > 0:
            gsi /= numValid
            
        return gsi
    def CalcRunningAvgHerbFM(self,SnowDay):
        GSI = self.CalcRunningAvgGSI()
        rescale = GSI / self.m_MaxGSI
        ret = self.m_MinLFMVal
        rescale = min(1.0, rescale)
        rescale = max(0.0, rescale)
        if rescale >= self.m_GreenupThreshold and not SnowDay == 1:
            ret = self.m_Slope * rescale + self.m_Intercept
            if not self.hasGreenedUpThisYear:
                self.hasGreenedUpThisYear = True;
                self.canIncreaseHerb = True;

        if not self.canIncreaseHerb and self.lastHerbFM >= 0:
            ret = min(ret, self.lastHerbFM)
        if not self.hasExceeded120ThisYear and ret >= 120:
            self.hasExceeded120ThisYear = True
        if self.hasExceeded120ThisYear and ret < 120.0 and self.m_IsAnnual:
            self.canIncreaseHerb = False
        
        self.lastHerbFM = ret
        return ret

    def ResetHerbState(self):
        self.hasGreenedUpThisYear = False
        self.canIncreaseHerb = False
        self.hasExceeded120ThisYear = False
        self.lastHerbFM = -1.0

    def CalcRunningAvgWoodyFM(self,SnowDay):
        GSI = self.CalcRunningAvgGSI()
        rescale = GSI / self.m_MaxGSI
        rescale = min(1.0, rescale)
        rescale = max(0.0, rescale)
        if(rescale >= self.m_GreenupThreshold and not SnowDay == 1):
            return self.m_Slope * rescale + self.m_Intercept
        return self.m_MinLFMVal

    def CalcDPT(self,tempF, RH):
        safeTemp = min(140.0, tempF)
        safeRH = max(5.0, RH)
        return -398.36 - 7428.6 / (-15.674 + log(safeRH / 100.0 * exp(-7482.6/(safeTemp + 398.36) + 15.675)))
        
    def CalcVP(self,tempF):
        tmpC =  (tempF - 32.0) / 1.8
        vp = 610.7 * exp((17.38 * tmpC)/(239 + tmpC))
        return vp
    def CalcVPD(self,RH, TempF):
        vp = self.CalcVP(TempF)
        vpd = vp - (RH / 100) * vp
        if(vpd < 0.0):
            vpd = 0.0;
        return vpd
        
    def CalcVPDavg(self,TempDewF,TempAvgF):
        vpDew = self.CalcVP(TempDewF)
        vpAvg = self.CalcVP(TempAvgF)
        vpd = vpAvg - vpDew
        if(vpd < 0.0):
            vpd = 0.0;
        return vpd


    def CalcDayl(self,lat,yday):
        # Daylength function from MT-CLIM */
        lat = lat * self.RADPERDEG
        if lat > 1.5707:
            lat = 1.5707
        if lat < -1.5707:
            lat = -1.5707
        coslat = cos(lat)
        sinlat = sin(lat)

        #* calculate cos and sin of declination */
        decl = self.MINDECL * cos((yday + self.DAYSOFF) * self.RADPERDAY)
        cosdecl = cos(decl)
        sindecl = sin(decl)
        cosegeom = coslat * cosdecl
        sinegeom = sinlat * sindecl
        coshss = -(sinegeom) / cosegeom
        if coshss < -1.0:
            coshss = -1.0  # 24-hr daylight */
        if coshss > 1.0:
            coshss = 1.0    # 0-hr daylight */
        hss = acos(coshss)                # hour angle at sunset (radians) */
        #* daylength (seconds) */
        return 2.0 * hss * self.SECPERRAD


if __name__ == "__main__":
    l = LiveFuelMoisture(45,False,False)
    print(l.CalcDayl(45,20))
     
   
    filen = r'S:\NFDRS\pyNFDRS\MSOShort.csv'
    dat = pd.read_csv(filen,sep='\s*,\s*',engine='python')  # doctest: +SKIP  
    
    datdt = []
    temp = []
    import csv
    with open(filen, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=',', quotechar='|')
        i = 0        
        for row in reader:
            if i > 0:
                dstr = "%s %s" % (row[0],row[1])
                dt = datetime.strptime(dstr, '%m/%d/%Y %H%M')
                print(dt)
                #DATE,Time,Temp, MinT, MaxT,   RH, MnRH, MxRH, Rain, RnDr, HrRa, Wind,
                t = float(row[2])
                if t < 140:
                    if(row[1] == "1300"):
                        JDay = float(dt.strftime("%j"))
                        print(JDay)
                        Temp = float(row[2])
                        MinT = float(row[3])
                        MaxT = float(row[4])
                        RH = float(row[5])
                        MinRH = float(row[6])
                        MaxRH = float(row[7])
                        Rain = float(row[10])
                        l.Update(Temp,MaxT,MinT,RH,MinRH,JDay)
                        temp.append(l.GetMoisture(False))
                        
                        datdt.append(dt)
                    

                
                
            i = i + 1
    plt.plot_date(datdt, temp,linestyle='solid', marker='None')
    plt.xlabel('Date')
    plt.ylabel('Temperature (Deg C)')
    plt.show()
            
    #dat.groupby('DATE')['MaxT'].mean().plot()
    #for row in dat:
    #    print (row)
        #dat.loc[0,'MaxT']
    #plt.scatter(df['DATE'], df['MaxT'])
    #plt.show() # Dependi
#    plt.plot(ind,val)
 #   plt.xlabel('Date')
  #  plt.ylabel('Daylength')
   # plt.show()
    