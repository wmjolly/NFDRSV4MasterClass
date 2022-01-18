__author__ = 'cbevins & mjolly'

from math import *
import math
from datetime import datetime, timedelta, time
import random
import numpy as np

# ------------------------------------------------------------------------------
#/*! \file NelsonDeadFuelMoisture.py
#    \brief DeadFuelMoisture class definition and implementation.
#    \author Copyright (C) 2005 by Collin D. Bevins.
#     \author W. Matt Jolly, ported to Python from CBevins code
#   \version 1.0.0 - Uses only standard Python libraries

#    \par Licensed under GNU GPL
#    This program is free software you can redistribute it and/or modify it
#    under the terms of the GNU General Public License as published by the
#    Free Software Foundation either version 2 of the License, or (at your
#    option) any later version.#

#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.#

#    If you have not received a copy of the GNU General Public License
#    along with this program write to the Free Software Foundation, Inc.,
#    51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#    The GNU GPL is also available in full from
#    http:#www.opensource.org/licenses/gpl-license.php
# */

Aks = 2.0e-13
Alpha = 0.25    #/*! \var Alpha \brief Fraction of cell length that overlaps adjacent cells (0.25 cm/cm).
Ap = 0.000772   #/*! \var Ap \brief Psychrometric constant (0.000772 / oC).
Aw = 0.8        #/*! \var Aw \brief Ratio of cell cavity to total cell width (0.8 cm/cm).
Eps = 0.85      #/*! \var Eps \brief Longwave emissivity of stick surface (0.85 dl).

#/*! \var Hfs
#    \brief Saturation value of the stick surface humidity (0.99 g/g).
Hfs = 0.99

#/*! \var Kelvin
#    \brief Celcius-to-Kelvin offset (273.2 oC).
Kelvin = 273.2

#/*! \var Pi
#    \brief A well-rounded number (dl).
Pi = 3.141592654

#/*! \var Pr
#    \brief Prandtl number (0.7 dl).
Pr = 0.7

#/*! \var Sbc
#    \brief Stefan-Boltzmann constant (1.37e-12 cal/cm2-s-K4).
Sbc = 1.37e-12

#/*! \var Sc
#    \brief Schmidt number (0.58 dl).
Sc = 0.58

#/*! \var Smv
#    \brief Factor to convert solar radiation from W/m2 to milliVolts
#    \f$mv = \frac {W / m^{2}} {94.743}\f$
Smv = 94.743

#/*! \var St
#\brief Surface tension (72.8).
St = 72.8

#/*! \var Tcd
#    \brief Day time clear sky temperature (6 oC).
Tcd = 6.

#/*! \var Tcn
#    \brief Night time clear sky temperature (3 oC).
Tcn = 3.

#/*! \var Thdiff
#    \brief Thermal diffusivity ( 8.0 cms/h).
Thdiff = 8.0

#/*! \var Wl
#\brief Diameter of interior cell cavity ( 0.0023 cm).
Wl = 0.0023

#/*! \var Srf
#    \brief Factor to derive "sr", solar radiation received (14.82052 cal/cm2-h).
Srf = 14.82052

#/*! \var Wsf
#    \brief Manifest constant equal to -log(1.0 - 0.99)
Wsf = 4.60517

#/*! \var Hrd
#    \brief Factor to derive daytime long wave radiative surface heat transfer.
#    \arg Let Sb = 4. * 3600. * Sbc * Eps = 1.67688e-08
#    \arg Let tsk = Tcd + Kelvin
#    \arg Then Hrd = Sb * tsk * tsk * tsk / Pi
#    \arg And Hrd = 0.116171
Hrd = 0.116171

#/*! \var Hrn
#    \brief Factor to derive nighttime long wave radiative surface heat transfer.
#    \arg Let Sb = 4. * 3600. * Sbc * Eps = 1.67688e-08
#    \arg Let tsk = Tcn + Kelvin
#    \arg Then Hrn = Sb * tsk * tsk * tsk / Pi
#    \arg And Hrn = 0.112467
Hrn = 0.112467

#/*! \var Sir
#    \brief Saturation value below which liquid water columns no longer exist.
#    \arg Sir = Aw * Alpha / (4. * (1.-(2.-Aw) * Alpha))
Sir = 0.0714285

#/*! \var Scr
#    \brief Saturation value at which liquid miniscus first enters the tapered
#    portion of wood cells.
#    \arg Scr = 4. * Sir = 0.285714
Scr = 0.285714

class DeadFuelMoisture:
    DFM_State_None = 0
    DFM_State_Adsorption = 1
    DFM_State_Desorption = 2
    DFM_State_Condensation1 = 3
    DFM_State_Condensation2 = 4
    DFM_State_Evaporation = 5
    DFM_State_Rainfall1 = 6
    DFM_State_Rainfall2 = 7
    DFM_State_Rainstorm = 8
    DFM_State_Stagnation = 9
    DFM_State_Error = 10

    m_density = 0.0  #!< Stick density (g/cm3).
    m_dSteps = 0  #!< Number of diffusivity computation steps per observation.
    m_hc = 0.0  #!< Stick planar heat transfer rate (cal/cm2-h-C).
    m_length = 0.0  #!< Stick length (cm).
    m_name = ""  #!< Stick name or other descriptive text.
    m_nodes = 0  #!< Number of stick nodes in the radial direction.
    m_radius = 0.0  #!< Stick radius (cm).
    m_rai0 = 0  #!< Rain runoff factor during the initial hour of rainfall (dl).
    m_rai1 = 0  #!< Rain runoff factor after the initial hour of rainfall (dl) [no longer used].
    m_stca = 0.0  #!< Adsorption surface mass transfer rate ((cm3/cm2)/h).
    m_stcd = 0.0  #!< Desorption surface mass transfer rate ((cm3/cm2)/h).
    m_mSteps = 0  #!< Number of moisture content computation steps per observation.
    m_stv = 0.0  #!< Storm transition value (cm/h) [no longer used].
    m_wfilmk = 0.0  #!< Water film contribution to stick moisture content (g water/g dry fuel).
    m_wmx = 0.0  #!< Stick maximum local moisture due to rain (g water/g dry fuel).
    m_init = False
    # Configuration parameters
    m_allowRainfall2 = False  # If TRUE, applies Nelson's logic for rainfall runoff after the first hour
    m_allowRainstorm = False  # If TRUE, applies Nelson's logic for rainstorm transition and state
    m_pertubateColumn = True  # If TRUE, the continuous liquid column condition get pertubated
    m_rampRai0 = True  # If TRUE, used Bevins' ramping of rainfall runoff factor rather than Nelsons rai0 *= 0.15
    m_aml = 0

    # Intermediate stick variables derived in initializeStick()
    m_dx = 0.0  #!< Internodal radial distance (cm).
    m_wmax = 0.0  #!< Maximum possible stick moisture content (g water/g dry fuel).
    m_x = []  #!< Array of nodal radial distances from stick center (cm).
    m_v = []  #!< Array of nodal volume weighting fractions (cm3 node/cm3 stick).

    # Optimization factors derived in initializeStick()
    m_amlf = 0.0  #!< \a aml optimization factor.
    m_capf = 0.0  #!< \a cap optimization factor.
    m_hwf = 0.0  #!< \a hw and \a aml computation factor.
    m_dx_2 = 0.0  #!< 2 times the internodal distance \a m_dx (cm).
    m_vf = 0.0  #!< optimization factor used in update().

    # Environmental variables provided to initializeEnvironment():
    m_bp0 = 0.0  #!< Previous observation's barometric presure (cal/cm3).
    m_ha0 = 0.0  #!< Previous observation's air humidity (dl).
    m_rc0 = 0.0  #!< Previous observation's cumulative rainfall amount (cm).
    m_sv0 = 0.0  #!< Previous observation's solar radiation (mV).
    m_ta0 = 0.0  #!< Previous observation's air temperature (oC).
    m_init = 0.0  #!< Flag set by initialize().

    # Environmental variables provided to update():
    m_bp1 = 0.0  #!< Current observation's barometric pressure (cal/cm3).
    m_et = 0.0  #!< Elapsed time since previous observation (h).
    m_ha1 = 0.0  #!< Current observation's air humidity (dl).
    m_rc1 = 0.0  #!< Current observation's cumulative rainfall amount (cm).
    m_sv1 = 0.0  #!< Current observation's solar radiation (mV).
    m_ta1 = 0.0  #!< Current observation's air temperature (oC).

    # Intermediate environmental variables derived by update():
    m_ddt = 0.0  #!< Stick diffusivity computation interval (h).
    m_mdt = 0.0  #!< Stick moisture content computation interval (h).
    m_mdt_2 = 0.0  #!< 2 times the moisture time step \a m_mdt (h).
    m_pptrate = 0.0  #!< Current observation period's rainfall rate (cm/h).
    m_ra0 = 0.0  #!< Previous observation period's rainfall amount (cm).
    m_ra1 = 0.0  #!< Current observation period's rainfall amount (cm).
    m_rdur = 0.0  #!< Current rainfall event duration (h)
    m_sf = 0.0  #!< optimization factor used in update().
    m_rcum = 0.0 # Cumulative rainfall
    # Stick moisture condition variables derived in update():
    m_hf = 0.0  #!< Stick surface humidity (g water/g dry fuel).
    m_wsa = 0.0  #!< Stick fiber saturation point (g water/g dry fuel).
    m_sem = 0.0  #!< Stick equilibrium moisture content (g water/g dry fuel).
    m_wfilm = 0.0  #!< Amount of water film (0 or \a m_wfilmk) (g water/g dry fuel).
    m_elapsed = 0.0  #!< Total simulation elapsed time (h).
    m_t = []  #!< Array of nodal temperatures (oC).
    m_s = []  #!< Array of nodal fiber saturation points (g water/g dry fuel).
    m_d = []  #!< Array of nodal bound water diffusivities (cm2/h).
    m_w = []  #!< Array of nodal moisture contents (g water/g dry fuel).
    m_updates = 0  #!< Number of calls made to update().
    m_state = ""  #!< Prevailing dead fuel moisture state.
    m_randseed = 0  #!< If not zero, nodal temperature, saturation, and moisture contents are pertubated by some small amount. If < 0, uses system clock for seed.
    m_Ttold = []  #!< Temporary array of nodal temperatures (oC).
    m_Tsold = []  #!< Temporary array of nodal fiber saturation points (g water/g dry fuel).
    m_Twold = []  #!< Temporary array of nodal moisture contents (g water/g dry fuel).
    m_Tv = []  #!< Temporary array used to redistribute nodal temperatures
    m_To = []  #!< Temporary array used to redistribute moisture contents
    m_Tg = []  #!< Temporary array of nodal free water transport coefficients
    m_datetime = datetime.now()



    #------------------------------------------------------------------------------
    #! \brief Default class constructor.

    #    Creates a dead fuel moisture stick with the passed \a radius and \a name,
    #    and derives all other parameters using interpolations from \ref bevins2005.

    #\param[in] radius Dead fuel stick radius (cm).
    #\param[in] name Name or description of the dead fuel stick.

    def __init__(self, radius, name):
        self.initializeParameters(radius, name)
        # 1-h time lag 0.2cm (.16in)
        # 10-h time lag 0.64 (0.5in)
        # 100-h time lag 2.0 cm (1.57in)
        # 1000-h time lag 6.4 cm (5in)

    #! \brief Access to the stick's adsorption surface mass transfer rate.

    #    \return Current surface mass transfer rate for adsorption
    #    (\f$\frac {cm^{3}} {cm^{2} \cdot h}\f$).
    def adsorptionRate(self):
        return self.m_stca


    #------------------------------------------------------------------------------
    # \brief Static convenience method to determine the adsorption surface mass
    #    transfer rate for a DeadFuelMoisture stick of given \a radius.

    #    The adsorption surface mass transfer rate is determined from
    #    \ref bevins2005
    #    \f[\alpha = 0.0004509 + \frac {0.006126}{r^2.6}\f]
    #    where
    #    - \f$\alpha\f$ is the adsorption surface mass transfer rate
    #        (\f$\frac {cm^{3}} {cm^{2} \cdot h}\f$), and
    #    - \f$r\f$ is the stick radius (\f$cm\f$).

    #    \param[in] radius Dead fuel stick radius (\f$cm\f$).

    #    \return Estimated surface mass transfer rate for adsorption
    #    (\f$\frac {cm^{3}} {cm^{2} \cdot h}\f$).

    def deriveAdsorptionRate(self, radius):
        alpha = 0.06 + 0.006126 / pow(radius, 2.6)
        return alpha

    # \brief Static convenience method to determine the minimum number of
    #    diffusivity computation time steps per observation for a DeadFuelMoisture stick of given \a radius.

    #    The minimum number of diffusivity computation time steps per observation
    #    must be large enough to provide stability in the model computations
    #    (\ref carlson2004a ), and is used to determine the internal diffusivity
    #    computation time step within update().
    #
    #   The minumum number of diffusivity computation time steps is determined from
    #   \ref bevins2005
    #    \f[$n_{d} = int( 4.777 + \frac {2.496}{r^1.3} )\f]
    #    where
    #    - \f$n_{d}\f$ is the minimum number of diffusivity computation time steps
    #    per observation, and
    #    - \f$r\f$ is the stick radius (\f$cm\f$).

    #    \param[in] radius Dead fuel stick radius (\f$cm\f$).

    #   \return Minimum number of diffusivity computation time steps per observation.

    def deriveDiffusivitySteps(self, radius):
        steps = int(4.777 + 2.496 / pow(radius, 1.3))
        return steps

    # \brief Static convenience method to determine the minimum number of
    #    moisture content computation time steps per observation
    #    for a DeadFuelMoisture stick of given \a radius.

    #   The number of moisture computation time steps per observation must be
    #    large enough to provide stability in the model computations
    #    (\ref carlson2004a ), and is used to determine the internal moisture content
    #    computation time step within update().

    #    The minimum number of moisture content computation time steps
    #    is determined from
    #    \ref bevins2005 \f[$n_{m} = int( 9.8202 + \frac {26.865}{r^1.4} )\f]
    #    where
    #    - \f$n_{m}\f$ is the minimum number of moisture content computation time
    #    steps per observation, and
    #    - \f$r\f$ is the stick radius (\f$cm\f$).

    #    \param[in] radius Dead fuel stick radius (\f$cm\f$).

    #    \return Minimum number of moisture computation time steps per observation.


    def deriveMoistureSteps(self,radius):
        steps = int(9.8202 + 26.865 / pow(radius, 1.4))
        return steps


    #/*! \brief Static convenience method to determine the planar heat transfer rate
    #    for a DeadFuelMoisture stick of given \a radius.
    #The planar heat transfer rate is determined from \ref bevins2005
    #\f[$h_{c} = 0.2195 + \frac {0.05260}{r^2.5}\f]
    #where
    #- \f$h_{c}\f$ is the planar heat transfer rate
    #    (\f$\frac {cal} {cm^{2} \cdot h \cdot C}\f$), and
    #- \f$r\f$ is the stick radius (\f$cm\f$).

    #\param[in] radius Dead fuel stick radius (\f$cm\f$).

    #\return Estimated planar heat transfer rate
    #    (\f$\frac {cal} {cm^{2} \cdot h \cdot C}\f$).
    def derivePlanarHeatTransferRate(self,radius):
        hc = 0.2195 + 0.05260 / pow(radius, 2.5)
        return ( hc )


    #/*! \brief Static convenience method to determine the rainfall runoff factor
    #for a DeadFuelMoisture stick of given \a radius.

    #The initial rainfall runoff factor is determined from \ref bevins2005
    #\f[$f_{0} = 0.02822 + \frac {0.1056}{r^2.2}\f]
    #where
    #- \f$f_{0}\f$ is the initial rainfall runoff rate (\f$dl\f$), and
    #- \f$r\f$ is the stick radius (\f$cm\f$).

    #\param[in] radius Dead fuel stick radius (\f$cm\f$).

    #\return Estimated initial rainfall runoff factor (\f$dl\f$).


    def deriveRainfallRunoffFactor(self,radius):
        rrf0 = 0.02822 + 0.1056 / pow(radius, 2.2)
        return (rrf0)

    #------------------------------------------------------------------------------
    #/*! \brief Static convenience method to determine the number of moisture content
    #    radial computation nodes for a DeadFuelMoisture stick of given \a radius.

    #    The number of moisture computation nodes must be large enough to provide
    #    stability in the model computations (\ref carlson2004a ).

    #    The minimum number of moisture content computation nodes
    #    is determined from
    #    \ref bevins2005 \f[$n_{m} = int( 10.727 + \frac {0.1746}{r} )\f]
    #    where
    #    - \f$n_{m}\f$ is the minimum number of moisture content computation
    #        radial nodes,
    #    - \f$r\f$ is the stick radius (\f$cm\f$).
    #    \param[in] radius Dead fuel stick radius (\f$cm\f$).
    #    \return Minimum number of moisture computation radial nodes.
    def deriveStickNodes(self,radius):
        nodes = (int)(10.727 + 0.1746 / radius)
        if ( ( nodes % 2 ) == 0 ):
            nodes = nodes + 1
        return ( nodes )

    #------------------------------------------------------------------------------
    #/*! \brief Access to the stick's desorption surface mass transfer rate.

    #    \return Current surface mass transfer rate for desorption
    #    (\f$\frac {cm^{3}} {cm^{2} \cdot h}\f$).
    def desorptionRate(self):
        return ( self.m_stcd )


    #------------------------------------------------------------------------------
    #/*! \brief Determines bound water diffusivity at each radial nodes.
    #    \param[in] bp Barometric pressure (cal/m3)
    def diffusivity(self, bp):
        #changed SB, 10/2009, moce variable declarations outside of loop
        #cpv = 0, dv = 0, ps1 = 0, c1 = 0, c2 = 0,  wc = 0
        #dhdm = 0, daw = 0, svaw = 0, vfaw = 0, vfcw = 0, rfcw = 0, fac = 0
        #con = 0, qw = 0, e = 0, dvpr = 0

        for i in range(0, self.m_nodes):
            # Stick temperature (oK)
            tk = self.m_t[i] + 273.2
            # Latent heat of vaporization of water (cal/mol
            qv = 13550. - 10.22 * tk
            # Specific heat of water vapor (cal/(mol*K))
            #cpv   = 7.22 + .002374 * tk + 2.67e-07 * tk * tk
            cpv = 7.22 + (tk * (0.002374 + 2.67e-07 * tk))
            # Sea level atmospheric pressure = 0.0242 cal/cm3
            dv = 0.22 * 3600. * ( 0.0242 / bp ) * (( tk / 273.2 ) ** 1.75 )
            # Water saturation vapor pressure at surface temp (cal/cm3)
            ps1 = 0.0000239 * math.exp(20.58 - (5205. / tk))
            # Emc sorption isotherm parameter (g/g)
            c1 = 0.1617 - 0.001419 * self.m_t[i]
            # Emc sorption isotherm parameter (g/g)
            c2 = 0.4657 + 0.003578 * self.m_t[i]
            # Lesser of nodal or fiber saturation moisture (g/g)
            #Reciprocal slope of the sorption isotherm
            dhdm = 0.0
            if ( self.m_w[i] < self.m_wsa ):
                wc = self.m_w[i]
                if ( c2 != 1. and self.m_hf < 1.0 and c1 != 0.0 and c2 != 0.0 ):
                    dhdm = (1.0 - self.m_hf) * ((-math.log(1.0 - self.m_hf)) ** (1.0 - c2)) / (c1 * c2)
            else:
                wc = self.m_wsa
                if ( c2 != 1. and Hfs < 1.0 and c1 != 0.0 and c2 != 0.0 ):
                    dhdm = (1.0 - Hfs) * pow(Wsf, (1.0 - c2)) / (c1 * c2)

            # Density of adsorbed water (g/cm3)
            daw = 1.3 - 0.64 * wc
            # Specific volume of adsorbed water (cm3/g)
            svaw = 1. / daw
            # Volume fraction of adborbed water (dl)
            vfaw = svaw * wc / (0.685 + svaw * wc)
            # Volume fraction of moist cell wall (dl)
            vfcw = (0.685 + svaw * wc) / ((1.0 / self.m_density) + svaw * wc)
            # Converts D from wood substance to whole wood basis
            rfcw = 1.0 - math.sqrt(1.0 - vfcw)
            # Converts D from wood substance to whole wood basis
            fac = 1.0 / (rfcw * vfcw)
            # Correction for tortuous paths in cell wall
            con = 1.0 / (2.0 - vfaw)
            # Differential heat of sorption of water (cal/mol)
            qw = 5040. * math.exp(-14.0 * wc)
            # Activation energy for bound water diffusion (cal/mol)
            e = (qv + qw - cpv * tk) / 1.2

            #----------------------------------------------------------------------
            # The factor 0.016 below is a correction for hindered water vapor
            # diffusion (this is 62.5 times smaller than the bulk vapor diffusion)
            #  0.0242 cal/cm3 = sea level atm pressure
            #      -- note from Ralph Nelson
            #----------------------------------------------------------------------

            dvpr = 18.0 * 0.016 * (1.0 - vfcw) * dv * ps1 * dhdm / ( self.m_density * 1.987 * tk )
            self.m_d[i] = dvpr + 3600. * 0.0985 * con * fac * math.exp(-e / (1.987 * tk))

    #/*! \brief Access to the stick's number of moisture diffusivity computation
    #    time steps per observation.

    #    \return Current number of moisture diffusivity computation time steps
    #    per observation.
    def diffusivitySteps(self):
        return self.m_dSteps

    #/*! \brief Access to the current total running elapsed time.
    #    \return The current total running elapsed time (h).
    def elapsedTime(self):
        return self.m_elapsed


    #/*! \brief Reports whether the client has called initializeEnvironment().
    #    \retval TRUE if client has called initializeEnvironment().
    #    \retval FALSE if client has not called initializeEnvironment().
    def initialized(self):
        return self.m_init



    #/*! \brief Initializes a dead fuel moisture stick's internal and external
    #environment.

    #Initializes the stick's internal and external environmental variables.
    #The stick's internal temperature and water content are presumed to be
    #uniformly distributed.

    #\param[in] ta Initial ambient air temperature (oC).
    #\param[in] ha Initial ambient air relative humidity (g/g).
    #\param[in] sr Initial solar radiation (W/m2).
    #\param[in] rc Initial cumulative rainfall amount (cm).
    #\param[in] ti Initial stick temperature (oC).
    #\param[in] hi Initial stick surface relative humidty (g/g).
    #\param[in] wi Initial stick fuel moisture fraction (g/g).
    #\param[in] bp Initial stick barometric pressure (cal/cm3).


    def initializeEnvironment (self, ta,ha,sr,rc,ti,hi,wi,bp):

        self.m_ta0 = self.m_ta1 = ta  # Previous and current ambient air temperature (oC)
        self.m_ha0 = self.m_ha1 = ha  # Previous and current ambient air relative humidity (g/g)
        self.m_sv0 = self.m_sv1 = sr / Smv  # Previous and current solar insolation (millivolts)
        self.m_rc0 = self.m_rc1 = rc  # Previous and current cumulative rainfall amount (cm)
        self.m_ra0 = self.m_ra1 = 0.  # Previous and current observation period's rainfall amount (cm)
        self.m_bp0 = self.m_bp1 = bp  # Previous and current observation barometric pressure (cal/m3)

        # Stick initialization
        self.m_hf = hi  # Relative humidity at fuel surface (g/g)
        self.m_wfilm = 0.0  # Water film moisture contribution (g/g)
        self.m_wsa = wi + .1  # Stick fiber saturation point (g/g)
        for i in range(0, self.m_nodes):
            self.m_t[i] = ti
            self.m_w[i] = wi
            self.m_s[i] = 0.0

        self.diffusivity(self.m_bp0)
        self.m_init = True
        return




    #/*! \brief Initializes a dead fuel moisture stick's internal and external
    #environment and its DateTime.

    #Initializes the stick's internal and external environmental variables.
    #The stick's internal temperature and water content are presumed to be
    #uniformly distributed.

    #This overloaded version also initializes the stick's Julian date from
    #the passed date and time parameters.  The first call to update()
    #calculates elapsed time from the date/time passed here.

    #\param[in] year     Observation year (4 digits).
    #\param[in] month    Observation month (Jan==1, Dec==12).
    #\param[in] day      Observation day-of-the-month [1..31].
    #\param[in] hour     Observation elapsed hours in the day [0..23].
    #\param[in] minute   Observation elapsed minutes in the hour (0..59].
    #\param[in] second   Observation elapsed seconds in the minute [0..59].
    #\param[in] ta Initial ambient air temperature (oC).
    #\param[in] ha Initial ambient air relative humidity (g/g).
    #\param[in] sr Initial solar radiation (W/m2).
    #\param[in] rc Initial cumulative rainfall amount (cm).
    #\param[in] ti Initial stick temperature (oC).
    #\param[in] hi Initial stick surface relative humidty (g/g).
    #\param[in] wi Initial stick fuel moisture fraction (g/g).
    #\param[in] bp Initial stick barometric pressure (cal/cm3).
    def initializeEnvironmentAndTime(self, year, month, day, hour, minute, second, ta, ha, sr, rc, ti, hi, wi, bp):

        self.m_datetime = datetime(year, month, day, hour, minute, second, tzinfo=TZ()).isoformat(' ')
        self.initializeEnvironment(ta,ha,sr,rc,ti,hi,wi,bp)
        return

    #/*! \brief Initializes a DeadFuelMoisture stick with model parameters inferred
    #    from the \a radius.

    #    \param[in] radius Dead fuel stick radius (cm).
    #    \param[in] name Name or description of the dead fuel stick.
    #Initializes a dead fuel moisture stick, applying the passed parameters.

    #\param[in] name Name or description of the dead fuel stick.
    #\param[in] radius Dead fuel stick radius (cm).
    #\param[in] stickNodes Number of stick nodes in the radial direction
    #            [required]
    #\param[in] moistureSteps Number of moisture content computation steps
    #            per observation [required].
    #\param[in] diffusivitySteps Number of diffusivity computation steps per
    #            observation [required].
    #\param[in] planarHeatTransferRate Stick planar heat transfer rate
    #            (cal/cm2-h-C) [required].
    #\param[in] adsorptionRate Adsorption surface mass transfer rate
    #            ((cm3/cm2)/h) [required].
    #\param[in] rainfallRunoffFactor Rain runoff factor (dl) [required].
    #\param[in] stickLength Stick length (cm) [optional, default = 41 cm].
    #\param[in] stickDensity Stick density (g/cm3) [optional, default = 0.40].
    #\param[in] waterFilmContribution Water film contribution to stick moisture
    #            content (g water/g dry fuel) [optional, default = 0].
    #\param[in] localMaxMc Stick maximum local moisture due to rain
    #            (g water/g dry fuel) [optional, default = 0.6].
    #\param[in] desorptionRate Desorption surface mass transfer rate
    #            ((cm3/cm2)/h) [optional, default = 0.06].
    #\param[in] rainfallAdjustmentFactor [not used].
    #\param[in] stormTransitionValue Storm transition value (cm/h) [not used].
    #\param[in] randseed If not zero, nodal temperature, saturation, and
    #            moisture contents are pertibated by some small amount.
    #            If > 0, this value is used as the seed.
    #            If < 0, uses system clock for seed.
    #\param[in] allowRainfall2 If TRUE, applies Nelson's logic for rainfall runoff
    #            after the first hour [default=true].
    #\param[in] allowRainstorm If TRUE, applies Nelson's logic for rainstorm transition
    #            and state [default=true].
    #\param[in] pertubateColumn If TRUE, whenever the continuous liquid column condition
    #            occurs, the nodal moistures get pertubated by a small amount [default=true].
    #\param[in] rampRai0 If TRUE, uses Bevins' ramping of rainfall runoff factor
    #            rather than Nelsons rai0 *= 0.15 [default=false].
    def initializeParameters(self, radius, name):

        # Start with everything set to zero
        self.zero()

        # Constrain and store the passed parameters
        self.m_name = name
        self.m_randseed = 0  # random number seed (none)
        self.m_radius = radius
        self.m_length = 41.0  # stick length (cm)
        self.m_density = 0.400  # stick density (g/cm3)
        self.m_dSteps = self.deriveDiffusivitySteps(radius)
        self.m_hc = self.derivePlanarHeatTransferRate(radius)
        self.m_nodes = self.deriveStickNodes(radius)
        self.m_rai0 = self.deriveRainfallRunoffFactor(radius)
        self.m_rai1 = 0.5  # rainfall adjustment factor "rai1" (not used)
        self.m_stca = self.deriveAdsorptionRate(radius)
        self.m_stcd = 0.06  # desorption rate "stcd"
        self.m_mSteps = self.deriveMoistureSteps(radius)


        self.m_stv = 9999.  # storm transition value "stv" (not used)
        self.m_wmx = 0.35  # maximum local moisture content "wmx"
        self.m_wfilmk = 0.0  # water film contribution "wfilmk"
        self.m_allowRainfall2 = False  # If TRUE, applies Nelson's logic for rainfall runoff after the first hour
        self.m_allowRainstorm = False  # If TRUE, applies Nelson's logic for rainstorm transition and state
        self.m_pertubateColumn = True  # If TRUE, the continuous liquid column condition get pertubated
        self.m_rampRai0 = True  # If TRUE, used Bevins' ramping of rainfall runoff factor rather than Nelsons rai0 *= 0.15
        # Initialize all other stick parameters and intermediates
        self.initializeStick()
        return


    #/*! \brief Initializes a dead fuel moisture stick's model parameters.
    #
    #    \note Should be called after and set<parameter>().
    def initializeStick(self):
        # Should we randomize nodal moisture, saturation, and temperatures by some
        # small, insignificant amount to introduce computational stability when
        # propagating these parameters within the stick?
        # If > 0, the value is used as a random generator seed.
        # If < 0, the system clock is used to get the seed.
        self.setRandomSeed(self.m_randseed)

        # Internodal distance (cm)
        self.m_dx = self.m_radius / float(self.m_nodes - 1)
        self.m_dx_2 = self.m_dx * 2.

        # Maximum possible stick moisture content (g/g)
        self.m_wmax = ( 1. / self.m_density ) - ( 1. / 1.53 )

        for i in range(0, self.m_nodes):
            # Initialize ambient air temperature to 20 oC
            self.m_t.append(20.0)

            # Initialize fiber saturation point to 0 g/g
            self.m_s.append(0.0)

            # Initialize bound water diffusivity to 0 cm2/h
            self.m_d.append(0.0)

            # Initialize moisture content to half the local maximum (g/g)
            #self.m_w.append(0.5 * self.m_wmx)
            self.m_w.append(0.2)

        # Derive nodal radial distances

        for i in range(0, self.m_nodes - 1):
            # Initialize radial distance from center of stick (cm)
            self.m_x.append(self.m_radius - self.m_dx * i)
        self.m_x.append(0.0)


        # Derive nodal volume fractions

        ro = self.m_radius
        ri = ro - 0.5 * self.m_dx
        a2 = self.m_radius ** 2
        self.m_v.append((ro ** 2 - ri ** 2 ) / a2)
        vwt = self.m_v[0]
        for i in range(1, self.m_nodes - 1):
            ro = ri
            ri = ro - self.m_dx
            self.m_v.append(( ro ** 2 - ri ** 2) / a2)
            vwt = vwt + self.m_v[i]

        self.m_v.append(ri ** 2 / a2)
        vwt = vwt + self.m_v[self.m_nodes - 1]

        # Added by Stuart Brittain on 1/14/2007
        # for performance improvement in update()
        for i in range(0, self.m_nodes):
            self.m_Twold.append(0.0)
            self.m_Ttold.append( 0.0)
            self.m_Tsold.append(0.0)
            self.m_Tv.append(0.0)
            self.m_To.append(0.0)
            self.m_Tg.append(0.0)

        # Initialize the environment, but set m_init to FALSE when done
        self.initializeEnvironment(
            20.,  # Ambient air temperature (oC)
            0.20,  # Ambient air relative humidity (g/g)
            0.0,  # Solar radiation (W/m2)
            0.0,  # Cumulative rainfall (cm)
            20.0,  # Initial stick temperature (oC)
            0.20,  # Initial stick surface humidity (g/g)
            0.2,  # Initial stick moisture content
            0.0218  # Default value for barometric pressure
        )

        self.m_init = False

        #-------------------------------------------------------------------------
        # Computation optimization parameters
        #-------------------------------------------------------------------------

        # m_hwf == hw and aml computation factor used in update()
        self.m_hwf = 0.622 * self.m_hc * pow(( Pr / Sc), 0.667)

        # m_amlf == aml optimization factor
        self.m_amlf = self.m_hwf / ( 0.24 * self.m_density * self.m_radius )

        # m_capf = cap optimization factor. */
        rcav = 0.5 * Aw * Wl
        self.m_capf = 3600. * Pi * St * rcav * rcav / ( 16. * self.m_radius * self.m_radius * self.m_length * self.m_density )

        # m_vf == optimization factor used in update()
        # WAS: m_vf = St / (Aw * Wl * Scr)
        # WAS: m_vf = St / (Wl * Scr)
        self.m_vf = St / ( self.m_density * Wl * Scr )
        return




    #/*! \brief Access to the stick's current maximum local moisture content.
    #       \return Stick's current maximum local moisture content (g/g).
    def maximumLocalMoisture(self):
        return self.m_wmx


    #/*! \brief Determines the mean moisture content of the stick's radial profile.

    #    \note
    #    The integral average of the stick's nodal moisture contents is calculated
    #    without consideration to the nodes' volumetric representation.

    #    \deprecated Use Fms_MeanWtdMoisture() for a volume-weighted mean
    #    moisture content.

    #    \return The mean moisture content of the stick's radial profile (g/g).
    def meanMoisture(self):
        wec = self.m_w[0]
        wei = self.m_dx / ( 3. * self.m_radius )
        for i in range(1, self.m_nodes - 1):
            wea = 4. * self.m_w[i]
            web = 2. * self.m_w[i + 1]
            if ( ( i + 1 ) == ( self.m_nodes - 1 ) ):
                web = self.m_w[self.m_nodes - 1]

            wec = wec + (web + wea)

        wbr = wei * wec
        if ( wbr > self.m_wmx ):
            wbr = self.m_wmx
        # Add water film
        wbr = wbr + m_wfilm
        return wbr



    #/*! \brief Determines the volume-weighted mean moisture content of the
    #    stick's radial profile.

    #    \return The volume-weighted mean moisture content of the stick's radial
    #    profile (g/g).
    def meanWtdMoisture(self):
        wbr = 0.0

        for i in range(0, self.m_nodes):
            wbr = wbr + self.m_w[i] * self.m_v[i]


        if (wbr >self.m_wmx):
            wbr = self.m_wmx
        else:
            wbr = wbr

        # Add water film
        wbr = wbr + self.m_wfilm
        return wbr
    def medianRadialMoisture(self):
        return np.median(self.m_w)
	

    def medianMoisture(self):
        My_m_w = sorted(self.m_w)
        n = len(My_m_w)
        n2 = int(n / 2)
        return My_m_w[n2]
    
    def setMoisture(self, initFM):
        for i in range(0,self.m_nodes):
            self.m_w[i] = initFM
	
	
#/*! \brief Determines the volume-weighted mean temperature of the stick's
    #    radial profile.
    #    \return The volume-weighted mean temperature of the stick's radial
    #    profile (oC).
    def meanWtdTemperature(self):
        wbr = 0.0

        for i in range(0, self.m_nodes):
            wbr = wbr + self.m_t[i] * self.m_v[i]

        return wbr


    ##/*! \brief Access to the stick's number of moisture content computation time
    #    steps per observation.
    #   \return Current number of moisture content computation time steps
    #     per observation.
    def moistureSteps(self):
        return self.m_mSteps


    #/*! \brief Access to the stick's name.#
    #
    #    \return The stick's name.
    # */

    def name(self):
        return ( self.m_name )


    #/*! \brief Access to the stick's current precipitation rate (cm/h).
    #    \return The current (most recent) precipitation rate (cm/h).
    def pptRate(self):
        if self.m_et > 0.00:
            return m_ra1 / m_et
        else:
            return 0.00


    #/*! \brief Access to the stick's current planar heat transfer rate..

    #    \return Stick's current planar heat transfer rate
    #        (\f$\frac {cal} {cm^{2} \cdot h \cdot C}\f$).
    def planarHeatTransferRate(self):
        return self.m_hc


    #/*! \brief Access to the stick's current rainfall runoff factor.
    #
    #    \return Stick's current rainfall runoff factor (dl).

    def rainfallRunoffFactor(self):
        return ( self.m_rai0 )


    #/*! \brief Updates the stick's adsorption rate.
    #    \param[in] adsorptionRate Adsorption surface mass transfer rate ((cm3/cm2)/h).
    def setAdsorptionRate(self, adsorptionRate):
        self.m_stca = adsorptionRate
        return
    
    #/*! \brief Updates the stick's desorption rate.
    #    \param[in] desorptionRate Desorption surface mass transfer rate ((cm3/cm2)/h).
    def setDesorptionRate(self, desorptionRate):
        self.m_stcD = desorptionRate
        return

    #/*! \brief Updates the stick's configuration to toggle Nelson's logic
    #    for rainfall runoff factor after the first hour of rain.
    #    \param[in]  allow (default=true).
    def setAllowRainfall2(self, allow):
        self.m_allowRainfall2 = allow
        return


    #/*! \brief Updates the stick's configuration to toggle Nelson's logic
    #    to allow the Rainstorm state and rainfall-to-rainstorm transition.
    #   \param[in]  allow (default=true).
    def setAllowRainstorm(self, allow):
        self.m_allowRainstorm = allow
        return


    #/*! \brief Updates the stick's desorption rate.
    #    \param[in] desorptionRate Desorption surface mass transfer rate ((cm3/cm2)/h).
    def setDesorptionRate(self, desorptionRate):
        self.m_stcd = desorptionRate
        return


    #/*! \brief Updates the stick's diffusivity computation steps per update().
    #    \param[in] diffusivitySteps Number of diffusivity computation steps per observation.
    def setDiffusivitySteps(self, diffusivitySteps):
        self.m_dSteps = diffusivitySteps
        return


    #/*! \brief Updates the stick's maximum local moisture content.
    #    \param[in] localMaxMc Stick maximum local moisture due to rain (g water/g dry fuel).
    def setMaximumLocalMoisture(self, localMaxMc):
        self.m_wmx = localMaxMc
        return


    #/*! \brief Updates the stick's moisture content computation steps per update().
    #    \param[in] moistureSteps Number of moisture content computation steps per observation.
    def setMoistureSteps(self, moistureSteps):
        self.m_mSteps = moistureSteps
        return


    #/*! \brief Updates the stick's planar heat transfer rate.
    #    \param[in] planarHeatTransferRate Stick planar heat transfer rate (cal/cm2-h-C).
    def setPlanarHeatTransferRate(self, planarHeatTransferRate):
        self.m_hc = planarHeatTransferRate
        return


    #/*! \brief Updates the stick's column pertubation configuration.
    #
    #    \param[in] pertubate (default=true).
    def setPertubateColumn(self, pertubate):
        self.m_pertubateColumn = pertubate
        return


    #/*! \brief Updates the stick's rainfall runoff factor.
    #
    #    \param[in] rainfallRunoffFactor Rain runoff factor during the initial hour of rainfall (dl).
    def setRainfallRunoffFactor(self, rainfallRunoffFactor):
        self.m_rai0 = rainfallRunoffFactor
        return


    #/*! \brief Updates the stick's configuration to toggle Bevins' logic
    #    that ramps the rainfall runoff factor during falling humidity, rather than
    #    use Nelson's rai0 *= 0.15.
    #    \param[in]  ramp (default=true).
    def setRampRai0(self, ramp):
        self.m_rampRai0 = ramp
        return


    #/*! \brief Sets the ransom seed behavior.
    #    \param[in] randseed If not zero, nodal temperature, saturation,
    #    and moisture contents are pertubated by some small amount.
    #    If > 0, this value is used as the seed. If < 0, uses system clock for seed.
    def setRandomSeed(self, randseed):
        self.m_randseed = randseed
        if ( self.m_randseed > 0 ):

            random.seed(self.m_randseed)
        else:
            random.seed(time())

        return


    #/*! \brief Updates the stick density.
    #    \param[in] stickDensity Stick density (g/cm3) [optional, default = 0.40].
    def setStickDensity(self, stickDensity):
        self.m_density = stickDensity
        return


    #/*! \brief Updates the stick length.
    #    \param[in] stickLength Stick length (cm) [optional, default= 41 cm].
    def setStickLength(self, stickLength):
        self.m_length = stickLength
        return


    #/*! \brief Updates the number of stick radial computation nodes.
    #    \param[in] stickNodes Number of stick nodes in the radial direction [optional, default = 11].
    def setStickNodes(self, stickNodes):
        self.m_nodes = stickNodes
        return


    #/*! \brief Updates the water film contribution to stick weight.
    #    \param[in] waterFilm Water film contribution to stick moisture content (g water/g dry fuel).
    def setWaterFilmContribution(self, waterFilm):
        self.m_wfilm = waterFilm
        return


    #/*! \brief Access to prevailing state for most recent update.
    #    \return Prevailing state for most recent update.
    def state(self):
        return self.m_state

    #/*! \brief Access to the stick's current state name.
    #    \return The current state name.
    def stateName(self):
        States = {0: "None", 1: "Adsorption", 2: "Desorption", 3: "Condensation1", 4: "Condensation2", 5: "Evaporation", 6: "Rainfall1", 7:"Rainfall2", 8:"Rainstorm", 9: "Stagnation", 10: "Error"}
        return States[self.m_state]

    #/*! \brief Access to the stick's particle density.
    #    \return Current stick density (g/cm3).
    def stickDensity(self):
        return self.m_density  #/*! \brief Access to the stick's length.
    #    \return Current stick length (cm).
    def stickLength(self):
        return self.m_length

    #/*! \brief Access to the stick's number of moisture content radial computation
    #    nodes.
    #    \return Current number of moisture content radial computation  nodes.
    #     per observation.
    def stickNodes(self):
        return self.m_nodes

    #/*! \brief Access to the stick's surface fuel moisture content.
    #    \return The current surface fuel moisture content (g/g).
    def surfaceMoisture(self):
        return self.m_w[0]


    #/*! \brief Access to the stick's surface fuel temperature.
    #    \return The current surface fuel temperature (oC).
    def surfaceTemperature(self):
        return ( self.m_t[0] )

    #/*! \brief Static convenience method to derive a random number uniformly
    #    distributed in the range [\a min .. \a max].
    #    \param[in] min  Minimum range value.
    #    \param[in] max  Maximum range value.
    #    Uses the system rand() to generate the number.
    #    \return A uniformly distributed random number within [\a min .. \a max].
    def uniformRandom( min, max ):
        return random.uniform(min, max)

    #/*! \brief Updates a dead moisture stick's internal and external environment
    #    based on the current weather observation values.
    #    This overloaded version accepts the current date and time as arguments,
    #    and automatically calculates elapsed time since the previous update().
    #    \note The client must have called the corresponding initializeEnvironment()
    #    method that accepts date and time arguments to ensure that the date and
    #    time has been initialized.
    #    \note The client program should try it's hardest to catch all bad input data
    #    and corrected or duplicated records \b before calling this method.
    #    \param[in] year     Observation year (4 digits).
    #    \param[in] month    Observation month (Jan==1, Dec==12).
    #    \param[in] day      Observation day-of-the-month [1..31].
    #    \param[in] hour     Observation elapsed hours in the day [0..23].
    #    \param[in] minute   Observation elapsed minutes in the hour (0..59].
    #   \param[in] second   Observation elapsed seconds in the minute [0..59].
    #    \param[in] at   Current observation's ambient air temperature (oC).
    #    \param[in] rh   Current observation's ambient air relative humidity (g/g).
    #    \param[in] sW   Current observation's solar radiation (W/m2).
    #    \param[in] rcum Current observation's total cumulative rainfall amount (cm).
    #    \param[in] bpr  Current observation's stick barometric pressure (cal/cm3).
    #    \retval TRUE if all inputs are ok and the stick is updated.
    #    \retval FALSE if inputs are out of range and the stick is \b not updated.

    def date_to_julian_day(my_date):
        """Returns the Julian day number of a date."""
        a = (14 - my_date.month) // 12
        y = my_date.year + 4800 - a
        m = my_date.month + 12 * a - 3
        return my_date.day + ((153 * m + 2) // 5) + 365 * y + y // 4 - y // 100 + y // 400 - 32045


    def updateAll(self, year, month, day, hour, minute, second, at, rh, sW, rcum, bpr):
        # Determine Julian date for this new observation
        #jd0 = date_to_julian_day(self.m_datetime)
        #self.m_datetime = datetime(year, month, day, hour, minute, second, tzinfo=TZ()).isoformat(' ')
        #jd1 = date_to_julian(self.m_datetime)

        # Determine elapsed time (h) between the current and previous dates
        #et = 24.0 * ( jd1 - jd0 )
        et = 1
        return ( self.update(et, at, rh, sW, rcum, bpr) )



    #/*! \brief Updates a dead moisture stick's internal and external environment
    #    based on the passed (current) weather observation values.

    # This overloaded version accepts the elapsed time since the previous
    # observation, and does not automatically update the Julian date.

    #\param[in] et   Elapsed time since the previous observation (h).
    #    \param[in] at   Current observation's ambient air temperature (oC).
    #    \param[in] rh   Current observation's ambient air relative humidity (g/g).
    #    \param[in] sW   Current observation's solar radiation (W/m2).
    #    \param[in] rcum Current observation's total cumulative rainfall amount (cm).
    #    \param[in] bpr  Current observation's stick barometric pressure (cal/cm3).

    #    \note The client program should try it's hardest to catch all bad input data
    #    and corrected or duplicated records \b before calling this method.

    #    \retval TRUE if all inputs are ok and the stick is updated.
    #    \retval FALSE if inputs are out of range and the stick is \b not updated.

    def update(self,et, at, rh, sW, ppt, bpr = 0.0218,sv_a = 0.0000239, sv_b = 20.58, sv_c = 5205):

        # Increment update counter
        self.m_updates = self.m_updates + 1
        self.m_elapsed = self.m_elapsed + et
        mpre = self.m_w[0]
        # If the elapsed time < 0.00027 hours (<0.01 sec), then treat this as
        # a duplicate or corrected observation and return
        if ( et < 0.0000027 ):
            print ("update has a regressive elapsed time of "+ str(et))
            return False

        # Cumulative rainfall must equal or exceed its previous value
        #if ( rcum < self.m_rc1 ):
        #    print ("DeadFuelMoisture::update() rcum < rc1")
        #    print (rcum, self.m_rc1)

            return False
        # Relative humidity must be reasonable
        if ( rh < 0.001 or rh > 1.0 ):
            print ("DeadFuelMoisutre::update RH out of range: %s" % (rh))
            rh = 0.5
            #return False

        # Ambient temperature must be reasonable
        if ( at < -60. or at > 60. ):
            print ("DeadFuelMoisture::update() Temp out of range: %s"% (at)) 
            return False
        # Insolation must be reasonable

        if sW < 0.0:
            sW = 0.0
        if ( sW > 2000. ):
            print ("DeadFuelMoisture::update() Radiation out of range: %s"% (sW)) 
        self.m_rcum = self.m_rcum + ppt
        # First save the previous weather observation values
        self.m_ta0 = self.m_ta1  # Previous air temperature (oC)
        self.m_ha0 = self.m_ha1  # Previous air relative humidity (g/g)
        self.m_sv0 = self.m_sv1  # Previous pyranometer voltage (millivolts)
        self.m_rc0 = self.m_rc1  # Previous cumulative rainfall (cm)
        self.m_ra0 = self.m_ra1  # Previous period's rainfall amount (cm)
        self.m_bp0 = self.m_bp1  # Previous barometric pressure (cal/m3)

        # Then save the current weather observation values
        self.m_ta1 = at  # Current air temperature (oC)
        self.m_ha1 = rh  # Current air relative humidity (g/g)
        self.m_sv1 = sW / Smv  # Current pyranometer voltage (millivolts)
        self.m_rc1 = self.m_rcum  # Current cumulative rainfall (cm)
        self.m_bp1 = bpr  # Current barometric pressure (cal/m3)
        self.m_et = et  # Current elapsed time since previous update (h)

        # Precipitation amount since last observation
        self.m_ra1 = self.m_rc1 - self.m_rc0
        self.m_ra1 = ppt
        # If no precipitation, reset the precipitation duration timer
        if (self.m_ra1 < 0.0001):
            self.m_rdur = 0.0


        # Precipitation rate since last observation adjusted by Pi (cm/h)
        self.m_pptrate = self.m_ra1 / et / Pi
        # Determine moisture computation time step interval (h)
        if(et > 1):
            self.m_mdt = (1 / (float(self.m_mSteps) * et/2))
        else:
            self.m_mdt = et / float(self.m_mSteps)
        self.m_mdt = et / float(self.m_mSteps)
        self.m_mdt_2 = self.m_mdt * 2.
        # Nelson's "s" factor used in update() loop
        self.m_sf = 3600. * self.m_mdt / ( self.m_dx_2 * self.m_density )
        # Determine bound water diffusivity time step interval (h)
        self.m_ddt = et / float(self.m_dSteps)

        # First hour runoff factor h-(g/(g-h))
        rai0 = self.m_mdt * self.m_rai0 * ( 1.0 - exp(-100. * self.m_pptrate) )

        # Adjustment for rainfall cases when humidity is dropping
        if ( self.m_ha1 < self.m_ha0 ):

            # If Bevin's rainfall ramping is used
            if ( self.m_rampRai0 ):
                rai0 *= ( 1.0 - ( ( self.m_ha0 - self.m_ha1 ) / self.m_ha0 ) )

            else:
                rai0 *= 0.15
        # Subsequent runoff factor h-(g/(g/h))
        rai1 = self.m_mdt * self.m_rai1 * self.m_pptrate

        # DFM state counter
        tstate = [11]
        for i in range(0, 11):
            tstate.append(0)

        # Next time (tt) to run diffusivity computations.
        ddtNext = self.m_ddt
        # Elapsed moisture computation time (h)
        tt = self.m_mdt

        continuousLiquid = False

        # Loop for each moisture time step between environmental inputs.
        nstep = 0
        # Somehwere here we need an et > 1 case

        while tt <= et:
            #for ( int nstep=1 tt <= et tt = nstep*m_mdt, nstep++ )

            nstep = nstep + 1
            # Fraction of time elapsed between previous and current obs (dl)
            tfract = tt / et
            # Air temperature interpolated between previous and current obs (oC)
            ta = self.m_ta0 + ( self.m_ta1 - self.m_ta0 ) * tfract
            # Air humidity interpolated between previous and current obs (dl)
            ha = self.m_ha0 + ( self.m_ha1 - self.m_ha0 ) * tfract
            # Solar radiation interpolated between previous and current obs (millivolts)
            sv = self.m_sv0 + ( self.m_sv1 - self.m_sv0 ) * tfract
            # Barometric pressure interpolated between previous and current obs (bal/m3)
            bp = self.m_bp0 + ( self.m_bp1 - self.m_bp0 ) * tfract
            # Fraction of the solar constant interpolated between obs (mv)
            # double fsc = 0.07 * sv
            fsc = sv / Srf
            # Ambient air temperature (oK)
            tka = ta + Kelvin
            # Dew point temperature (oK)
            tdw = 5205. / ( ( 5205. / tka ) - log(ha) )
            # Dew point temperature (oC)
            tdp = tdw - Kelvin
            # Sky temperature (oK)
            # Long wave radiative surface heat transfer coefficient (cal/cm2-h-C)
            if (fsc < 0.000001 ):
                tsk = Tcn + Kelvin
                hr = Hrn
                sr = 0.0
            else:
                tsk = Tcd + Kelvin
                hr = Hrd
                sr = Srf * fsc
            # Water saturation vapor pressure in ambient air (cal/cm3)
            psa = sv_a * exp(sv_b - ( sv_c / tka ))
            #psa = 0.0000239 * exp(20.4 - ( 5205. / tka ))
            # From Campbell and Norman
            #psa = 0.611*exp((17.502*ta)/(240.97+ta)) / 4184
            # Water saturation vapor pressure in air (cal/cm3)
            pa = ha * psa
            # Water saturation vapor pressure at dewpoint (cal/cm3)
            psd = sv_a * exp(sv_b - ( sv_c / tdw ))
            #psd = 0.611*exp((17.502*(tdw-Kelvin))/(240.97+(tdw-Kelvin))) / 4184
            #psd = 0.0000239 * exp(20.4 - ( 5205. / tdw ))
            # Rainfall duration (h)
            if (self.m_ra1 > 0.0001):
                self.m_rdur = self.m_rdur + self.m_mdt
            else:
                self.m_rdur = 0.0

            #----------------------------------------------------------------------
            # Stick surface temperature and humidity
            #----------------------------------------------------------------------
            # Intermediate stick surface temperature (oC)
            tfd = ta + ( sr - hr * ( ta - tsk + Kelvin ) ) / ( hr + self.m_hc )
            # Latent heat of vaporization of water (cal/mole)
            qv = 13550. - 10.22 * ( tfd + Kelvin )
            # Differential heat of sorption of water (cal/mole)
            qw = 5040. * exp(-14. * self.m_w[0])
            # Stick heat transfer coefficient for vapor diffusion above FSP
            hw = ( self.m_hwf * Ap / 0.24 ) * qv / 18.
            # Stick surface temperature (oC)
            self.m_t[0] = tfd - ( hw * ( tfd - ta ) / ( hr + self.m_hc + hw ) )
            # Stick surface temperature (oK)
            tkf = self.m_t[0] + Kelvin
            #print "TKF:" ,self.m_t[0],tkf
            # Kinematic viscosity of liquid water (cm2/s)
            gnu = 0.00439 + 0.00000177 * ( 338.76 - tkf )**2.1237

            # EMC sorption isotherm parameter (g/g)
            c1 = 0.1617 - 0.001419 * self.m_t[0]
            # EMC sorption isotherm parameter (g/g)
            c2 = 0.4657 + 0.003578 * self.m_t[0]
            # Stick fiber saturation point (g/g)
            self.m_wsa = c1 * pow(Wsf, c2)

            # Maximum minus current fiber saturation (g/g)
            wdiff = self.m_wmax - self.m_wsa
            if (wdiff < 0.000001):
                wdiff = 0.000001

            # Water saturation vapor pressure at surface temp (cal/cm3)
            ps1 = 0.0000239 * exp(20.58 - ( 5205. / tkf ))
            # Water vapor pressure at the stick surface (cal/cm3)
            p1 = pa + Ap * bp * ( qv / (qv + qw) ) * ( tka - tkf )
            if (p1 < 0.000001):
                p1 = 0.000001

            # Stick surface humidity (g/g)
            self.m_hf = p1 / ps1
            if (self.m_hf > Hfs):
                self.m_hf = Hfs

            # Stick equilibrium moisture content (g/g). */
            hf_log = -math.log(1. - self.m_hf)
            self.m_sem = c1 * pow(hf_log, c2)

            #----------------------------------------------------------------------
            # Stick surface moisture content
            #----------------------------------------------------------------------

            # Initialize state for this m_mdt
            self.m_state = self.DFM_State_None
            # Start with no water film contribution
            self.m_wfilm = 0.
            # Factor related to rate of evaporation or condensation ((g/g)/h)
            self.aml = 0.0
            # Mass transfer biot number (dl)
            self.bi = 0.0
            # Previous and new value of m_w[0] (g/g) and m_s[0]
            s_new = self.m_s[0]
            w_new = self.m_w[0]
            w_old = self.m_w[0]
            #print ps1,psd,(ps1-psd)
            #......1: If it is RAINING:
            if ( self.m_ra1 > 0.0 ):

                #..........1a: If this is a RAINSTORM:
                if ( self.m_allowRainstorm and self.m_pptrate >= self.m_stv ):
                    self.m_state = self.DFM_State_Rainstorm
                    self.m_wfilm = self.m_wfilmk
                    w_new = self.m_wmx

                #..........1b: Else this is RAINFALL:
                else:
                   # print "Rainfall"
                   # First hour of rainfall
                    if ( self.m_rdur < 1.0 or not self.m_allowRainfall2 ):

                        self.m_state = self.DFM_State_Rainfall1
                        w_new = w_old + rai0
                        #print w_old,rai0,w_new


                    else:

                        # Note: if m_allowRainfall2 is False, this else will never be reached
                        self.m_state = self.DFM_State_Rainfall2
                        w_new = w_old + rai1

                self.m_wfilm = self.m_wfilmk
                s_new = ( w_new - self.m_wsa ) / wdiff
                self.m_t[0] = tfd   # Stick surface temperature

                self.m_hf = Hfs     # Stick surface humidity



            #......2: Else it is not raining:
            else:

                #.........2a: If moisture content exceeds the fiber saturation point:
                if ( w_old > self.m_wsa ):
                    p1 = ps1
                    self.m_hf = Hfs

                    # Factor related to evaporation or condensation rate ((g/g)/h)
                    aml = self.m_amlf * (ps1 - psd) / bp

                    self.m_aml = aml
                    if ( self.m_t[0] <= tdp and p1 > psd ):
                        aml = 0.

                    w_new = w_old - aml * self.m_mdt_2
                    #print w_old,w_new,aml,self.m_mdt_2,tt,et
                    # Evaporation
                    if ( aml > 0. ):
                        w_new = w_new - ( self.m_mdt * self.m_capf / gnu )
                        #print (tkf,gnu,self.m_capf,self.m_mdt,self.m_mdt * self.m_capf/gnu)
                    if (w_new > self.m_wmx):
                        w_new = self.m_wmx

                    s_new = ( w_new - self.m_wsa ) / wdiff

                    #..............2a1: if moisture content is rising: CONDENSATION
                    if ( w_new > w_old ):
                        self.m_state = self.DFM_State_Condensation1

                    #..............2a2: else if moisture content is steady: STAGNATION
                    elif ( w_new == w_old ):
                        self.m_state = self.DFM_State_Stagnation

                        #..............2a3: else if moisture content is falling: EVAPORATION
                    elif ( w_new < w_old ):
                        self.m_state = self.DFM_State_Evaporation


                #..........2b: else if fuel temperature is less than dewpoint: CONDENSATION
                elif ( self.m_t[0] <= tdp ):

                    self.m_state = self.DFM_State_Condensation2
                    # Factor related to evaporation or condensation rate ((g/g)/h)
                    # If the stick surface vapor pressure is greater than the
                    # vapor pressure at dewpoint temperature
                    if (p1 > psd):
                        aml = 0.0
                    else:
                        aml = self.m_amlf * (p1 - psd) / bp

                    w_new = w_old - aml * self.m_mdt_2
                    if(aml > 0):
                        print(p1,psd,aml,self.m_mdt_2)
                    s_new = ( w_new - self.m_wsa ) / wdiff

                #..........2c: else surface moisture content less than fiber saturation point
                #              and stick temperature greater than dewpoint ...
                else:

                    #..............2c1: if surface moisture greater than equilibrium: DESORPTION
                    if ( w_old >= self.m_sem ):
                        self.m_state = self.DFM_State_Desorption
                        self.bi = self.m_stcd * self.m_dx / self.m_d[0]

                    #..............2c2: else surface moisture less than equilibrium: ADSORPTION
                    else:
                        self.m_state = self.DFM_State_Adsorption
                        self.bi = self.m_stca * self.m_dx / self.m_d[0]

                    w_new = ( self.m_w[1] + self.bi * self.m_sem ) / ( 1. + self.bi )
                    #print self.bi * self.m_sem,1+self.bi,self.m_stca,self.m_stcd,self.m_dx,self.m_d[0],self.m_sem,w_old,self.m_w[1] - w_new

                    s_new = 0.

            # end of not raining
            # Store the new surface moisture and saturation
            if (w_new > self.m_wmx):
                self.m_w[0] = self.m_wmx
            else:
                self.m_w[0] = w_new

            if (s_new < 0):
                self.m_s[0] = 0.0
            else:
                self.m_s[0] = s_new

            tstate[self.m_state] = tstate[self.m_state] + 1

            #----------------------------------------------------------------------
            # Compute interior nodal moisture content values.
            #----------------------------------------------------------------------
            for i in range(0, self.m_nodes):
                self.m_Twold[i] = self.m_w[i]
                self.m_Tsold[i] = self.m_s[i]
                self.m_Ttold[i] = self.m_t[i]
                self.m_Tv[i] = Thdiff * self.m_x[i]
                self.m_To[i] = self.m_d[i] * self.m_x[i]

            # Propagate the moisture content changes
            if ( self.m_state != self.DFM_State_Stagnation ):

                for i in range(0, self.m_nodes):
                    self.m_Tg[i] = 0.0
                    svp = ( self.m_w[i] - self.m_wsa ) / wdiff
                    if ( svp >= Sir and svp <= Scr ): #Sir-Saturation value below which liquid water columns no longer exist and Scr - Saturation value at which liquid miniscus first enters the tapered
                        # Permeability of stick when nonsaturated (cm2)
                        ak = Aks * ( 2. * sqrt(svp / Scr) - 1. )
                        # Free water transport coefficient (cm2/h)
                        self.m_Tg[i] = ( ak / ( gnu * wdiff ) ) * self.m_x[i] * self.m_vf * ( ( Scr / svp ) ** 1.5 )

                # Propagate the fiber saturation moisture content changes
                for i in range(1, self.m_nodes - 1):
                    ae = self.m_Tg[i + 1] / self.m_dx
                    aw = self.m_Tg[i - 1] / self.m_dx
                    ar = self.m_x[i] * self.m_dx / self.m_mdt
                    ap = ae + aw + ar
                    self.m_s[i] = ( ae * self.m_Tsold[i + 1] + aw * self.m_Tsold[i - 1] + ar * self.m_Tsold[i] ) / ap
                    if (self.m_s[i] > 1.):
                        self.m_s[i] = 1.

                    if (self.m_s[i] < 0):
                        self.m_s[i] = 0

                self.m_s[self.m_nodes - 1] = self.m_s[self.m_nodes - 2]

                # Check if m_s[] is less than Sir (limit of continuous liquid columns) at ANY stick node.
                continuousLiquid = True
                for i in range(1, self.m_nodes - 1):
                    if( self.m_s[i] < Sir ):
                        continuousLiquid = False
                        break

                # If all nodes have continuous liquid columns (s >= Sir) ...
                # This never happens for the 1-h or 10-h test data!
                if continuousLiquid:

                    for i in range(1, self.m_nodes):

                        self.m_w[i] = self.m_wsa + self.m_s[i] * wdiff
                        if ( self.m_pertubateColumn ):
                            rn = 0#random.uniform(-.0001, 0.0001)
                            self.m_w[i] = self.m_w[i] + rn
                        if (self.m_w[i] > self.m_wmx):
                            self.m_w[i] = self.m_wmx
                        if (self.m_w[i] < 0):
                            self.m_w[i] = 0.0
                # ... else at least one node has s < Sir.
                else:
                    # Propagate the moisture content changes
                    for i in range(1, self.m_nodes - 1):
                        ae = self.m_To[i + 1] / self.m_dx
                        aw = self.m_To[i - 1] / self.m_dx
                        ar = self.m_x[i] * self.m_dx / self.m_mdt
                        ap = ae + aw + ar
                        self.m_w[i] = ( ae * self.m_Twold[i + 1] + aw * self.m_Twold[i - 1] + ar * self.m_Twold[i] ) / ap
                        randunif = 0#random.uniform(-0.0001,0.0001)
                        self.m_w[i] = self.m_w[i] + randunif
                        if (self.m_w[i] > self.m_wmx):
                            self.m_w[i] = self.m_wmx
                        if (self.m_w[i] < 0):
                            self.m_w[i] = 0.0

            self.m_w[self.m_nodes - 1] = self.m_w[self.m_nodes - 2]

            # Propagate the fuel temperature changes
            for i in range(1, self.m_nodes - 1):
                ae = self.m_Tv[i + 1] / self.m_dx
                aw = self.m_Tv[i - 1] / self.m_dx
                ar = self.m_x[i] * self.m_dx / self.m_mdt
                ap = ae + aw + ar
                self.m_t[i] = ( ae * self.m_Ttold[i + 1] + aw * self.m_Ttold[i - 1] + ar * self.m_Ttold[i] ) / ap
                randunif = random.uniform(-0.0001,0.0001)
                self.m_t[i] = self.m_t[i] + randunif
                if (self.m_t[i] > 71.):
                    self.m_t[i] = 71.
            self.m_t[self.m_nodes - 1] = self.m_t[self.m_nodes - 2]
            # Update the moisture diffusivity if within less than half a time step
            if ( ( ddtNext - tt ) < ( 0.5 * self.m_mdt ) ):
                self.diffusivity(bp)
                ddtNext = ddtNext + self.m_ddt
            tt = nstep * self.m_mdt
            # Next moisture time step
        #print mpre-w_new,self.m_d[0]
        # Store prevailing state
        self.m_state = self.DFM_State_None
        self.max = tstate[0]

        for i in range(1, 11):
            if ( tstate[i] > self.max ):
                self.m_state = i
                self.max = tstate[i]

        return True

    #/ * ! \brief Access to the current number of observation updates.
    #\return The current number of observation updates. * /
    def updates(self):
        return ( self.m_updates )


    #/ * ! \brief Access to the stick's current water film contribution to the moisture content.
    #\return Stick's current water film contribution to the moisture content (g / g).
    def waterFilmContribution(void):
        return ( self.m_wfilmk )



    #/ * ! \brief  Sets everything to zero. * /
    def zero(self):

        self.m_density = 0.0
        self.m_dSteps = 0
        self.m_hc = 0.0
        self.m_length = 0.0
        self.m_name = ""
        self.m_nodes = 0
        self.m_radius = 0.0
        self.m_rai0 = 0.0
        self.m_rai1 = 0.0
        self.m_stca = 0.0
        self.m_stcd = 0.0
        self.m_mSteps = 0
        self.m_stv = 0.0
        self.m_wfilmk = 0.0
        self.m_wmx = 0.0
        self.m_dx = 0.0
        self.m_wmax = 0.0
        self.m_x = []
        self.m_v = []
        self.m_amlf = 0.0
        self.m_capf = 0.0
        self.m_hwf = 0.0
        self.m_dx_2 = 0.0
        self.m_vf = 0.0
        self.m_bp0 = 0.0
        self.m_ha0 = 0.0
        self.m_rc0 = 0.0
        self.m_sv0 = 0.0
        self.m_ta0 = 0.0
        self.m_init = False
        self.m_bp1 = 0.0
        self.m_et = 0.0
        self.m_ha1 = 0.0
        self.m_rc1 = 0.0
        self.m_sv1 = 0.0
        self.m_ta1 = 0.0
        self.m_ddt = 0.0
        self.m_mdt = 0.0
        self.m_mdt_2 = 0.0
        self.m_pptrate = 0.0
        self.m_ra0 = 0.0
        self.m_ra1 = 0.0
        self.m_rdur = 0.0
        self.m_sf = 0.0
        self.m_hf = 0.0
        self.m_wsa = 0.0
        self.m_sem = 0.0
        self.m_wfilm = 0.0
        self.m_elapsed = 0.0
        self.m_t = []
        self.m_s = []
        self.m_d = []
        self.m_w = []
        self.m_updates = 0
        self.m_state = self.DFM_State_None
        self.m_randseed = 0
        return


if __name__ == "__main__":
    import matplotlib.pyplot as plt
    import numpy


    T = []
    I = []
    MyStick = DeadFuelMoisture(2.0,"Little Stick")
    MyStick.setMaximumLocalMoisture(0.35)

    for i in range(1,100):
        MyStick.update(2,20,0.2,random.uniform(0,0),0,0.0218)
        T.append(MyStick.medianRadialMoisture())
        I.append(i)
        

    for i in range(1,100):
        MyStick.update(1,20,0.2,random.uniform(0,0),0.001*i,0.0218)
        T.append(MyStick.medianRadialMoisture())
        I.append(i+100)

    for i in range(1,30):
        MyStick.update(1,20,0.8,random.uniform(0,0),0.0001,0.0218)
        T.append(MyStick.medianRadialMoisture())
        I.append(i+200)

    plt.plot(I,T )

   # plt.plot(mxval2, a *mxval2, 'r', label='Fitted line')
    #plt.ylabel('Measured fine dead fuel moisture (% dry wt)')
    #plt.xlabel('Modeled fine dead fuel moisture (% dry wt)')
    #plt.ylim([0,35])
    #plt.xlim([0,35])
    plt.show()
    
    