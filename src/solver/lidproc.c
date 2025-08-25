//-----------------------------------------------------------------------------
//   lidproc.c
//
//   Project:  EPA SWMM5
//   Version:  5.2
//   Date:     07/13/23   (Build 5.2.4)
//   Author:   L. Rossman
//
//   This module computes the hydrologic performance of an LID (Low Impact
//   Development) unit at a given point in time.
//
//   Update History
//   ==============
//   Build 5.1.007:
//   - Euler integration now applied to all LID types except Vegetative
//     Swale which continues to use successive approximation.
//   - LID layer flux routines were re-written to more accurately model
//     flooded conditions.
//   Build 5.1.008:
//   - MAX_STATE_VARS replaced with MAX_LAYERS.
//   - Optional soil layer added to Porous Pavement LID.
//   - Rooftop Disconnection added to types of LIDs.
//   - Separate accounting of drain flows added.
//   - Indicator for currently wet LIDs added.
//   - Detailed reporting procedure fixed.
//   - Possibile negative head on Bioretention Cell drain avoided.
//   - Bug in computing flow through Green Roof drainage mat fixed.
//   Build 5.1.009:
//   - Fixed typo in net flux rate for vegetative swale LID.
//   Build 5.1.010:
//   - New modified version of Green-Ampt used for surface layer infiltration.
//   Build 5.1.011:
//   - Re-named STOR_INFIL to STOR_EXFIL and StorageInfil to StorageExfil to
//     better reflect their meaning.
//   - Evaporation rates from sub-surface layers reduced by fraction of
//     surface that is pervious (applies to block paver systems)
//   - Flux rate routines for LIDs with underdrains modified to produce more
//     physically meaningful results.
//   - Reporting of detailed results re-written.
//   Build 5.1.012:
//   - Modified upper limit for soil layer percolation.
//   - Modified upper limit on surface infiltration into rain gardens.
//   - Modified upper limit on drain flow for LIDs with storage layers.
//   - Used re-defined wasDry variable for LID reports to fix duplicate lines.
//   Build 5.1.013:
//   - Support added for open/closed head levels and multiplier v. head curve
//     to control underdrain flow.
//   - Support added for regenerating pavement permeability at fixed intervals.
//   Build 5.1.014:
//   - Fixed failure to initialize all LID layer moisture volumes to 0 before
//     computing LID unit performance in lidproc_getOutflow.
//   Build 5.2.0:
//   - Fixed failure to account for effect of Impervious Surface Fraction on
//     pavement permeability for Permeable Pavement LID
//   - Fixed units conversion for pavement depth in detailed report file.
//   Build 5.2.4:
//   - Modified flux limits in biocellFluxRates, pavementFluxRates and
//     trenchFluxRates.
//   - Corrected head calculation in getStorageDrainRate when unit has both
//     a soil and pavement layer.
//-----------------------------------------------------------------------------
#define _CRT_SECURE_NO_DEPRECATE

#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "lid.h"
#include "headers.h"
#include "odesolve.h"
#include "rosenbrock.h"
#include "funcs.h"

//-----------------------------------------------------------------------------
//  Constants
//-----------------------------------------------------------------------------
#define STOPTOL  0.00328     // integration error tolerance in ft (= 1 mm)
#define MINFLOW  2.3e-8      // flow cutoff for dry conditions (= 0.001 in/hr)
static const double TREEPITTOL = 0.0001;    // ODE solver tolerance
static const double XTOL  = 0.001;     // tolerance for moisture & depth
static const double EPSILON = 0.000001;  // minimum value for exfiltration

//-----------------------------------------------------------------------------
//  Enumerations
//-----------------------------------------------------------------------------
enum LidLayerTypes {
    SURF,                    // surface layer
    SOIL,                    // soil layer
    STOR,                    // storage layer
    PAVE,                    // pavement layer
    ROOT,                    // rooted soil layer
    DRAIN};                  // underdrain system

enum LidRptVars {
    SURF_INFLOW,             // inflow to surface layer
    TOTAL_EVAP,              // evaporation rate from all layers
    SURF_INFIL,              // infiltration into surface layer
    PAVE_PERC,               // percolation through pavement layer
    SOIL_PERC,               // percolation through soil layer
    ROOT_PERC,               // percolation through root layer
    STOR_EXFIL,              // exfiltration out of storage layer
    SURF_OUTFLOW,            // outflow from surface layer
    STOR_DRAIN,              // outflow from storage layer
    SURF_DEPTH,              // ponded depth on surface layer
    PAVE_DEPTH,              // water level in pavement layer
    SOIL_MOIST,              // moisture content of soil layer
    STOR_DEPTH,              // water level in storage layer
    TOTAL_TRANSP,            // transpiration rate from all layers
    SOIL_MOIST2,             // moisture content of rooted soil layer
    UNROOT_PERC,              // infiltration rate into storage layer
    MAX_RPT_VARS};

//-----------------------------------------------------------------------------
//  Imported variables
//-----------------------------------------------------------------------------
extern char HasWetLids;      // TRUE if any LIDs are wet (declared in runoff.c)

//-----------------------------------------------------------------------------
//  Local Variables
//-----------------------------------------------------------------------------
static TLidUnit*  theLidUnit;     // ptr. to a subcatchment's LID unit
static TLidProc*  theLidProc;     // ptr. to a LID process

static double     Tstep;          // current time step (sec)
static double     EvapRate;       // evaporation rate (ft/s)
static double     MaxNativeInfil; // native soil infil. rate limit (ft/s)

static double     SurfaceInflow;  // precip. + runon to LID unit (ft/s)
static double     SurfaceInfil;   // infil. rate from surface layer (ft/s)
static double     RootedInfil;    // infil. rate from surface into rooted layer (ft/s)
static double     UnrootedInfil;  // infil. rate from surface into unrooted layer (ft/s)
static double     SurfaceEvap;    // evap. rate from surface layer (ft/s)
static double     SurfaceOutflow; // outflow from surface layer (ft/s)
static double     SurfaceVolume;  // volume in surface storage (ft)
static double     MaxInfil;       // maximum infiltration rate (ft/s)
static double     MaxInfilRooted; // maximum rooted infiltration rate (ft/s)
static double     MaxSurfaceExfil;// maximum infiltratin rate from available surface water (ft/s)
static double     MaxSurfaceOutflow;// maximum surface outflow rate from available surface water (ft/s)
static double     MaxInfilUnrooted; // maximum unrooted infiltration rate (ft/s)

static double     PaveEvap;       // evap. from pavement layer (ft/s)
static double     PavePerc;       // percolation from pavement layer (ft/s)
static double     PaveVolume;     // volume stored in pavement layer  (ft)

static double     SoilEvap;       // evap. from soil layer (ft/s)
static double     RootedEvap;     // evap. from rooted soil layer (ft/s)
static double     UnrootedEvap;   // evap. from unrooted soil layer (ft/s)
static double     SoilTransp;     // transpiration rate from rooted soil layer (ft/s)
static double     SoilPerc;       // percolation from soil layer (ft/s)
static double     UnrootedPerc;   // percolation from unrooted soil layer (ft/s)
static double     RootedPerc;     // percolation from rooted soil layer (ft/s)
static double     SoilDrain;      // drainage rate from soil underdrain (ft/s)
static double     SoilVolume;     // volume in soil/pavement storage (ft)
static double     RootedVolume;   // volume in rooted soil storage (ft)
static double     MaxSoilPerc;    // maximum percolation from upper to lower zone (ft/s)
static double     MaxSoilPercUnrooted;  // maximum unrooted percolation from upper to lower zone (ft/s)
static double     MaxSoilPercRooted;    // maximum rooted percolation from upper to lower zone (ft/s)
static double     MaxSoilDrain;   // maximum drainage rate from soil layer (ft/s)
static double     MaxSatExfil;    // maximum exfiltration rate from saturated zone in soil layer (ft/s)

static double     StorageInflow;  // inflow rate to storage layer (ft/s)
static double     StorageExfil;   // exfil. rate from storage layer (ft/s)
static double     StorageEvap;    // evap.rate from storage layer (ft/s)
static double     StorageDrain;   // underdrain flow rate layer (ft/s)
static double     MaxStorageDrain;// maximum drain rate (ft/s)
static double     StorageVolume;  // volume in storage layer (ft)
static double     MaxPavePerc;    // maximum infil. rate into storage layer (ft/s)
static double     MaxStorageExfil;// maximum exfil. rate from storage layer (ft/s)

static double     Xold[MAX_LAYERS];  // previous moisture level in LID layers

//-----------------------------------------------------------------------------
//  External Functions (declared in lid.h)
//-----------------------------------------------------------------------------
// lidproc_initWaterBalance  (called by lid_initState)
// lidproc_initWaterRate     (called by lid_initState, OWA addition)
// lidproc_getOutflow        (called by evalLidUnit in lid.c)
// lidproc_saveResults       (called by evalLidUnit in lid.c)

//-----------------------------------------------------------------------------
// Local Functions
//-----------------------------------------------------------------------------
static void   barrelFluxRates(double x[], double f[]);
static void   biocellFluxRates(double x[], double f[]);
static void   greenRoofFluxRates(double x[], double f[]);
static void   pavementFluxRates(double x[], double f[]);
static void   trenchFluxRates(double x[], double f[]);
static void   swaleFluxRates(double x[], double f[]);
static void   roofFluxRates(double x[], double f[]);
static void   treepitFluxRates(double x[], double f[]);

static double getSurfaceOutflowRate(double depth);
static double getAdaptiveSurfaceOutflowRate(double depth);
static double getSurfaceOverflowRate(double* surfaceDepth);
static double getPavementPermRate(void);
static double getSoilPercRate(double theta);
static double getStorageExfilRate(void);
static double getTreepitExfilRate(double storageDepth);
static double getStorageDrainRate(double storageDepth, double soilTheta,
              double paveDepth, double surfaceDepth);
static double getTreepitDrainRate(double storageDepth, double satDepth, double meanTheta,
                                  double distzoneDepth);
static double getDrainMatOutflow(double depth);
static void   getEvapRates(double surfaceVol, double paveVol,
              double soilVol, double storageVol, double pervFrac);
static double getWaterStressResponse(double theta);
static void   getTreepitFluxes(double distzoneDepth, double soilTheta, double rootedTheta,
                               double satDepth, double storageDepth);

static void   updateWaterBalance(TLidUnit *lidUnit, double inflow,
                                 double evap, double infil, double surfFlow,
                                 double drainFlow, double storage);

// OWA EDIT ##################################################################################
// function to store additional data variables used to compute the water balance of LID Units.
static void updateWaterRate(TLidUnit *lidUnit, double evap, double maxNativeInfil,
                            double surfaceInflow, double surfInfil, double surfaceEvap,
                            double surfaceOutflow, double paveEvap, double pavePerc,
                            double soilEvap, double soilPerc, double storageInflow,
                            double storageExfil, double storageEvap, double storageDrain);
// ###########################################################################################

static int    modpuls_solve(int n, double* x, double* xOld, double* xPrev,
                            double* xMin, double* xMax, double* xTol,
                            double* qOld, double* q, double dt, double omega,
                            void (*derivs)(double*, double*));

//=============================================================================

void lidproc_initWaterBalance(TLidUnit *lidUnit, double initVol)
//
//  Purpose: initializes the water balance components of a LID unit.
//  Input:   lidUnit = a particular LID unit
//           initVol = initial water volume stored in the unit (ft)
//  Output:  none
//
{
    lidUnit->waterBalance.inflow = 0.0;
    lidUnit->waterBalance.evap = 0.0;
    lidUnit->waterBalance.infil = 0.0;
    lidUnit->waterBalance.surfFlow = 0.0;
    lidUnit->waterBalance.drainFlow = 0.0;
    lidUnit->waterBalance.initVol = initVol;
    lidUnit->waterBalance.finalVol = initVol;
}

// OWA EDIT ##################################################################################
// initilize struct that stores additional data variables used to compute the water balance
// of LID Units. These variables are used in EPA SWMM, but are not stored for each LID Unit
void lidproc_initWaterRate(TLidUnit* lidUnit)
//
//  Purpose: initializes the water balance components of a LID unit.
//  Input:   lidUnit = a particular LID unit
//           initVol = initial water volume stored in the unit (ft)
//  Output:  none
//
{
    lidUnit->waterRate.evap = 0.0;
    lidUnit->waterRate.maxNativeInfil = 0.0;
    lidUnit->waterRate.surfaceInflow = 0.0;
    lidUnit->waterRate.surfaceInfil = 0.0;
    lidUnit->waterRate.surfaceEvap = 0.0;
    lidUnit->waterRate.surfaceOutflow = 0.0;
    lidUnit->waterRate.paveEvap = 0.0;
    lidUnit->waterRate.pavePerc = 0.0;
    lidUnit->waterRate.soilEvap = 0.0;
    lidUnit->waterRate.soilPerc = 0.0;
    lidUnit->waterRate.storageInflow = 0.0;
    lidUnit->waterRate.storageExfil = 0.0;
    lidUnit->waterRate.storageEvap = 0.0;
    lidUnit->waterRate.storageDrain = 0.0;
}
// ###########################################################################################

//=============================================================================

double lidproc_getOutflow(TLidUnit* lidUnit, TLidProc* lidProc, double inflow,
                          double evap, double infil, double maxInfil,
                          double tStep, double* lidEvap,
                          double* lidInfil, double* lidDrain)
//
//  Purpose: computes runoff outflow from a single LID unit.
//  Input:   lidUnit  = ptr. to specific LID unit being analyzed
//           lidProc  = ptr. to generic LID process of the LID unit
//           inflow   = runoff rate captured by LID unit (ft/s)
//           evap     = potential evaporation rate (ft/s)
//           infil    = infiltration rate to native soil (ft/s)
//           maxInfil = max. infiltration rate to native soil (ft/s)
//           tStep    = time step (sec)
//  Output:  lidEvap  = evaporation rate for LID unit (ft/s)
//           lidInfil = infiltration rate for LID unit (ft/s)
//           lidDrain = drain flow for LID unit (ft/s)
//           returns surface runoff rate from the LID unit (ft/s)
//
{
    int    i;
    double x[MAX_LAYERS];        // layer moisture levels
    double xOld[MAX_LAYERS];     // work vector
    double xPrev[MAX_LAYERS];    // work vector
    double xMin[MAX_LAYERS];     // lower limit on moisture levels
    double xMax[MAX_LAYERS];     // upper limit on moisture levels
    double fOld[MAX_LAYERS];     // previously computed flux rates
    double f[MAX_LAYERS];        // newly computed flux rates

    // convergence tolerance on moisture levels (ft, moisture fraction , ft)
    double xTol[MAX_LAYERS] = {STOPTOL, STOPTOL, STOPTOL, STOPTOL, STOPTOL};


    double omega = 0.0;          // integration time weighting

    //... define a pointer to function that computes flux rates through the LID
    void (*fluxRates) (double *, double *) = NULL;

    //... save references to the LID process and LID unit
    theLidProc = lidProc;
    theLidUnit = lidUnit;

    //... save evap, max. infil. & time step to shared variables
    EvapRate = evap;
    MaxNativeInfil = maxInfil;
    Tstep = tStep;

    //... store current moisture levels in vector x
    x[SURF] = theLidUnit->surfaceDepth;
    x[SOIL] = theLidUnit->soilMoisture;
    x[STOR] = theLidUnit->storageDepth;
    x[PAVE] = theLidUnit->paveDepth;
    x[ROOT] = theLidUnit->soilMoisture2;

    //... initialize layer moisture volumes, flux rates and moisture limits
    SurfaceVolume  = 0.0;
    PaveVolume     = 0.0;
    SoilVolume     = 0.0;
    RootedVolume   = 0.0;
    StorageVolume  = 0.0;
    SurfaceInflow  = inflow;
    SurfaceInfil   = 0.0;
    RootedInfil    = 0.0;
    SurfaceEvap    = 0.0;
    SurfaceOutflow = 0.0;
    PaveEvap       = 0.0;
    PavePerc       = 0.0;
    SoilEvap       = 0.0;
    SoilPerc       = 0.0;
    SoilTransp     = 0.0;
    SoilDrain      = 0.0;
    RootedPerc     = 0.0;
    StorageInflow  = 0.0;
    StorageExfil   = 0.0;
    StorageEvap    = 0.0;
    StorageDrain   = 0.0;
    for (i = 0; i < MAX_LAYERS; i++)
    {
        f[i] = 0.0;
        fOld[i] = theLidUnit->oldFluxRates[i];
        xMin[i] = 0.0;
        xMax[i] = BIG;
        Xold[i] = x[i];
    }

    //... find Green-Ampt infiltration from surface layer
    if ( theLidProc->lidType == POROUS_PAVEMENT ) SurfaceInfil = 0.0;
    else if ( theLidUnit->soilInfil.Ks > 0.0 )
    {
        SurfaceInfil =
            grnampt_getInfil(&theLidUnit->soilInfil, Tstep,
                             SurfaceInflow, theLidUnit->surfaceDepth,
                             MOD_GREEN_AMPT);
    }
    else SurfaceInfil = infil;

    //... set moisture limits for soil & storage layers
    if ( theLidProc->soil.thickness > 0.0 )
    {
        xMin[SOIL] = theLidProc->soil.wiltPoint;
        xMax[SOIL] = theLidProc->soil.porosity;
    }
    if ( theLidProc->pavement.thickness > 0.0 )
    {
        xMax[PAVE] = theLidProc->pavement.thickness;
    }
    if ( theLidProc->storage.thickness > 0.0 )
    {
//        printf("storage > 0.0: %.2f", theLidProc->storage.thickness);
        xMax[STOR] = theLidProc->storage.thickness;
    }
    if ( theLidProc->lidType == GREEN_ROOF )
    {
        xMax[STOR] = theLidProc->drainMat.thickness;
    }
    if ( theLidProc->lidType == TREEPIT )
    {
        xMin[ROOT] = theLidProc->soil.wiltPoint;
        xMax[ROOT] = theLidProc->soil.porosity;
    }

    //... update moisture levels and flux rates over the time step
    // Condition to use a different solver for TREE_PIT
    if (theLidProc->lidType == TREEPIT)
    {
        // Assuming rungekutta45 is defined elsewhere and ready to be used
        treepitFluxRates(x, f);
    }
    else
    {
        //... determine which flux rate function to use
        switch (theLidProc->lidType)
        {
        case BIO_CELL:
        case RAIN_GARDEN:     fluxRates = &biocellFluxRates;   break;
        case GREEN_ROOF:      fluxRates = &greenRoofFluxRates; break;
        case INFIL_TRENCH:    fluxRates = &trenchFluxRates;    break;
        case POROUS_PAVEMENT: fluxRates = &pavementFluxRates;  break;
        case RAIN_BARREL:     fluxRates = &barrelFluxRates;    break;
        case ROOF_DISCON:     fluxRates = &roofFluxRates;      break;
        case VEG_SWALE:       fluxRates = &swaleFluxRates;
                              omega = 0.5;
                              break;
        default:              return 0.0;
        }
        i = modpuls_solve(MAX_LAYERS, x, xOld, xPrev, xMin, xMax, xTol,
                          fOld, f, tStep, omega, fluxRates);
    }

/** For debugging only ********************************************
    if  (i == 0)
    {
        fprintf(Frpt.file,
        "\n  WARNING 09: integration failed to converge at %s %s",
            theDate, theTime);
        fprintf(Frpt.file,
        "\n              for LID %s placed in subcatchment %s.",
            theLidProc->ID, theSubcatch->ID);
    }
*******************************************************************/

    //... add any surface overflow to surface outflow
    if ( theLidProc->surface.canOverflow || theLidUnit->fullWidth == 0.0 )
    {
        SurfaceOutflow += getSurfaceOverflowRate(&x[SURF]);
    }

    //... save updated results
    theLidUnit->surfaceDepth    = x[SURF];
    theLidUnit->soilMoisture    = x[SOIL];
    theLidUnit->soilMoisture2   = x[ROOT];
    theLidUnit->paveDepth       = x[PAVE];
    theLidUnit->storageDepth    = x[STOR];
    for (i = 0; i < MAX_LAYERS; i++) theLidUnit->oldFluxRates[i] = f[i];

    //... assign values to LID unit evaporation, infiltration & drain flow
    *lidEvap = SurfaceEvap + PaveEvap + SoilEvap + StorageEvap + SoilTransp;
    *lidInfil = StorageExfil;
    *lidDrain = StorageDrain;

    //... return surface outflow (per unit area) from unit
    return SurfaceOutflow;
}

//=============================================================================

void lidproc_saveResults(TLidUnit* lidUnit, double ucfRainfall, double ucfRainDepth, double ucfVolume)
//
//  Purpose: updates the mass balance for an LID unit and saves
//           current flux rates to the LID report file.
//  Input:   lidUnit = ptr. to LID unit
//           ucfRainfall = units conversion factor for rainfall rate
//           ucfDepth = units conversion factor for rainfall depth
//  Output:  none
//
{
    double ucf;                        // units conversion factor
    double totalEvap;                  // total evaporation rate (ft/s)
    double totalVolume;                // total volume stored in LID (ft)
    double rptVars[MAX_RPT_VARS];      // array of reporting variables
    int    isDry = FALSE;              // true if current state of LID is dry
    char   timeStamp[TIME_STAMP_SIZE + 1]; // date/time stamp
    double elapsedHrs;                 // elapsed hours

    //... find total evap. rate and stored volume
    totalEvap = SurfaceEvap + PaveEvap + SoilEvap + StorageEvap + SoilTransp;
    totalVolume = SurfaceVolume + PaveVolume + SoilVolume + StorageVolume;

    //... update mass balance totals
    updateWaterBalance(theLidUnit, SurfaceInflow, totalEvap, StorageExfil,
                       SurfaceOutflow, StorageDrain, totalVolume);

    // OWA EDIT ###############################################################
    //... update water rate structs
    updateWaterRate(theLidUnit, EvapRate, MaxNativeInfil, SurfaceInflow,
                    SurfaceInfil, SurfaceEvap, SurfaceOutflow, PaveEvap,
                    PavePerc, SoilEvap+SoilTransp, SoilPerc, StorageInflow, StorageExfil,
                    StorageEvap, StorageDrain);
    // ########################################################################

    //... check if dry-weather conditions hold
    if ( SurfaceInflow  < MINFLOW &&
         SurfaceOutflow < MINFLOW &&
         StorageDrain   < MINFLOW &&
         StorageExfil   < MINFLOW &&
         totalEvap      < MINFLOW
       ) isDry = TRUE;

    //... update status of HasWetLids
    if ( !isDry ) HasWetLids = TRUE;

    //... write results to LID report file
    if ( lidUnit->rptFile )
    {
        //... convert rate results to original units (in/hr or mm/hr)
        ucf = ucfRainfall;
        rptVars[SURF_INFLOW]    = SurfaceInflow*ucf;
        rptVars[TOTAL_EVAP]     = (totalEvap-SoilTransp)*ucf;
        rptVars[SURF_INFIL]     = SurfaceInfil*ucf;
        rptVars[PAVE_PERC]      = PavePerc*ucf;
        rptVars[SOIL_PERC]      = SoilPerc*ucf;
        rptVars[ROOT_PERC]      = RootedPerc*ucf;
        rptVars[STOR_EXFIL]     = StorageExfil*ucf;
        rptVars[SURF_OUTFLOW]   = SurfaceOutflow*ucf;
        rptVars[STOR_DRAIN]     = StorageDrain*ucf;
        rptVars[TOTAL_TRANSP]   = SoilTransp*ucf;
        rptVars[UNROOT_PERC]     = UnrootedPerc*ucf;

        //... convert storage results to original units (in or mm)
        ucf = ucfRainDepth;
        rptVars[SURF_DEPTH] = theLidUnit->surfaceDepth*ucf;
        rptVars[PAVE_DEPTH] = theLidUnit->paveDepth*ucf;
        rptVars[SOIL_MOIST] = theLidUnit->soilMoisture;
        rptVars[STOR_DEPTH] = theLidUnit->storageDepth*ucf;
        rptVars[SOIL_MOIST2] = theLidUnit->soilMoisture2;

        //... if the current LID state is wet but the previous state was dry
        //    for more than one period then write the saved previous results
        //    to the report file thus marking the end of a dry period
        if ( !isDry && theLidUnit->rptFile->wasDry > 1)
        {
            fprintf(theLidUnit->rptFile->file, "%s",
                theLidUnit->rptFile->results);
        }

        //... write the current results to a string which is saved between
        //    reporting periods
        elapsedHrs = NewRunoffTime / 1000.0 / 3600.0;
        datetime_getTimeStamp(
            M_D_Y, getDateTime(NewRunoffTime), TIME_STAMP_SIZE, timeStamp);
        snprintf(theLidUnit->rptFile->results, sizeof(theLidUnit->rptFile->results),
             "\n%20s\t %8.3f\t %8.3f\t %8.4f\t %8.3f\t %8.3f\t %8.3f\t %8.3f\t %8.3f\t"
             "%8.3f\t %8.3f\t %8.3f\t %8.3f\t %8.3f\t %8.3f\t %9.3f\t %9.3f\t %9.3f",
             timeStamp, elapsedHrs, rptVars[0], rptVars[1], rptVars[2],
             rptVars[3], rptVars[4], rptVars[5], rptVars[6], rptVars[7],
             rptVars[8], rptVars[9], rptVars[10], rptVars[11], rptVars[12],
             rptVars[13], rptVars[14], rptVars[15]);

        //... if the current LID state is dry
        if ( isDry )
        {
            //... if the previous state was wet then write the current
            //    results to file marking the start of a dry period
            if ( theLidUnit->rptFile->wasDry == 0 )
            {
                fprintf(theLidUnit->rptFile->file, "%s",
                    theLidUnit->rptFile->results);
            }

            //... increment the number of successive dry periods
            theLidUnit->rptFile->wasDry++;
        }

        //... if the current LID state is wet
        else
        {
            //... write the current results to the report file
            fprintf(theLidUnit->rptFile->file, "%s",
                theLidUnit->rptFile->results);

            //... re-set the number of successive dry periods to 0
            theLidUnit->rptFile->wasDry = 0;
        }
    }
}

//=============================================================================

void roofFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates for roof disconnection.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double surfaceDepth = x[SURF];

    getEvapRates(surfaceDepth, 0.0, 0.0, 0.0, 1.0);
    SurfaceVolume = surfaceDepth;
    SurfaceInfil = 0.0;
    if ( theLidProc->surface.alpha > 0.0 )
      SurfaceOutflow = getSurfaceOutflowRate(surfaceDepth);
    else getSurfaceOverflowRate(&surfaceDepth);
    StorageDrain = MIN(theLidProc->drain.coeff/UCF(RAINFALL), SurfaceOutflow);
    SurfaceOutflow -= StorageDrain;
    f[SURF] = (SurfaceInflow - SurfaceEvap - StorageDrain - SurfaceOutflow);
}

//=============================================================================

void greenRoofFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates from the layers of a green roof.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // Moisture level variables
    double surfaceDepth;
    double soilTheta;
    double storageDepth;

    // Intermediate variables
    double availVolume;
    double maxRate;

    // Green roof properties
    double soilThickness    = theLidProc->soil.thickness;
    double storageThickness = theLidProc->storage.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double storageVoidFrac  = theLidProc->storage.voidFrac;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double soilWiltPoint    = theLidProc->soil.wiltPoint;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    soilTheta    = x[SOIL];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    SurfaceVolume = surfaceDepth * theLidProc->surface.voidFrac;
    SoilVolume = soilTheta * soilThickness;
    StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = SoilVolume - soilWiltPoint * soilThickness;
    getEvapRates(SurfaceVolume, 0.0, availVolume, StorageVolume, 1.0);
    if ( soilTheta >= soilPorosity ) StorageEvap = 0.0;

    //... soil layer perc rate
    SoilPerc = getSoilPercRate(soilTheta);

    //... limit perc rate by available water
    availVolume = (soilTheta - soilFieldCap) * soilThickness;
    maxRate = MAX(availVolume, 0.0) / Tstep - SoilEvap;
    SoilPerc = MIN(SoilPerc, maxRate);
    SoilPerc = MAX(SoilPerc, 0.0);

    //... storage (drain mat) outflow rate
    StorageExfil = 0.0;
    StorageDrain = getDrainMatOutflow(storageDepth);

    //... unit is full
    if ( soilTheta >= soilPorosity && storageDepth >= storageThickness )
    {
        //... outflow from both layers equals limiting rate
        maxRate = MIN(SoilPerc, StorageDrain);
        SoilPerc = maxRate;
        StorageDrain = maxRate;

        //... adjust inflow rate to soil layer
        SurfaceInfil = MIN(SurfaceInfil, maxRate);
    }

    //... unit not full
    else
    {
        //... limit drainmat outflow by available storage volume
        maxRate = storageDepth * storageVoidFrac / Tstep - StorageEvap;
        if ( storageDepth >= storageThickness ) maxRate += SoilPerc;
        maxRate = MAX(maxRate, 0.0);
        StorageDrain = MIN(StorageDrain, maxRate);

        //... limit soil perc inflow by unused storage volume
        maxRate = (storageThickness - storageDepth) * storageVoidFrac / Tstep +
                  StorageDrain + StorageEvap;
        SoilPerc = MIN(SoilPerc, maxRate);

        //... adjust surface infil. so soil porosity not exceeded
        maxRate = (soilPorosity - soilTheta) * soilThickness / Tstep +
                  SoilPerc + SoilEvap;
        SurfaceInfil = MIN(SurfaceInfil, maxRate);
    }

    // ... find surface outflow rate
    SurfaceOutflow = getSurfaceOutflowRate(surfaceDepth);

    // ... compute overall layer flux rates
    f[SURF] = (SurfaceInflow - SurfaceEvap - SurfaceInfil - SurfaceOutflow) /
              theLidProc->surface.voidFrac;
    f[SOIL] = (SurfaceInfil - SoilEvap - SoilPerc) /
              theLidProc->soil.thickness;
    f[STOR] = (SoilPerc - StorageEvap - StorageDrain) /
              theLidProc->storage.voidFrac;
}

//=============================================================================

void biocellFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates from the layers of a bio-retention cell LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // Moisture level variables
    double surfaceDepth;
    double soilTheta;
    double storageDepth;

    // Intermediate variables
    double availVolume;
    double maxRate;

    // LID layer properties
    double soilThickness    = theLidProc->soil.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double soilWiltPoint    = theLidProc->soil.wiltPoint;
    double storageThickness = theLidProc->storage.thickness;
    double storageVoidFrac  = theLidProc->storage.voidFrac;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    soilTheta    = x[SOIL];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    SurfaceVolume = surfaceDepth * theLidProc->surface.voidFrac;
    SoilVolume    = soilTheta * soilThickness;
    StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = SoilVolume - soilWiltPoint * soilThickness;
    getEvapRates(SurfaceVolume, 0.0, availVolume, StorageVolume, 1.0);
    if ( soilTheta >= soilPorosity ) StorageEvap = 0.0;

    //... soil layer perc rate
    SoilPerc = getSoilPercRate(soilTheta);

    //... limit perc rate by available water
    availVolume =  (soilTheta - soilFieldCap) * soilThickness;
    maxRate = MAX(availVolume, 0.0) / Tstep - SoilEvap;
    SoilPerc = MIN(SoilPerc, maxRate);
    SoilPerc = MAX(SoilPerc, 0.0);

    //... exfiltration rate out of storage layer
    StorageExfil = getStorageExfilRate();

    //... underdrain flow rate
    StorageDrain = 0.0;
    if ( theLidProc->drain.coeff > 0.0 )
    {
        StorageDrain = getStorageDrainRate(storageDepth, soilTheta, 0.0,
                                           surfaceDepth);
    }

    //... special case of no storage layer present
    if ( storageThickness == 0.0 )
    {
        StorageEvap = 0.0;
        maxRate = MIN(SoilPerc, StorageExfil);
        SoilPerc = maxRate;
        StorageExfil = maxRate;

        //... limit surface infil. by unused soil volume
        maxRate = (soilPorosity - soilTheta) * soilThickness / Tstep +
                  SoilPerc + SoilEvap;
        SurfaceInfil = MIN(SurfaceInfil, maxRate);
    }

    else
    {
        //... storage & soil layers are full
        if ( soilTheta >= soilPorosity && storageDepth >= storageThickness )
        {
            //... limiting rate is smaller of soil perc and storage outflow
            maxRate = StorageExfil + StorageDrain;
            if ( SoilPerc < maxRate )
            {
                maxRate = SoilPerc;
                if ( maxRate > StorageExfil ) StorageDrain = maxRate - StorageExfil;
                else
                {
                    StorageExfil = maxRate;
                    StorageDrain = 0.0;
                }
            }
            else SoilPerc = maxRate;

            //... apply limiting rate to surface infil.
            SurfaceInfil = MIN(SurfaceInfil, maxRate);
        }

        //... either layer not full
        else
        {
            //... limit storage exfiltration by available storage volume
            maxRate = SoilPerc - StorageEvap + storageDepth*storageVoidFrac/Tstep;
            StorageExfil = MIN(StorageExfil, maxRate);
            StorageExfil = MAX(StorageExfil, 0.0);

            //... limit underdrain flow by volume above drain offset
            if ( StorageDrain > 0.0 )
            {
                maxRate = -StorageExfil - StorageEvap;
                if ( storageDepth >= storageThickness) maxRate += SoilPerc;
                if ( theLidProc->drain.offset <= storageDepth )
                {
                    maxRate += (storageDepth - theLidProc->drain.offset) *
                               storageVoidFrac/Tstep;
                }
                maxRate = MAX(maxRate, 0.0);
                StorageDrain = MIN(StorageDrain, maxRate);
            }

            //... limit soil perc by unused storage volume
            maxRate = StorageExfil + StorageDrain + StorageEvap +
                      (storageThickness - storageDepth) *
                      storageVoidFrac/Tstep;
            SoilPerc = MIN(SoilPerc, maxRate);

            //... limit surface infil. by unused soil volume
            maxRate = (soilPorosity - soilTheta) * soilThickness / Tstep +
                      SoilPerc + SoilEvap;
            SurfaceInfil = MIN(SurfaceInfil, maxRate);
        }
    }

    //... find surface layer outflow rate
    SurfaceOutflow = getSurfaceOutflowRate(surfaceDepth);

    //... compute overall layer flux rates
    f[SURF] = (SurfaceInflow - SurfaceEvap - SurfaceInfil - SurfaceOutflow) /
              theLidProc->surface.voidFrac;
    f[SOIL] = (SurfaceInfil - SoilEvap - SoilPerc) /
              theLidProc->soil.thickness;
    if ( storageThickness == 0.0 ) f[STOR] = 0.0;
    else f[STOR] = (SoilPerc - StorageEvap - StorageExfil - StorageDrain) /
                   theLidProc->storage.voidFrac;
}

//=============================================================================

void trenchFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates from the layers of an infiltration trench LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // Moisture level variables
    double surfaceDepth;
    double storageDepth;

    // Intermediate variables
    double availVolume;
    double maxRate;

    // Storage layer properties
    double storageThickness = theLidProc->storage.thickness;
    double storageVoidFrac = theLidProc->storage.voidFrac;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    SurfaceVolume = surfaceDepth * theLidProc->surface.voidFrac;
    SoilVolume = 0.0;
    StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = (storageThickness - storageDepth) * storageVoidFrac;
    getEvapRates(SurfaceVolume, 0.0, 0.0, StorageVolume, 1.0);

    //... no storage evap if surface ponded
    if ( surfaceDepth > 0.0 ) StorageEvap = 0.0;

    //... nominal storage inflow
    StorageInflow = SurfaceInflow + SurfaceVolume / Tstep;

    //... exfiltration rate out of storage layer
   StorageExfil = getStorageExfilRate();

    //... underdrain flow rate
    StorageDrain = 0.0;
    if ( theLidProc->drain.coeff > 0.0 )
    {
        StorageDrain = getStorageDrainRate(storageDepth, 0.0, 0.0, surfaceDepth);
    }

    //... limit storage exfiltration by available storage volume
    maxRate = StorageInflow - StorageEvap + storageDepth*storageVoidFrac/Tstep;
    StorageExfil = MIN(StorageExfil, maxRate);
    StorageExfil = MAX(StorageExfil, 0.0);

    //... limit underdrain flow by volume above drain offset
    if ( StorageDrain > 0.0 )
    {
        maxRate = -StorageExfil - StorageEvap;
        if (storageDepth >= storageThickness ) maxRate += StorageInflow;
        if ( theLidProc->drain.offset <= storageDepth )
        {
            maxRate += (storageDepth - theLidProc->drain.offset) *
                       storageVoidFrac/Tstep;
        }
        maxRate = MAX(maxRate, 0.0);
        StorageDrain = MIN(StorageDrain, maxRate);
    }

    //... limit storage inflow to not exceed storage layer capacity
    maxRate = (storageThickness - storageDepth)*storageVoidFrac/Tstep +
              StorageExfil + StorageEvap + StorageDrain;
    StorageInflow = MIN(StorageInflow, maxRate);

    //... equate surface infil to storage inflow
    SurfaceInfil = StorageInflow;

    //... find surface outflow rate
    SurfaceOutflow = getSurfaceOutflowRate(surfaceDepth);

    // ... find net fluxes for each layer
    f[SURF] = (SurfaceInflow - SurfaceEvap - StorageInflow - SurfaceOutflow) /
              theLidProc->surface.voidFrac;
    f[STOR] = (StorageInflow - StorageEvap - StorageExfil - StorageDrain) /
              theLidProc->storage.voidFrac;
    f[SOIL] = 0.0;
}

//=============================================================================

void pavementFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates for the layers of a porous pavement LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    //... Moisture level variables
    double surfaceDepth;
    double paveDepth;
    double soilTheta;
    double storageDepth;

    //... Intermediate variables
    double pervFrac = (1.0 - theLidProc->pavement.impervFrac);
    double storageInflow;    // inflow rate to storage layer (ft/s)
    double availVolume;
    double maxRate;

    //... LID layer properties
    double paveVoidFrac     = theLidProc->pavement.voidFrac * pervFrac;
    double paveThickness    = theLidProc->pavement.thickness;
    double soilThickness    = theLidProc->soil.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double soilWiltPoint    = theLidProc->soil.wiltPoint;
    double storageThickness = theLidProc->storage.thickness;
    double storageVoidFrac  = theLidProc->storage.voidFrac;

    //... retrieve moisture levels from input vector
    surfaceDepth = x[SURF];
    paveDepth    = x[PAVE];
    soilTheta    = x[SOIL];
    storageDepth = x[STOR];

    //... convert moisture levels to volumes
    SurfaceVolume = surfaceDepth * theLidProc->surface.voidFrac;
    PaveVolume = paveDepth * paveVoidFrac;
    SoilVolume = soilTheta * soilThickness;
    StorageVolume = storageDepth * storageVoidFrac;

    //... get ET rates
    availVolume = SoilVolume - soilWiltPoint * soilThickness;
    getEvapRates(SurfaceVolume, PaveVolume, availVolume, StorageVolume,
                 pervFrac);

    //... no storage evap if soil or pavement layer saturated
    if ( paveDepth >= paveThickness ||
       ( soilThickness > 0.0 && soilTheta >= soilPorosity )
       ) StorageEvap = 0.0;

    //... find nominal rate of surface infiltration into pavement layer
    SurfaceInfil = SurfaceInflow + (SurfaceVolume / Tstep);

    //... find perc rate out of pavement layer
    PavePerc = getPavementPermRate() * pervFrac;

    //... surface infiltration can't exceed pavement permeability
    SurfaceInfil = MIN(SurfaceInfil, PavePerc);

    //... limit pavement perc by available water
    maxRate = PaveVolume/Tstep + SurfaceInfil - PaveEvap;
    maxRate = MAX(maxRate, 0.0);
    PavePerc = MIN(PavePerc, maxRate);

    //... find soil layer perc rate
    if ( soilThickness > 0.0 )
    {
        SoilPerc = getSoilPercRate(soilTheta);
        availVolume = (soilTheta - soilFieldCap) * soilThickness;
        maxRate = MAX(availVolume, 0.0) / Tstep - SoilEvap;
        SoilPerc = MIN(SoilPerc, maxRate);
        SoilPerc = MAX(SoilPerc, 0.0);
    }
    else SoilPerc = PavePerc;

    //... exfiltration rate out of storage layer
    StorageExfil = getStorageExfilRate();

    //... underdrain flow rate
    StorageDrain = 0.0;
    if ( theLidProc->drain.coeff > 0.0 )
    {
        StorageDrain = getStorageDrainRate(storageDepth, soilTheta, paveDepth,
                                           surfaceDepth);
    }

    //... check for adjacent saturated layers

    //... no soil layer, pavement & storage layers are full
    if ( soilThickness == 0.0 &&
         storageDepth >= storageThickness &&
         paveDepth >= paveThickness )
    {
        //... pavement outflow can't exceed storage outflow
        maxRate = StorageEvap + StorageDrain + StorageExfil;
        if ( PavePerc > maxRate ) PavePerc = maxRate;

        //... storage outflow can't exceed pavement outflow
        else
        {
            //... use up available exfiltration capacity first
            StorageExfil = MIN(StorageExfil, PavePerc);
            StorageDrain = PavePerc - StorageExfil;
        }

        //... set soil perc to pavement perc
        SoilPerc = PavePerc;

        //... limit surface infil. by pavement perc
        SurfaceInfil = MIN(SurfaceInfil, PavePerc);
    }

    //... pavement, soil & storage layers are full
    else if ( soilThickness > 0 &&
              storageDepth >= storageThickness &&
              soilTheta >= soilPorosity &&
              paveDepth >= paveThickness )
    {
        //... find which layer has limiting flux rate
        maxRate = StorageExfil + StorageDrain;
        if ( SoilPerc < maxRate) maxRate = SoilPerc;
        else maxRate = MIN(maxRate, PavePerc);

        //... use up available storage exfiltration capacity first
        if ( maxRate > StorageExfil ) StorageDrain = maxRate - StorageExfil;
        else
        {
            StorageExfil = maxRate;
            StorageDrain = 0.0;
        }
        SoilPerc = maxRate;
        PavePerc = maxRate;

        //... limit surface infil. by pavement perc
        SurfaceInfil = MIN(SurfaceInfil, PavePerc);
    }

    //... storage & soil layers are full
    else if ( soilThickness > 0.0 &&
              storageDepth >= storageThickness &&
              soilTheta >= soilPorosity )
    {
        //... soil perc can't exceed storage outflow
        maxRate = StorageDrain + StorageExfil;
        if ( SoilPerc > maxRate ) SoilPerc = maxRate;

        //... storage outflow can't exceed soil perc
        else
        {
            //... use up available exfiltration capacity first
            StorageExfil = MIN(StorageExfil, SoilPerc);
            StorageDrain = SoilPerc - StorageExfil;
        }
        PavePerc = MIN(PavePerc, SoilPerc);

        //... limit surface infil. by available pavement volume
        availVolume = (paveThickness - paveDepth) * paveVoidFrac;
        maxRate = availVolume / Tstep + PavePerc + PaveEvap;
        SurfaceInfil = MIN(SurfaceInfil, maxRate);
    }

    //... soil and pavement layers are full
    else if ( soilThickness > 0.0 &&
              paveDepth >= paveThickness &&
              soilTheta >= soilPorosity )
    {
        PavePerc = MIN(PavePerc, SoilPerc);
        SoilPerc = PavePerc;
        SurfaceInfil = MIN(SurfaceInfil,PavePerc);
        maxRate = MAX(StorageVolume / Tstep + SoilPerc - StorageEvap, 0.0);
	    StorageExfil = MIN(StorageExfil, maxRate);
    }

    //... no adjoining layers are full
    else
    {
        //... limit storage exfiltration by available storage volume
        //    (if no soil layer, SoilPerc is same as PavePerc)
        maxRate = SoilPerc - StorageEvap + StorageVolume / Tstep;
        maxRate = MAX(0.0, maxRate);
        StorageExfil = MIN(StorageExfil, maxRate);

        //... limit underdrain flow by volume above drain offset
        if ( StorageDrain > 0.0 )
        {
            maxRate = -StorageExfil - StorageEvap;
            if (storageDepth >= storageThickness ) maxRate += SoilPerc;
            if ( theLidProc->drain.offset <= storageDepth )
            {
                maxRate += (storageDepth - theLidProc->drain.offset) *
                           storageVoidFrac/Tstep;
            }
            maxRate = MAX(maxRate, 0.0);
            StorageDrain = MIN(StorageDrain, maxRate);
        }

        //... limit soil & pavement outflow by unused storage volume
        availVolume = (storageThickness - storageDepth) * storageVoidFrac;
        maxRate = availVolume/Tstep + StorageEvap + StorageDrain + StorageExfil;
        maxRate = MAX(maxRate, 0.0);
        if ( soilThickness > 0.0 )
        {
            SoilPerc = MIN(SoilPerc, maxRate);
            maxRate = (soilPorosity - soilTheta) * soilThickness / Tstep +
                      SoilPerc;
        }
        PavePerc = MIN(PavePerc, maxRate);

        //... limit surface infil. by available pavement volume
        availVolume = (paveThickness - paveDepth) * paveVoidFrac;
        maxRate = availVolume / Tstep + PavePerc + PaveEvap;
        SurfaceInfil = MIN(SurfaceInfil, maxRate);
    }

    //... surface outflow
    SurfaceOutflow = getSurfaceOutflowRate(surfaceDepth);

    //... compute overall layer flux rates
    f[SURF] = SurfaceInflow - SurfaceEvap - SurfaceInfil - SurfaceOutflow;
    f[PAVE] = (SurfaceInfil - PaveEvap - PavePerc) / paveVoidFrac;
    if ( theLidProc->soil.thickness > 0.0)
    {
        f[SOIL] = (PavePerc - SoilEvap - SoilPerc) / soilThickness;
        storageInflow = SoilPerc;
    }
    else
    {
        f[SOIL] = 0.0;
        storageInflow = PavePerc;
        SoilPerc = 0.0;
    }
    f[STOR] = (storageInflow - StorageEvap - StorageExfil - StorageDrain) /
              storageVoidFrac;
}

//=============================================================================

void swaleFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates from a vegetative swale LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double depth;            // depth of surface water in swale (ft)
    double topWidth;         // top width of full swale (ft)
    double botWidth;         // bottom width of swale (ft)
    double length;           // length of swale (ft)
    double surfInflow;       // inflow rate to swale (cfs)
    double surfWidth;        // top width at current water depth (ft)
    double surfArea;         // surface area of current water depth (ft2)
    double flowArea;         // x-section flow area (ft2)
    double lidArea;          // surface area of full swale (ft2)
    double hydRadius;        // hydraulic radius for current depth (ft)
    double slope;            // slope of swale side wall (run/rise)
    double volume;           // swale volume at current water depth (ft3)
    double dVdT;             // change in volume w.r.t. time (cfs)
    double dStore;           // depression storage depth (ft)
    double xDepth;           // depth above depression storage (ft)

    //... retrieve state variable from work vector
    depth = x[SURF];
    depth = MIN(depth, theLidProc->surface.thickness);

    //... depression storage depth
    dStore = 0.0;

    //... get swale's bottom width
    //    (0.5 ft minimum to avoid numerical problems)
    slope = theLidProc->surface.sideSlope;
    topWidth = theLidUnit->fullWidth;
    topWidth = MAX(topWidth, 0.5);
    botWidth = topWidth - 2.0 * slope * theLidProc->surface.thickness;
    if ( botWidth < 0.5 )
    {
        botWidth = 0.5;
        slope = 0.5 * (topWidth - 0.5) / theLidProc->surface.thickness;
    }

    //... swale's length
    lidArea = theLidUnit->area;
    length = lidArea / topWidth;

    //... top width, surface area and flow area of current ponded depth
    surfWidth = botWidth + 2.0 * slope * depth;
    surfArea = length * surfWidth;
    flowArea = (depth * (botWidth + slope * depth)) *
               theLidProc->surface.voidFrac;

    //... wet volume and effective depth
    volume = length * flowArea;

    //... surface inflow into swale (cfs)
    surfInflow = SurfaceInflow * lidArea;

    //... ET rate in cfs
    SurfaceEvap = EvapRate * surfArea;
    SurfaceEvap = MIN(SurfaceEvap, volume/Tstep);

    //... infiltration rate to native soil in cfs
    StorageExfil = SurfaceInfil * surfArea;

    //... no surface outflow if depth below depression storage
    xDepth = depth - dStore;
    if ( xDepth <= ZERO ) SurfaceOutflow = 0.0;

    //... otherwise compute a surface outflow
    else
    {
        //... modify flow area to remove depression storage,
        flowArea -= (dStore * (botWidth + slope * dStore)) *
                     theLidProc->surface.voidFrac;
        if ( flowArea < ZERO ) SurfaceOutflow = 0.0;
        else
        {
            //... compute hydraulic radius
            botWidth = botWidth + 2.0 * dStore * slope;
            hydRadius = botWidth + 2.0 * xDepth * sqrt(1.0 + slope*slope);
            hydRadius = flowArea / hydRadius;

            //... use Manning Eqn. to find outflow rate in cfs
            SurfaceOutflow = theLidProc->surface.alpha * flowArea *
                             pow(hydRadius, 2./3.);
        }
    }

    //... net flux rate (dV/dt) in cfs
    dVdT = surfInflow - SurfaceEvap - StorageExfil - SurfaceOutflow;

    //... when full, any net positive inflow becomes spillage
    if ( depth == theLidProc->surface.thickness && dVdT > 0.0 )
    {
        SurfaceOutflow += dVdT;
        dVdT = 0.0;
    }

    //... convert flux rates to ft/s
    SurfaceEvap /= lidArea;
    StorageExfil /= lidArea;
    SurfaceOutflow /= lidArea;
    f[SURF] = dVdT / surfArea;
    f[SOIL] = 0.0;
    f[STOR] = 0.0;

    //... assign values to layer volumes
    SurfaceVolume = volume / lidArea;
    SoilVolume = 0.0;
    StorageVolume = 0.0;
}

//=============================================================================

void barrelFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates for a rain barrel LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double storageDepth = x[STOR];
    double head;
    double maxValue;

    //... assign values to layer volumes
    SurfaceVolume = 0.0;
    SoilVolume = 0.0;
    StorageVolume = storageDepth;

    //... initialize flows
    SurfaceInfil = 0.0;
    SurfaceOutflow = 0.0;
    StorageDrain = 0.0;

    //... compute outflow if time since last rain exceeds drain delay
    //    (dryTime is updated in lid.evalLidUnit at each time step)
    if ( theLidProc->drain.delay == 0.0 ||
        theLidUnit->dryTime >= theLidProc->drain.delay )
    {
        head = storageDepth - theLidProc->drain.offset;
        if ( head > 0.0 )
        {
            StorageDrain = getStorageDrainRate(storageDepth, 0.0, 0.0, 0.0);
            maxValue = (head/Tstep);
            StorageDrain = MIN(StorageDrain, maxValue);
        }
    }

    //... limit inflow to available storage
    StorageInflow = SurfaceInflow;
    maxValue = (theLidProc->storage.thickness - storageDepth) / Tstep +
        StorageDrain;
    StorageInflow = MIN(StorageInflow, maxValue);
    SurfaceInfil = StorageInflow;

    //... assign values to layer flux rates
    f[SURF] = SurfaceInflow - StorageInflow;
    f[STOR] = StorageInflow - StorageDrain;
    f[SOIL] = 0.0;
}

//=============================================================================

void getTreepitFluxes(double surfaceDepth, double soilTheta, double rootedTheta, double satDepth, double storageDepth)
//
//  Input:   upperVolume = vol. depth of upper zone (ft)
//           upperDepth  = depth of upper zone (ft)
//  Output:  none
//  Purpose: computes water fluxes into/out of upper/lower GW zones.
//
{
    // LID layer properties
    double soilThickness    = theLidProc->soil.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double soilWiltPoint    = theLidProc->soil.wiltPoint;
    double surfaceArea      = theLidUnit->area;
    double fracRooted       = theLidProc->tree.fracRooted;
    double storageVoidFrac  = theLidProc->storage.voidFrac;
    double storageThickness = theLidProc->storage.thickness;

    double fracUnrooted = 1.0 - fracRooted;
    double unsatDepth = soilThickness - satDepth;
    double meanTheta = fracRooted * rootedTheta + fracUnrooted * soilTheta;
    double maxRate;

    //... unrooted perc rate
    UnrootedPerc = getSoilPercRate(soilTheta);
    UnrootedPerc = fracUnrooted * UnrootedPerc;
    //... limit unrooted perc rate by available water
    UnrootedPerc = MIN(UnrootedPerc, MaxSoilPercUnrooted);
    UnrootedPerc = MAX(UnrootedPerc, 0.0);

    //... rooted perc rate
    RootedPerc = getSoilPercRate(rootedTheta);
    RootedPerc = fracRooted * RootedPerc;
    //... limit rooted perc rate by available water
    maxRate = MaxSoilPercRooted;
    RootedPerc = MIN(RootedPerc, maxRate);
    RootedPerc = MAX(RootedPerc, 0.0);

    //... storage infil rate
    PavePerc = theLidProc->soil.kSat;
    //... limit storage infil rate by available water in saturated zone
    PavePerc = MIN(PavePerc, MaxSatExfil);
    PavePerc = MAX(PavePerc, 0.0);

    //... exfiltration rate out of storage layer
    StorageExfil = getTreepitExfilRate();

    //... underdrain flow rate
    SoilDrain = 0.0;
    StorageDrain = 0.0;

    if ( theLidProc->drain.coeff > 0.0 )
    {
        if ( theLidProc->drain.offset < storageThickness ) {
            StorageDrain = getTreepitDrainRate(storageDepth, satDepth,
                                               meanTheta, surfaceDepth);
        } else {
            SoilDrain = getTreepitDrainRate(storageDepth, satDepth,
                                            meanTheta, surfaceDepth);
        }
    }

    //... special case of no IWS (storage) layer present - structural soil + susp. pavement
    if ( storageThickness == 0.0 )
    {
        StorageEvap = 0.0;

        maxRate = MIN(PavePerc, StorageExfil);
        PavePerc = maxRate;
        StorageExfil = maxRate;

        //... limit soildrain by volume above drain
        SoilDrain = MIN(SoilDrain, MaxSoilDrain);
        SoilDrain = MAX(SoilDrain, 0.0);

        //... limit surface infil. by unused soil volume
        RootedInfil = MIN(MaxInfilRooted, fracRooted*SurfaceInfil);
        UnrootedInfil = MIN(MaxInfilUnrooted, fracUnrooted*SurfaceInfil);

        SurfaceInfil = RootedInfil + UnrootedInfil;
        if (SurfaceInfil > MaxInfil) printf("Mehr Infil als Platz");
    }

    //... special case of IWS layer present - other systems
    else
    {
        //... water constraint - limit storage exfiltration by water in storage
        maxRate = MaxStorageExfil + PavePerc;
        StorageExfil = MIN(StorageExfil, maxRate);
        StorageExfil = MAX(StorageExfil, 0.0);

        //... water constraint - limit storage-underdrain flow by volume above drain offset
        maxRate = MaxStorageDrain + PavePerc;
        maxRate = MAX(maxRate, 0.0);
        StorageDrain = MIN(StorageDrain, maxRate);

        //... space constraint - limit store infil by unused storage volume
        maxRate = MaxPavePerc + StorageDrain + StorageExfil + StorageEvap;
        maxRate = MIN(maxRate, MaxSatExfil);
        PavePerc = MIN(PavePerc, maxRate);
        PavePerc = MAX(PavePerc, 0.0);

        //... water constraint - limit soil-underdrain flow by volume above drain offset
        SoilDrain = MIN(SoilDrain, MaxSoilDrain);
        SoilDrain = MAX(SoilDrain, 0.0);

        //... space constraint - limit surface infil. by unused soil volume
        RootedInfil = MIN(MaxInfilRooted, fracRooted*SurfaceInfil);
        UnrootedInfil = MIN(MaxInfilUnrooted, fracUnrooted*SurfaceInfil);

        SurfaceInfil = RootedInfil + UnrootedInfil;
        if (SurfaceInfil > MaxInfil) printf("Mehr Infil als Platz");
    }
    //... find surface layer outflow rate
    SurfaceOutflow = getAdaptiveSurfaceOutflowRate(surfaceDepth);
}

//=============================================================================

void  getTreepitDxDt(double t, double* x, double* dxdt)
//
//  Input:   t    = current time (not used)
//           x    = array of state variables
//  Output:  dxdt = array of time derivatives of state variables
//  Purpose: computes time derivatives of upper moisture content
//           and lower depth.
//
//
//  Purpose: computes flux rates from the layers of a bio-retention cell LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    double fUnsatRooted;    // inflow - outflow for upper rooted zone (ft/sec)
    double fUnsatUnrooted;  // inflow - outflow for upper unrooted zone (ft/sec)
    double fSat;            // inflow - outflow for lower zone (ft/sec)
    double fStorage;        // inflow - outflow of storage layer (ft/sec)
    double denom;

    // LID layer properties
    double soilFieldCap  = theLidProc->soil.fieldCap;
    double soilWiltPoint = theLidProc->soil.wiltPoint;
    double satConduct    = theLidProc->soil.kSat;
    double fracRooted    = theLidProc->tree.fracRooted;

    double fracUnrooted = 1.0 - fracRooted;
    double meanTheta    = fracUnrooted * x[SOIL] + fracRooted * x[ROOT];

    //... calculate fluxes
    getTreepitFluxes(x[SURF], x[SOIL], x[ROOT], x[PAVE], x[STOR]);

    //... calculate net fluxes
    fUnsatRooted = RootedInfil - RootedEvap - SoilTransp - RootedPerc;
    fUnsatUnrooted = UnrootedInfil - UnrootedEvap - UnrootedPerc;

    fSat = RootedPerc + UnrootedPerc - PaveEvap - SoilDrain - PavePerc;
    fStorage = PavePerc - StorageDrain - StorageExfil - StorageEvap;

    dxdt[SURF] = SurfaceInflow - RootedInfil - UnrootedInfil - SurfaceEvap - SurfaceOutflow;

    // --- d(upper zone moisture)/dt = (net upper zone flow) /
    //                                 (upper zone depth)
    denom = theLidProc->soil.thickness - x[PAVE];
    if (denom > 0.0) {
        dxdt[SOIL] = 0.0;
        dxdt[ROOT] = 0.0;
        if ( fracUnrooted > 0.0 ) dxdt[SOIL] = fUnsatUnrooted / (fracUnrooted * denom);
        if ( fracRooted > 0.0 ) dxdt[ROOT] = fUnsatRooted / (fracRooted * denom);
    }
    else {
        dxdt[SOIL] = 0.0;
        dxdt[ROOT] = 0.0;
    }
    // --- d(lower zone depth)/dt = (net lower zone flow) /
    //                              (upper zone moisture deficit)
    denom = theLidProc->soil.porosity - meanTheta;
    if (denom > 0.0)
        dxdt[PAVE] = fSat / denom;
    else
        dxdt[PAVE] = 0.0;

    if ( theLidProc->storage.thickness == 0.0) dxdt[STOR] = 0.0;
    else dxdt[STOR] = fStorage / theLidProc->storage.voidFrac;
}

//=============================================================================

void treepitFluxRates(double x[], double f[])
//
//  Purpose: computes flux rates for a rain barrel LID.
//  Input:   x = vector of storage levels
//  Output:  f = vector of flux rates
//
{
    // LID layer properties
    double soilThickness    = theLidProc->soil.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double soilWiltPoint    = theLidProc->soil.wiltPoint;
    double storageThickness = theLidProc->storage.thickness;
    double storagePorosity  = theLidProc->storage.voidFrac;
    double drainageOffset   = theLidProc->drain.offset;
    int pattern             = theLidProc->tree.LAICurve;
    double lai              = theLidProc->tree.LAI;
    double crownArea        = theLidProc->tree.crownArea;
    double fracRooted       = theLidProc->tree.fracRooted;

    // Intermediate variables
    double vAvail;
    double vUnsat;
    double tempRate;
    double g;
    int month;

    //... retrieve moisture levels from input vector
    double distzoneDepth = x[SURF];  // distribution zone
    double soilTheta     = x[SOIL];  // unrooted SMC rooting zone
    double rootedTheta   = x[ROOT];  // rooted SMC rooting zone
    double satDepth      = x[PAVE];  // saturated depth in rooting zone
    double storageDepth  = x[STOR];  // IWS

    //... calculate helper variables
    double fracUnrooted = 1.0 - fracRooted;
    double meanTheta = fracRooted * rootedTheta + fracUnrooted * soilTheta;
    double unsatDepth = soilThickness - satDepth;

    //... convert moisture levels to volumes
    SurfaceVolume = distzoneDepth * theLidProc->surface.voidFrac;
    SoilVolume    = meanTheta * unsatDepth;
    PaveVolume    = satDepth * soilPorosity;
    StorageVolume = storageDepth * storagePorosity;

    //... get ET rates
    // for order of calculations, paveVol ^= SMC and soilVol ^= sat. zone in soil layer
    vAvail = MAX(SoilVolume - soilWiltPoint * unsatDepth, 0.0);
    getEvapRates(SurfaceVolume, vAvail, PaveVolume, StorageVolume, 1.0);
    if ( soilTheta >= soilPorosity ) StorageEvap = 0.0;
    // switch PaveEvap and SoilEvap for better names
    tempRate = PaveEvap;
    PaveEvap = SoilEvap;
    SoilEvap = tempRate;

    // split SoilEvap into RootedEvap and UnrootedEvap
    vAvail = fracRooted * unsatDepth * (rootedTheta - soilWiltPoint);
    RootedEvap = MIN(fracRooted * SoilEvap, vAvail / Tstep);
    vAvail = fracUnrooted * unsatDepth * (soilTheta - soilWiltPoint);
    UnrootedEvap = MIN(SoilEvap - RootedEvap, vAvail / Tstep);
    UnrootedEvap = MAX(UnrootedEvap, 0);
    SoilEvap = RootedEvap + UnrootedEvap;

    // --- calculate potential maximum transpiration
    // --- apply user-supplied LAI pattern
    g = 1.0;
    if ( pattern >= 0 )
    {
        month = datetime_monthOfYear(getDateTime(NewRunoffTime));
        g = Pattern[pattern].factor[month-1];
    }
    lai *= g;
    // --- transpiration rate from rooted fraction but calculted for entire LID area
    tempRate = EvapRate * lai * crownArea / theLidUnit->area;
    SoilTransp = tempRate * getWaterStressResponse(rootedTheta);
    // --- set limit on transpiration rate from rooted fraction of upper GW zone
    vUnsat = fracRooted * (rootedTheta - soilWiltPoint) * unsatDepth;
    tempRate = vUnsat / Tstep - RootedEvap;
    tempRate = MAX(0.0, tempRate);
    SoilTransp = MIN(SoilTransp, tempRate);

    // --- set limit on SurfaceOutflow based on exceeded ponding depth
    MaxSurfaceOutflow = MAX((distzoneDepth - theLidProc->surface.thickness) / Tstep, 0);
    MaxSurfaceOutflow += SurfaceInflow;
    // --- set limit on SurfaceInfiltration based on available volume on surface
    MaxSurfaceExfil = distzoneDepth / Tstep;
    // --- set limit on SurfaceInfiltration based on available volume in unsaturated Layer
    vUnsat = fracRooted * (soilPorosity - rootedTheta) * unsatDepth;
    MaxInfilRooted = vUnsat / Tstep;

    vUnsat = fracUnrooted * (soilPorosity - soilTheta) * unsatDepth;
    MaxInfilUnrooted = vUnsat / Tstep;

    MaxInfil = MaxInfilUnrooted + MaxInfilRooted;

    // --- set limit on percolation rate from upper to lower rooted GW zone
    vUnsat = fracRooted * unsatDepth * (rootedTheta - soilFieldCap);
    MaxSoilPercRooted = vUnsat / Tstep;
    MaxSoilPercRooted = MAX(0.0, MaxSoilPercRooted);

    // --- set limit on percolation rate from upper to lower unrooted GW zone
    vUnsat = fracUnrooted * unsatDepth * (soilTheta - soilFieldCap);
    MaxSoilPercUnrooted = vUnsat / Tstep;
    MaxSoilPercUnrooted = MAX(0.0, MaxSoilPercUnrooted);

    //... limit storage infil by avail saturated volume
    MaxSatExfil = ( satDepth * soilPorosity ) / Tstep;

    //... limit storage infil by avail storage volume
    MaxPavePerc = ( storageThickness - storageDepth ) * storagePorosity / Tstep;

    //... limit storage exfiltration by available storage volume
    MaxStorageExfil = storageDepth * storagePorosity / Tstep;

    //... limit underdrain flow by volume above drain offset - case underdrain is in storage layer
    if ( theLidProc->drain.offset < storageThickness )
    {
        MaxSoilDrain = 0.0;
        vAvail = MAX(0.0, (storageDepth - theLidProc->drain.offset) * storagePorosity);
        MaxStorageDrain =  vAvail / Tstep;
    }
    //... case underdrain is in soil layer
    else
    {
        MaxStorageDrain = 0.0;
        vAvail = MAX(0.0, (storageThickness + satDepth - theLidProc->drain.offset) * soilPorosity);
        MaxSoilDrain = vAvail / Tstep;
    }

    // --- integrate eqns. for d(Theta)/dt and d(LowerDepth)/dt
//    odesolve_integrate(x, 5, 0, Tstep, TREEPITTOL, Tstep, getTreepitDxDt);
    int status = rosenbrock_integrate(x, 5, 0, Tstep, TREEPITTOL, Tstep, getTreepitDxDt);
    if (status != 0) {
        printf("Integration failed.\n");
        // Handle the error or exit
    }
    // --- keep state variables within allowable bounds
    x[SURF] = MAX(x[SURF], 0.0);
    x[SOIL] = MAX(x[SOIL], soilWiltPoint);
    if ( x[SOIL] >= soilPorosity )
    {
        printf("Unrooted SMC MAX SOLL: %.2f, IST: %.2f, Depth und ThetaUR werden zurueckgesetzt\n", soilPorosity, x[SOIL]);
        x[SOIL] = soilPorosity - XTOL;
        x[PAVE] = soilThickness - XTOL;
    }
    x[ROOT] = MAX(x[ROOT], soilWiltPoint);
    if ( x[ROOT] >= soilPorosity )
    {
        printf("Rooted SMC MAX SOLL: %.2f, IST: %.2f, Depth und ThetaRt werden zurueckgesetzt\n", soilPorosity, x[ROOT]);
        x[ROOT] = soilPorosity - XTOL;
        x[PAVE] = soilThickness - XTOL;
    }
    x[PAVE] = MAX(x[PAVE],  0.0);
    if ( x[PAVE] >= soilThickness )
    {
        x[PAVE] = soilThickness - XTOL;
        printf("SaturatedDepth wird zurückgesetzt\n");
    }
    x[STOR] = MAX(x[STOR],  0.0);
    if ( x[STOR] > storageThickness )
    {
        x[STOR] = storageThickness;
        printf("IWS Depth wird zurückgesetzt\n");
    }
    //... for result write soil drain to storage drain
    StorageDrain += SoilDrain;
    //... for result write SoilPerc as sum of RootedPerc and UnrootedPerc
    SoilPerc = RootedPerc + UnrootedPerc;
}

//=============================================================================

double getSurfaceOutflowRate(double depth)
//
//  Purpose: computes outflow rate from a LID's surface layer.
//  Input:   depth = depth of ponded water on surface layer (ft)
//  Output:  returns outflow from surface layer (ft/s)
//
//  Note: this function should not be applied to swales or rain barrels.
//
{
    double delta;
    double outflow;

    //... no outflow if ponded depth below storage depth
    delta = depth - theLidProc->surface.thickness;
    if ( delta < 0.0 ) return 0.0;

    //... compute outflow from overland flow Manning equation
    outflow = theLidProc->surface.alpha * pow(delta, 5.0/3.0) *
              theLidUnit->fullWidth / theLidUnit->area;
    outflow = MIN(outflow, delta / Tstep);
    return outflow;
}

//=============================================================================

double getAdaptiveSurfaceOutflowRate(double depth)
//
//  Purpose: computes outflow rate from a LID's surface layer.
//  Input:   depth = depth of ponded water on surface layer (ft)
//  Output:  returns outflow from surface layer (ft/s)
//
//  Note: this function should not be applied to swales or rain barrels.
//
{
    double delta;
    double outflow;

    //... no outflow if ponded depth below storage depth
    delta = depth - theLidProc->surface.thickness;
    if ( delta < 0.0 ) return 0.0;

    //... compute outflow from overland flow Manning equation
    outflow = theLidProc->surface.alpha * pow(delta, 5.0/3.0) *
              2 * pow(theLidUnit->area, 0.5) / theLidUnit->area;
    outflow = MIN(outflow, MaxSurfaceOutflow-SurfaceInfil);
    return outflow;
}

//=============================================================================

double getPavementPermRate()
//
//  Purpose: computes reduced permeability of a pavement layer due to
//           clogging.
//  Input:   none
//  Output:  returns the reduced permeability of the pavement layer (ft/s).
//
{
    double permReduction = 0.0;
    double clogFactor= theLidProc->pavement.clogFactor;
    double regenDays = theLidProc->pavement.regenDays;

    // ... find permeability reduction due to clogging
    if ( clogFactor > 0.0 )
    {
        // ... see if permeability regeneration has occurred
        //     (regeneration is assumed to reduce the total
        //      volumetric loading that the pavement has received)
        if ( regenDays > 0.0 )
        {
            if ( OldRunoffTime / 1000.0 / SECperDAY >= theLidUnit->nextRegenDay )
            {
                // ... reduce total volume treated by degree of regeneration
                theLidUnit->volTreated *=
                    (1.0 - theLidProc->pavement.regenDegree);

                // ... update next day that regenration occurs
                theLidUnit->nextRegenDay += regenDays;
            }
        }

        // ... find permeabiity reduction factor
        permReduction = theLidUnit->volTreated / clogFactor;
        permReduction = MIN(permReduction, 1.0);
    }

    // ... return the effective pavement permeability
    return theLidProc->pavement.kSat * (1.0 - permReduction);
}

//=============================================================================

double getSoilPercRate(double theta)
//
//  Purpose: computes percolation rate of water through a LID's soil layer.
//  Input:   theta = moisture content (fraction)
//  Output:  returns percolation rate within soil layer (ft/s)
//
{
    double delta;            // moisture deficit

    // ... no percolation if soil moisture <= field capacity
    if ( theta <= theLidProc->soil.fieldCap ) return 0.0;

    // ... perc rate = unsaturated hydraulic conductivity
    delta = theLidProc->soil.porosity - theta;
    return theLidProc->soil.kSat * exp(-delta * theLidProc->soil.kSlope);

}

//=============================================================================

double getStorageExfilRate()
//
//  Purpose: computes exfiltration rate from storage zone into
//           native soil beneath a LID.
//  Input:   depth = depth of water storage zone (ft)
//  Output:  returns infiltration rate (ft/s)
//
{
    double infil = 0.0;
    double clogFactor = 0.0;

    if ( theLidProc->storage.kSat == 0.0 ) return 0.0;
    if ( MaxNativeInfil == 0.0 ) return 0.0;

    //... reduction due to clogging
    clogFactor = theLidProc->storage.clogFactor;
    if ( clogFactor > 0.0 )
    {
        clogFactor = theLidUnit->waterBalance.inflow / clogFactor;
        clogFactor = MIN(clogFactor, 1.0);
    }

    //... infiltration rate = storage Ksat reduced by any clogging
    infil = theLidProc->storage.kSat * (1.0 - clogFactor);

    //... limit infiltration rate by any groundwater-imposed limit
    return MIN(infil, MaxNativeInfil);
}

//=============================================================================

double getTreepitExfilRate(double storageDepth)
//
//  Purpose: computes exfiltration rate from storage zone into
//           native soil beneath a LID.
//  Input:   depth = depth of water storage zone (ft)
//  Output:  returns infiltration rate (ft/s)
//
{
    double infil = 0.0;
    double clogFactor = 0.0;

    if ( theLidProc->storage.kSat == 0.0 ) return 0.0;
    if ( MaxNativeInfil == 0.0 ) return 0.0;

    //... reduction due to clogging
    clogFactor = theLidProc->storage.clogFactor;
    if ( clogFactor > 0.0 )
    {
        clogFactor = theLidUnit->waterBalance.inflow / clogFactor;
        clogFactor = MIN(clogFactor, 1.0);
    }

    //... infiltration rate = storage Ksat reduced by any clogging
    infil = theLidProc->storage.kSat * (1.0 - clogFactor) * pow(storageDepth, theLidProc->storage.expon);;

    //... limit infiltration rate by any groundwater-imposed limit
    return MIN(infil, MaxNativeInfil);
}

//=============================================================================

double  getStorageDrainRate(double storageDepth, double soilTheta,
                            double paveDepth, double surfaceDepth)
//
//  Purpose: computes underdrain flow rate in a LID's storage layer.
//  Input:   storageDepth = depth of water in storage layer (ft)
//           soilTheta    = moisture content of soil layer
//           paveDepth    = effective depth of water in pavement layer (ft)
//           surfaceDepth = depth of ponded water on surface layer (ft)
//  Output:  returns flow in underdrain (ft/s)
//
//  Note:    drain eqn. is evaluated in user's units.
//  Note:    head on drain is water depth in storage layer plus the
//           layers above it (soil, pavement, and surface in that order)
//           minus the drain outlet offset.
{
    int    curve = theLidProc->drain.qCurve;
    double head = storageDepth;
    double outflow = 0.0;
    double paveThickness    = theLidProc->pavement.thickness;
    double soilThickness    = theLidProc->soil.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double storageThickness = theLidProc->storage.thickness;

    // --- storage layer is full
    if ( storageDepth >= storageThickness )
    {
        // --- a soil layer exists
        if ( soilThickness > 0.0 )
        {
            // --- increase head by fraction of soil layer saturated
            if ( soilTheta > soilFieldCap )
            {
                head += (soilTheta - soilFieldCap) /
                        (soilPorosity - soilFieldCap) * soilThickness;

                // --- soil layer is saturated, increase head by water
                //     depth in layer above it
                if ( soilTheta >= soilPorosity )
                {
                    if ( paveThickness > 0.0 )
                    {
                        head += paveDepth;
                        if ( paveDepth >= paveThickness ) head += surfaceDepth;
                    }
                    else head += surfaceDepth;
                }
            }
        }

        // --- no soil layer so increase head by water level in pavement
        //     layer and possibly surface layer
        else if ( paveThickness > 0.0 )
        {
            head += paveDepth;
            if ( paveDepth >= paveThickness ) head += surfaceDepth;
        }
    }

    // --- no outflow if:
    //     a) no prior outflow and head below open threshold
    //     b) prior outflow and head below closed threshold
    if ( theLidUnit->oldDrainFlow == 0.0 &&
         head <= theLidProc->drain.hOpen ) return 0.0;
    if ( theLidUnit->oldDrainFlow > 0.0 &&
         head <= theLidProc->drain.hClose ) return 0.0;

    // --- make head relative to drain offset
    head -= theLidProc->drain.offset;

    // --- compute drain outflow from underdrain flow equation in user units
    //     (head in inches or mm, flow rate in in/hr or mm/hr)
    if ( head > ZERO )
    {
        // --- convert head to user units
        head *= UCF(RAINDEPTH);

        // --- compute drain outflow in user units
        outflow = theLidProc->drain.coeff *
                  pow(head, theLidProc->drain.expon);

        // --- apply user-supplied control curve to outflow
        if (curve >= 0)  outflow *= table_lookup(&Curve[curve], head);

        // --- convert outflow to ft/s
        outflow /= UCF(RAINFALL);
    }
    return outflow;
}

//=============================================================================

double  getTreepitDrainRate(double storageDepth, double satDepth, double meanTheta,
                            double distzoneDepth)
//
//  Purpose: computes underdrain flow rate in a TREEPIT LID's storage layer.
//  Input:   storageDepth  = depth of water in storage layer (ft)
//           soilTheta     = moisture content of soil layer
//           distzoneDepth = effective depth of water in pavement layer (ft)
//  Output:  returns flow in underdrain (ft/s)
//
//  Note:    drain eqn. is evaluated in user's units.
//  Note:    head on drain is water depth in storage layer plus the
//           layers above it (soil, pavement, and surface in that order)
//           minus the drain outlet offset.
{
    int    curve = theLidProc->drain.qCurve;
    double head;
    double outflow = 0.0;
    double distzoneThickness= theLidProc->surface.thickness;
    double soilThickness    = theLidProc->soil.thickness;
    double soilPorosity     = theLidProc->soil.porosity;
    double soilFieldCap     = theLidProc->soil.fieldCap;
    double storageThickness = theLidProc->storage.thickness;
    double drainOffset      = theLidProc->drain.offset;

    //... special case drain is located in storage (IWS) layer
    if (drainOffset < storageThickness)
    {
        head = storageDepth;

        // --- storage layer is full
        if ( storageDepth >= storageThickness )
        {
            head += satDepth;

            // --- soil layer is full
            if ( satDepth >= soilThickness ) head += distzoneDepth;
            // --- saturated zone ends within soil layer
//            else {
//                // --- increase head by fraction of soil layer saturated
//                if ( meanTheta > soilFieldCap )
//                {
//                    head += (meanTheta - soilFieldCap) /
//                            (soilPorosity - soilFieldCap) * (soilThickness - satDepth);
//                }
//            }
        }
    }
    //... underdrain is located in soil layer
    else
    {
        head = storageThickness + satDepth;
        // --- soil layer is full
        if ( satDepth >= soilThickness ) head += distzoneDepth;
            // --- saturated zone ends within soil layer
//        else {
//            // --- increase head by fraction of soil layer saturated
//            if ( meanTheta > soilFieldCap )
//            {
//                head += ((meanTheta - soilFieldCap) /
//                        (soilPorosity - soilFieldCap)) * (soilThickness - satDepth);
//            }
//        }
    }
    // --- no outflow if:
    //     a) no prior outflow and head below open threshold
    //     b) prior outflow and head below closed threshold
    if ( theLidUnit->oldDrainFlow == 0.0 &&
         head <= theLidProc->drain.hOpen ) return 0.0;
    if ( theLidUnit->oldDrainFlow > 0.0 &&
         head <= theLidProc->drain.hClose ) return 0.0;

    // --- make head relative to drain offset
    head -= theLidProc->drain.offset;

    // --- compute drain outflow from underdrain flow equation in user units
    //     (head in inches or mm, flow rate in in/hr or mm/hr)
    if ( head > ZERO )
    {
        // --- convert head to user units
        head *= UCF(RAINDEPTH);

        // --- compute drain outflow in user units
        outflow = theLidProc->drain.coeff *
                  pow(head, theLidProc->drain.expon);

        // --- apply user-supplied control curve to outflow
        if (curve >= 0)  outflow *= table_lookup(&Curve[curve], head);

        // --- convert outflow to ft/s
        outflow /= UCF(RAINFALL);
    }
    return outflow;
}

//=============================================================================

double getDrainMatOutflow(double depth)
//
//  Purpose: computes flow rate through a green roof's drainage mat.
//  Input:   depth = depth of water in drainage mat (ft)
//  Output:  returns flow in drainage mat (ft/s)
//
{
    //... default is to pass all inflow
    double result = SoilPerc;

    //... otherwise use Manning eqn. if its parameters were supplied
    if ( theLidProc->drainMat.alpha > 0.0 )
    {
        result = theLidProc->drainMat.alpha * pow(depth, 5.0/3.0) *
                 theLidUnit->fullWidth / theLidUnit->area *
                 theLidProc->drainMat.voidFrac;
    }
    return result;
}

//=============================================================================

void getEvapRates(double surfaceVol, double paveVol, double soilVol,
    double storageVol, double pervFrac)
//
//  Purpose: computes surface, pavement, soil, and storage evaporation rates.
//  Input:   surfaceVol = volume/area of ponded water on surface layer (ft)
//           paveVol    = volume/area of water in pavement pores (ft)
//           soilVol    = volume/area of water in soil (or pavement) pores (ft)
//           storageVol = volume/area of water in storage layer (ft)
//           pervFrac   = fraction of surface layer that is pervious
//  Output:  none
//
{
    double availEvap;

    //... surface evaporation flux
    availEvap = EvapRate;
    SurfaceEvap = MIN(availEvap, surfaceVol/Tstep);
    SurfaceEvap = MAX(0.0, SurfaceEvap);
    availEvap = MAX(0.0, (availEvap - SurfaceEvap));
    availEvap *= pervFrac;

    //... no subsurface evap if water is infiltrating
    if ( SurfaceInfil > 0.0 )
    {
        PaveEvap = 0.0;
        SoilEvap = 0.0;
        StorageEvap = 0.0;
    }
    else
    {
        //... pavement evaporation flux
        PaveEvap = MIN(availEvap, paveVol / Tstep);
        availEvap = MAX(0.0, (availEvap - PaveEvap));

        //... soil evaporation flux
        SoilEvap = MIN(availEvap, soilVol / Tstep);
        availEvap = MAX(0.0, (availEvap - SoilEvap));

        //... storage evaporation flux
        StorageEvap = MIN(availEvap, storageVol / Tstep);
    }
}

//=============================================================================

double getWaterStressResponse(double theta)
{
    double h1 = theLidProc->soil.porosity;
    double h2 = theLidProc->tree.h2;
    double h3 = theLidProc->tree.h3;
    double h4 = theLidProc->soil.wiltPoint;

    if (theta <= h4) {
        return 0.0;
    } else if (theta < h1 && theta > h2) {
        return table_interpolate(theta, h2, 1.0, h1, 0.0);
    } else if (theta <= h2 && theta >= h3) {
        return 1.0;
    } else if (theta < h3 && theta > h4) {
        return table_interpolate(theta, h4, 0.0, h3, 1.0);
    } else return 0.0; // (theta >= h1)
}

//=============================================================================

double getSurfaceOverflowRate(double* surfaceDepth)
//
//  Purpose: finds surface overflow rate from a LID unit.
//  Input:   surfaceDepth = depth of water stored in surface layer (ft)
//  Output:  returns the overflow rate (ft/s)
//
{
    double delta = *surfaceDepth - theLidProc->surface.thickness;
    if (  delta <= 0.0 ) return 0.0;
    *surfaceDepth = theLidProc->surface.thickness;
    return delta * theLidProc->surface.voidFrac / Tstep;
}

//=============================================================================

void updateWaterBalance(TLidUnit *lidUnit, double inflow, double evap,
    double infil, double surfFlow, double drainFlow, double storage)
//
//  Purpose: updates components of the water mass balance for a LID unit
//           over the current time step.
//  Input:   lidUnit   = a particular LID unit
//           inflow    = runon + rainfall to the LID unit (ft/s)
//           evap      = evaporation rate from the unit (ft/s)
//           infil     = infiltration out the bottom of the unit (ft/s)
//           surfFlow  = surface runoff from the unit (ft/s)
//           drainFlow = underdrain flow from the unit
//           storage   = volume of water stored in the unit (ft)
//  Output:  none
//
{
    lidUnit->volTreated += inflow * Tstep;
    lidUnit->waterBalance.inflow += inflow * Tstep;
    lidUnit->waterBalance.evap += evap * Tstep;
    lidUnit->waterBalance.infil += infil * Tstep;
    lidUnit->waterBalance.surfFlow += surfFlow * Tstep;
    lidUnit->waterBalance.drainFlow += drainFlow * Tstep;
    lidUnit->waterBalance.finalVol = storage;
}
// OWA EDIT ##################################################################################
void updateWaterRate(TLidUnit *lidUnit, double evap, double maxNativeInfil,
    double surfaceInflow, double surfaceInfil, double surfaceEvap, double surfaceOutflow,
    double paveEvap, double pavePerc, double soilEvap, double soilPerc,
    double storageInflow, double storageExfil, double storageEvap, double storageDrain)
    //
    //  Purpose: updates components of the water rate for a LID unit
    //           over the current time step
    //  Input:   lidUnit   = a particular LID unit
    //           evap              = evaporation rate (ft/s)
    //           maxNativeInfil    = native soil infil. rate limit (ft/s)
    //           surfaceInflow     = precip. + runon to LID unit (ft/s)
    //           surfaceInfil      = infil. rate from surface layer (ft/s)
    //           surfaceEvap       = evap. rate from surface layer (ft/s
    //           surfaceOutflow    = outflow from surface layer (ft/s)
    //           paveEvap          = evap. from pavement layer (ft/s)
    //           pavePerc          = percolation from pavement layer (ft/s)
    //           soilEvap          = evap. from soil layer (ft/s)
    //           soilPerc          = percolation from soil layer (ft/s)
    //           storageInflow     = inflow rate to storage layer (ft/s)
    //           storageExfil      = exfil. rate from storage layer (ft/s)
    //           storageEvap       = evap.rate from storage layer (ft/s)
    //           storageDrain      = underdrain flow rate layer (ft/s)
    //  Output:  none
    //
{
    lidUnit->waterRate.evap = evap;
    lidUnit->waterRate.maxNativeInfil = maxNativeInfil;
    lidUnit->waterRate.surfaceInflow = surfaceInflow;
    lidUnit->waterRate.surfaceInfil = surfaceInfil;
    lidUnit->waterRate.surfaceEvap = surfaceEvap;
    lidUnit->waterRate.surfaceOutflow = surfaceOutflow;
    lidUnit->waterRate.paveEvap = paveEvap;
    lidUnit->waterRate.pavePerc = pavePerc;
    lidUnit->waterRate.soilEvap = soilEvap;
    lidUnit->waterRate.soilPerc = soilPerc;
    lidUnit->waterRate.storageInflow = storageInflow;
    lidUnit->waterRate.storageExfil = storageExfil;
    lidUnit->waterRate.storageEvap = storageEvap;
    lidUnit->waterRate.storageDrain = storageDrain;
}
// ###########################################################################################

//=============================================================================

int modpuls_solve(int n, double* x, double* xOld, double* xPrev,
                  double* xMin, double* xMax, double* xTol,
                  double* qOld, double* q, double dt, double omega,
                  void (*derivs)(double*, double*))
//
//  Purpose: solves system of equations dx/dt = q(x) for x at end of time step
//           dt using a modified Puls method.
//  Input:   n = number of state variables
//           x = vector of state variables
//           xOld = state variable values at start of time step
//           xPrev = state variable values from previous iteration
//           xMin = lower limits on state variables
//           xMax = upper limits on state variables
//           xTol = convergence tolerances on state variables
//           qOld = flux rates at start of time step
//           q = flux rates at end of time step
//           dt = time step (sec)
//           omega = time weighting parameter (use 0 for Euler method
//                   or 0.5 for modified Puls method)
//           derivs = pointer to function that computes flux rates q as a
//                    function of state variables x
//  Output:  returns number of steps required for convergence (or 0 if
//           process doesn't converge)
//
{
    int i;
    int canStop;
    int steps = 1;
    int maxSteps = 20;

    //... initialize state variable values
    for (i=0; i<n; i++)
    {
        xOld[i] = x[i];
        xPrev[i] = x[i];
    }

    //... repeat until convergence achieved
    while (steps < maxSteps)
    {
        //... compute flux rates for current state levels
        canStop = 1;
        derivs(x, q);

        //... update state levels based on current flux rates
        for (i=0; i<n; i++)
        {
            x[i] = xOld[i] + (omega*qOld[i] + (1.0 - omega)*q[i]) * dt;
            x[i] = MIN(x[i], xMax[i]);
            x[i] = MAX(x[i], xMin[i]);

            if ( omega > 0.0 &&
                 fabs(x[i] - xPrev[i]) > xTol[i] ) canStop = 0;
            xPrev[i] = x[i];
        }

        //... return if process converges
        if (canStop) return steps;
        steps++;
    }

    //... no convergence so return 0
    return 0;
}
