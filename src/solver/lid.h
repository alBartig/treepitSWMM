//-----------------------------------------------------------------------------
//   lid.h
//
//   Project: EPA SWMM5
//   Version: 5.2
//   Date:    11/01/21   (Build 5.2.0)
//   Author:  L. Rossman
//
//   Public interface for LID functions.
//
//   Update History
//   ==============
//   Build 5.1.008:
//   - Support added for Roof Disconnection LID.
//   - Support added for separate routing of LID drain flows.
//   - Detailed LID reporting modified.
//   Build 5.1.011:
//   - Water depth replaces moisture content for LID's pavement layer. 
//   - Arguments for lidproc_saveResults() modified.
//   Build 5.1.012:
//   - Redefined meaning of wasDry in TLidRptFile structure.
//   Build 5.1.013:
//   - New member fromPerv added to TLidUnit structure to allow LID
//     units to also treat pervious area runoff.
//   - New members hOpen and hClose addded to TDrainLayer to open/close
//     drain when certain heads are reached.
//   - New member qCurve added to TDrainLayer to allow underdrain flow to
//     be adjusted by a curve of multiplier v. head.
//   - New array drainRmvl added to TLidProc to allow for underdrain
//     pollutant removal values.
//   - New members added to TPavementLayer and TLidUnit to support
//     unclogging permeable pavement at fixed intervals.
//   Build 5.2.0:
//   - Covered property added to RAIN_BARREL parameters
//-----------------------------------------------------------------------------

#ifndef LID_H
#define LID_H

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "infil.h"
#include "objects.h"

//-----------------------------------------------------------------------------
//  Enumerations
//-----------------------------------------------------------------------------
enum LidTypes {
    BIO_CELL,                // bio-retention cell
    RAIN_GARDEN,             // rain garden 
    GREEN_ROOF,              // green roof 
    INFIL_TRENCH,            // infiltration trench
    POROUS_PAVEMENT,         // porous pavement
    RAIN_BARREL,             // rain barrel
    VEG_SWALE,               // vegetative swale
    ROOF_DISCON,             // roof disconnection
    TREEPIT};                // tree pit

enum TimePeriod {
    PREVIOUS,                // previous time period
    CURRENT};                // current time period

//-----------------------------------------------------------------------------
//  Data Structures
//-----------------------------------------------------------------------------
#define MAX_LAYERS 5

// LID Surface Layer
typedef struct
{
    double    thickness;          // depression storage or berm ht. (ft)
    double    voidFrac;           // available fraction of storage volume
    double    roughness;          // surface Mannings n 
    double    surfSlope;          // land surface slope (fraction)
    double    sideSlope;          // swale side slope (run/rise)
    double    alpha;              // slope/roughness term in Manning eqn.
    char      canOverflow;        // 1 if immediate outflow of excess water
}  TSurfaceLayer;

// LID Pavement Layer
typedef struct
{
    double   thickness;           // layer thickness (ft)
    double   voidFrac;            // void volume / total volume
    double   impervFrac;          // impervious area fraction
    double   kSat;                // permeability (ft/sec)
    double   clogFactor;          // clogging factor
    double   regenDays;           // clogging regeneration interval (days)
    double   regenDegree;         // degree of clogging regeneration 
}  TPavementLayer;

// LID Soil Layer
typedef struct
{
    double    thickness;          // layer thickness (ft)
    double    porosity;           // void volume / total volume
    double    fieldCap;           // field capacity
    double    wiltPoint;          // wilting point
    double    suction;            // suction head at wetting front (ft)
    double    kSat;               // saturated hydraulic conductivity (ft/sec)
    double    kSlope;             // slope of log(K) v. moisture content curve
}  TSoilLayer;

// LID Storage Layer
typedef struct
{
    double    thickness;          // layer thickness (ft)
    double    voidFrac;           // void volume / total volume
    double    kSat;               // saturated hydraulic conductivity (ft/sec)
    double    clogFactor;         // clogging factor
    int       covered;            // TRUE if rain barrel is covered
}  TStorageLayer;

// Underdrain System (part of Storage Layer)
typedef struct
{
    double    coeff;              // underdrain flow coeff. (in/hr or mm/hr)
    double    expon;              // underdrain head exponent (for in or mm)
    double    offset;             // offset height of underdrain (ft)
    double    delay;              // rain barrel drain delay time (sec)
    double    hOpen;              // head when drain opens (ft)
    double    hClose;             // head when drain closes (ft)
    int       qCurve;             // curve controlling flow rate (optional)
}  TDrainLayer;

// Drainage Mat Layer (for green roofs)
typedef struct
{
    double    thickness;          // layer thickness (ft)
    double    voidFrac;           // void volume / total volume
    double    roughness;          // Mannings n for green roof drainage mats
    double    alpha;              // slope/roughness term in Manning equation
}  TDrainMatLayer;

// Distribution Pipe Layer (for treepits)
typedef struct
{
    double    diameter;           // diameter of dist. pipe
    double    length;             // length of dist. pipe
    double    coeff;              // dist pipe flow coeff. (in/hr or mm/hr)
    double    expon;              // dist pipe head exponent (for in or mm)
    double    offset;             // offset height of dist. pipe (ft)
    double    clogFactor;         // clogging factor
    int       qCurve;             // curve controlling flow rate (optional)
    int       type;               // 24, number indicating crossectional shape
    TXsect    xsect;              // crosssection object of distpipe
    double    yFull;              // water depth of filled pipe
    double    aFull;              // area of filled pipe
    double    vFull;              // volume of filled pipe
}  TDistPipeLayer;

// Tree Layer (for treepits)
typedef struct
{
    double      h2;                 // h2 for water stress response
    double      h3;                 // h3 for water stress response
    double      LAI;                // Leaf Area Index
    double      crownArea;          // projected crown area
    double      fracRooted;         // fraction of the soil volume containing roots
    int         LAICurve;           // curve controlling LAI (optional)
}  TTreeLayer;

// LID Process - generic LID design per unit of area
typedef struct
{
    char*          ID;            // identifying name
    int            lidType;       // type of LID
    TSurfaceLayer  surface;       // surface layer parameters
    TPavementLayer pavement;      // pavement layer parameters
    TSoilLayer     soil;          // soil layer parameters
    TStorageLayer  storage;       // storage layer parameters
    TDrainLayer    drain;         // underdrain system parameters
    TDrainMatLayer drainMat;      // drainage mat layer
    TDistPipeLayer distPipe;      // dist. pipe layer
    TTreeLayer     tree;          // tree layer
    double*        drainRmvl;     // underdrain pollutant removals
}  TLidProc;

// Water Balance Statistics
typedef struct
{
    double         inflow;        // total inflow (ft)
    double         evap;          // total evaporation (ft)
    double         infil;         // total infiltration (ft)
    double         surfFlow;      // total surface runoff (ft)
    double         drainFlow;     // total underdrain flow (ft)
    double         initVol;       // initial stored volume (ft)
    double         finalVol;      // final stored volume (ft)
}  TWaterBalance;

// OWA EDIT ##################################################################################
// OWA SWMM exposes additional data variables used to compute the water balance of LID Units.
// Those variables may be found in lidproc_getOutflow in lidproc.c and are stored in the waterRate
// prop of lidUnits
typedef struct
{
    double         evap;           // evaporation rate (ft/s)
    double         maxNativeInfil; // native soil infil. rate limit (ft/s)
    double         surfaceInflow;  // percip. + runon to LID unit (ft/s)
    double         surfaceInfil;   // infil. rate from surface layer (ft/s)
    double         surfaceEvap;    // evap. rate from surface layer (ft/s)
    double         surfaceOutflow; // outflow from surface layer (ft/s)
    double         paveEvap;       // evap. from pavement layer (ft/s)
    double         pavePerc;       // percolation from pavement layer (ft/s)
    double         soilEvap;       // evap. from soil layer (ft/s)
    double         soilPerc;       // percolation from soil layer (ft/s)
    double         storageInflow;  // inflow rate to storage layer (ft/s)
    double         storageExfil;   // exfil. rate from storage layer (ft/s)
    double         storageEvap;    // evap. rate from storage layer (ft/s)
    double         storageDrain;   // underdrain flow rate layer (ft/s)
} TWaterRate;
// ###########################################################################################

// LID Report File
typedef struct
{
    FILE*     file;               // file pointer
    int       wasDry;             // number of successive dry periods
    char      results[256];       // results for current time period
}   TLidRptFile;

// LID Unit - specific LID process applied over a given area
typedef struct
{
    int      lidIndex;       // index of LID process
    int      number;         // number of replicate units
    double   area;           // area of single replicate unit (ft2)
    double   fullWidth;      // full top width of single unit (ft)
    double   botWidth;       // bottom width of single unit (ft)
    double   initSat;        // initial saturation of soil & storage layers
    double   fromImperv;     // fraction of impervious area runoff treated
    double   fromPerv;       // fraction of pervious area runoff treated
    int      toPerv;         // 1 if outflow sent to pervious area; 0 if not
    int      drainSubcatch;  // subcatchment receiving drain flow
    int      drainNode;      // node receiving drain flow
    TLidRptFile* rptFile;    // pointer to detailed report file

    TGrnAmpt soilInfil;      // infil. object for biocell soil layer
    TGrnAmpt rootedInfil;    // infil. object for biocell tree layer
    double   surfaceDepth;   // depth of ponded water on surface layer (ft)
    double   paveDepth;      // depth of water in porous pavement layer
    double   soilMoisture;   // moisture content of biocell soil layer
    double   soilMoisture2;   // moisture content of rooted zone in soil layer
    double   storageDepth;   // depth of water in storage layer (ft)
    double   distpipeVol;    // volume stored in distribution pipe

    // net inflow - outflow from previous time step for each LID layer (ft/s)
    double   oldFluxRates[MAX_LAYERS];
                                     
    double   dryTime;        // time since last rainfall (sec)
    double   oldDrainFlow;   // previous drain flow (cfs)
    double   newDrainFlow;   // current drain flow (cfs)
    double   volTreated;     // total volume treated (ft)
    double   nextRegenDay;   // next day when unit regenerated
    TWaterBalance  waterBalance;     // water balance quantites
    TWaterRate     waterRate;       // OWA Addition - water rate within lid layers
}  TLidUnit;

// OWA EDIT ##################################################################################
// LidList and LidGroup struct defs moved to lid.h from lid.c to be shared by toolkit.c

// LID List - list of LID units contained in an LID group
struct LidList
{
    TLidUnit*        lidUnit;     // ptr. to a LID unit
    struct LidList*  nextLidUnit;
};
typedef struct LidList TLidList;

// LID Group - collection of LID units applied to a specific subcatchment
////  Re-defined for release 5.1.008. ////                                     //(5.1.008)
struct LidGroup
{
    double         pervArea;      // amount of pervious area in group (ft2)
    double         flowToPerv;    // total flow sent to pervious area (cfs)
    double         oldDrainFlow;  // total drain flow in previous period (cfs)
    double         newDrainFlow;  // total drain flow in current period (cfs)
    TLidList*      lidList;       // list of LID units in the group
};
typedef struct LidGroup* TLidGroup;
// ###########################################################################################
//-----------------------------------------------------------------------------
//   LID Methods
//-----------------------------------------------------------------------------
void     lid_create(int lidCount, int subcatchCount);
void     lid_delete(void);

int      lid_readProcParams(char* tok[], int ntoks);
int      lid_readGroupParams(char* tok[], int ntoks);

void     lid_validate(void);
void     lid_initState(void);
void     lid_setOldGroupState(int subcatch);

double   lid_getPervArea(int subcatch);
double   lid_getFlowToPerv(int subcatch);
double   lid_getDrainFlow(int subcatch, int timePeriod);
double   lid_getStoredVolume(int subcatch);
void     lid_addDrainLoads(int subcatch, double c[], double tStep);
void     lid_addDrainRunon(int subcatch);
void     lid_addDrainInflow(int subcatch, double f);
void     lid_getRunoff(int subcatch, double tStep);
void     lid_writeSummary(void);
void     lid_writeWaterBalance(void);

// OWA EDIT #########################################################
// additional setter and getter functions for toolkit api
int         lid_getLidUnitCount(int index);
TLidUnit*   lid_getLidUnit(int index, int lidIndex, int* errcode);
TLidProc*   lid_getLidProc(int index);
TLidGroup   lid_getLidGroup(int index);
void        lid_validateLidProc(int index);
void        lid_validateLidGroup(int index);
void        lid_updateLidUnit(TLidUnit* lidUnit, int subIndex);
void        lid_updateAllLidUnit(int lidIndex);
void        lid_updateLidGroup(int index);
void        lidproc_initWaterRate(TLidUnit *lidUnit);
// ##################################################################
//-----------------------------------------------------------------------------

void     lidproc_initWaterBalance(TLidUnit *lidUnit, double initVol);

double   lidproc_getOutflow(TLidUnit* lidUnit, TLidProc* lidProc,
         double inflow, double evap, double infil, double maxInfil,
         double tStep, double* lidEvap, double* lidInfil, double* lidDrain);

void     lidproc_saveResults(TLidUnit* lidUnit, double ucfRainfall,
         double ucfRainDepth, double ucfVolume);

#endif
