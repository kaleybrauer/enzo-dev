/***********************************************************************
/
/  CONVERT SHINE PARTICLE INTO RADIATIVE TRANSFER PARTICLE
/
/  written by: John Wise
/  date:       November, 2005
/  modified1:
/
/ PURPOSE: This routine converts particles that shone by a 1/r^2 law
/          into particles that utilize an adaptive 3D ray tracing
/          scheme.
/
************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"
#include "Hierarchy.h"
#include "TopGridData.h"
#include "LevelHierarchy.h"
#include "phys_constants.h"

#define NO_MEAN_ENERGY

int GetUnits(float *DensityUnits, float *LengthUnits,
	     float *TemperatureUnits, float *TimeUnits,
	     float *VelocityUnits, FLOAT Time);

RadiationSourceEntry* DeleteRadiationSource(RadiationSourceEntry *RS);


int StarParticleRadTransfer(LevelHierarchyEntry *LevelArray[], int level,
			    Star *AllStars)
{

  /* If photon test simulation, don't change the radiation sources. */
  if (ProblemType == 50 || ProblemType == 27)
    return SUCCESS;

  int i, j, nShine, nbins;
  double Q[MAX_ENERGY_BINS], QTotal;
  float qFrac, lramp, energies[MAX_ENERGY_BINS];
  float XRayLuminosityFraction = 0.43;
  Star *cstar;

  for (i = 0; i < MAX_ENERGY_BINS; i++){
    Q[i] = 0.0;
    energies[i] = 0.0;
  }

  /* If sources exist, delete them */

  RadiationSourceEntry *dummy;
  if (GlobalRadiationSources!=NULL){
    dummy = GlobalRadiationSources->NextSource;
    while (dummy != NULL) {
      dummy = DeleteRadiationSource(dummy);
    }
    delete GlobalRadiationSources;
    GlobalRadiationSources = NULL;
  }

  GlobalRadiationSources = new RadiationSourceEntry;
  GlobalRadiationSources->NextSource = NULL;
  GlobalRadiationSources->PreviousSource = NULL;

  if (AllStars == NULL)
    return SUCCESS;

  /* Retrieve the units */

  //FLOAT Time = LevelArray[level]->GridData->ReturnTime();
  float DensityUnits, LengthUnits, TemperatureUnits, TimeUnits,
    VelocityUnits;
  GetUnits(&DensityUnits, &LengthUnits, &TemperatureUnits,
	   &TimeUnits, &VelocityUnits, PhotonTime);

  // Convert from #/s to RT units
  double LConv = (double) TimeUnits / pow(LengthUnits,3);

  // Convert to years
  float TimeInYears = yr_s / TimeUnits;

  // bins for computing the luminosity properties (IR, FUV, LW)
  // for rad sources. this is used in OT methods.
  //  making as named variables here for clarity
  const int IRbin=4, FUVbin=7, LWbin=3;

  for (cstar = AllStars; cstar; cstar = cstar->NextStar) {

    // Check the rules if this star particle is radiative
    if (cstar->IsARadiationSource(PhotonTime)) {

      // Calculate photon luminosity
      if (cstar->ComputePhotonRates(TimeUnits, nbins, energies, Q) == FAIL) {
	ENZO_FAIL("Error in ComputePhotonRates.\n");
      }

      QTotal = 0;
      for (j = 0; j < nbins; j++) QTotal += Q[j];
      for (j = 0; j < nbins; j++) Q[j] /= QTotal;
      if (QTotal < tiny_number) continue;

#ifdef USE_MEAN_ENERGY
      double meanEnergy = 0;
      nbins = 1;
      for (j = 0; j < nbins; j++)
	meanEnergy += energies[j] * Q[j];
      meanEnergy /= QTotal;
      energies[0] = meanEnergy;
      Q[0] = QTotal;
#endif /* USE_MEAN_ENERGY */

      /* (TODO) If requested, calculate ramping time for the luminosity */

      float ramptime = 0.0;   // zero for no ramp
      float tdyn, ti;
      if (cstar->ReturnType() == PopII) {
	if (StarClusterUnresolvedModel) {  // Cen & Ostriker
	  ramptime = cstar->ReturnLifetime();
	} else {  // Wise & Cen
	  ramptime = TimeInYears * StarClusterMinDynamicalTime;
	}
      } else if (cstar->ReturnType() == PopIII)
	// should be an parameter or determined from the data
	ramptime = TimeInYears * 10e3;
      else if (cstar->ReturnType() == SimpleSource)
	ramptime = TimeInYears * 1e6 * SimpleRampTime;

      /* Transfer the shining particle properties to the radiative
	 transfer source particle */

      RadiationSourceEntry *RadSource;
      RadSource = cstar->RadiationSourceInitialize();
      RadSource->Luminosity     = QTotal * LConv;
      RadSource->RampTime       = ramptime;
      RadSource->EnergyBins     = nbins;
      RadSource->Energy         = new float[nbins];
      RadSource->SED            = new float[nbins];

      // AJE: new
      RadSource->LifeTime     = cstar->ReturnLifetime();  // in code (should be?)
      RadSource->CreationTime = cstar->ReturnBirthTime(); // in code

      // NOTE: RadSource->LWLuminosity and RadSource->FUVLuminosity
      //       are currently only used when applying optically thin
      //       (1/r^2) radiation profiles in Grid_AddH2Dissociation
      //       and Grid_AddPeHeating, AND only when running with
      //       RadiativeTransferOpticallyThinSourceClustering ON.
      //       (see CreateSourceClusteringTree). Otherwise,
      //       energy / photon bins are used for RT methods to create
      //       photon packages, and star particle lists are used
      //       (instead of RadSource) to compute luminosities in
      //       optically thin case without source clustering:

      /* AJE: 3/31/20 - I may not really need to zero these if they're ON */
      if (RadiativeTransferOpticallyThinH2){
	RadSource->LWLuminosity = Q[LWbin] * LConv;
        RadSource->IRLuminosity = Q[IRbin] * LConv;
      } else{
	RadSource->LWLuminosity = 0.0;
        RadSource->IRLuminosity = 0.0;
      }

      if (RadiativeTransferOpticallyThinFUV){
        RadSource->FUVLuminosity = Q[FUVbin] * LConv; // AJE: 3/31 make photon! * energies[7]; // NOTE: actual Luminosity (erg/s) - not photon luminosity
      } else{
        RadSource->FUVLuminosity = 0.0;
      }

      for (j = 0; j < nbins; j++) {
	RadSource->Energy[j] = energies[j];
	RadSource->SED[j]    = Q[j];
      }

      // if the source needs a beaming direction, define it here
      RadSource->Orientation    = NULL;

      if (GlobalRadiationSources->NextSource != NULL)
	GlobalRadiationSources->NextSource->PreviousSource = RadSource;
      GlobalRadiationSources->NextSource = RadSource;

    } // ENDIF is a radiation source?

  } // ENDFOR stars

  return SUCCESS;

}
