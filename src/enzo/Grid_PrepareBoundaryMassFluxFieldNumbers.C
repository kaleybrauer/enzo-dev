/***********************************************************************
/
/  GRID CLASS (INITIALIZE CONTAINERS FOR COUNTING DOMAIN BOUNDARY MASS FLUX)
/
/  written by: Andrew Emerick
/  date:       2018-ish
/
/  PURPOSE:
/           Sets up the container that keeps track of how much mass
/           in each of the mass-density fields move across the external
/           boundary
/
/  RETURNS: FAIL or SUCCESS
/
************************************************************************/

#include <stdio.h>
#include "ErrorExceptions.h"
#include "macros_and_parameters.h"
#include "typedefs.h"
#include "global_data.h"
#include "Fluxes.h"
#include "GridList.h"
#include "ExternalBoundary.h"
#include "Grid.h"

int FindField(int f, int farray[], int n);


int grid::PrepareBoundaryMassFluxFieldNumbers(){

  if (ProcessorNumber != MyProcessorNumber)
    return SUCCESS;

  if (!StoreDomainBoundaryMassFlux) return SUCCESS;

  if (BoundaryMassFluxFieldNumbers[0] >= 0) return SUCCESS; // we've already initialized

  /* obtain baryon field indexes */
  int DensNum, GENum, TENum, Vel1Num, Vel2Num, Vel3Num, CRNum, B1Num, B2Num, B3Num;
  this->DebugCheck("StarParticleHandler");
  if (this->IdentifyPhysicalQuantities(DensNum, GENum, Vel1Num, Vel2Num,
                                       Vel3Num, TENum, B1Num, B2Num, B3Num) == FAIL) {
    ENZO_FAIL("Error in IdentifyPhysicalQuantities.");
  }

  /* identify species fields if they exist */
  int DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum, HMNum, H2INum, H2IINum,
      DINum, DIINum, HDINum;
  if ( MultiSpecies ){
    IdentifySpeciesFields(DeNum, HINum, HIINum, HeINum, HeIINum, HeIIINum,
                          HMNum, H2INum, H2IINum, DINum, DIINum, HDINum);
  }

  /* identify colour fields if they exist */
  int SNColourNum, MetalNum, MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum,
    MetalIaNum, MetalIINum;

  if (this->IdentifyColourFields(SNColourNum, MetalNum, MetalIaNum, MetalIINum,
              MBHColourNum, Galaxy1ColourNum, Galaxy2ColourNum) == FAIL)
    ENZO_FAIL("Error in grid->IdentifyColourFields.\n");

  /* pre-compute field numbers so we only need to do this once */
  int count = 0;
  BoundaryMassFluxFieldNumbers[count++] = DensNum;

  /* H, He and ionization states */
  if (MultiSpecies > 0){
    BoundaryMassFluxFieldNumbers[count++] = HINum;
    BoundaryMassFluxFieldNumbers[count++] = HIINum;
    BoundaryMassFluxFieldNumbers[count++] = HeINum;
    BoundaryMassFluxFieldNumbers[count++] = HeIINum;
    BoundaryMassFluxFieldNumbers[count++] = HeIIINum;
  }

  /* H-, H2, H2- */
  if (MultiSpecies > 1){
    BoundaryMassFluxFieldNumbers[count++] = HMNum;
    BoundaryMassFluxFieldNumbers[count++] = H2INum;
    BoundaryMassFluxFieldNumbers[count++] = H2IINum;
  }

  if (MultiSpecies > 2){
    BoundaryMassFluxFieldNumbers[count++] = DINum;
    BoundaryMassFluxFieldNumbers[count++] = DIINum;
    BoundaryMassFluxFieldNumbers[count++] = HDINum;
  }

  if (MetalNum != -1){
    BoundaryMassFluxFieldNumbers[count++] = MetalNum;     // metallicity
    if (MultiMetals || TestProblemData.MultiMetals) {
      BoundaryMassFluxFieldNumbers[count++] = MetalNum+1; // ExtraType0
      BoundaryMassFluxFieldNumbers[count++] = MetalNum+2; // ExtraType1
    }
  }

  if (MetalIaNum       != -1) BoundaryMassFluxFieldNumbers[count++] = MetalIaNum;
  if (MetalIINum       != -1) BoundaryMassFluxFieldNumbers[count++] = MetalIINum;
  if (SNColourNum      != -1) BoundaryMassFluxFieldNumbers[count++] = SNColourNum;
  if (MBHColourNum     != -1) BoundaryMassFluxFieldNumbers[count++] = MBHColourNum;
  if (Galaxy1ColourNum != -1) BoundaryMassFluxFieldNumbers[count++] = Galaxy1ColourNum;
  if (Galaxy2ColourNum != -1) BoundaryMassFluxFieldNumbers[count++] = Galaxy2ColourNum;

  /* stellar yields tracers */
  if (MultiMetals > 1 && STARMAKE_METHOD(INDIVIDUAL_STAR)){
    for (int yield = 0; yield < StellarYieldsNumberOfSpecies; yield++){
      int anum = StellarYieldsAtomicNumbers[yield];
      if (anum == 1 || anum == 2) continue; // handled above

      int field_num;
      this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, anum);
      BoundaryMassFluxFieldNumbers[count++] = field_num;
    }

    int PopIIIMetalNum, PopIIIPISNeMetalNum, AGBMetalNum, RProcMetalNum,
        WindMetalNum, WindMetalNum2; //, SNIaMetalNum, SNIIMetalNum;
    int ExtraMetalNum0, ExtraMetalNum1, ExtraMetalNum2;

    AGBMetalNum         = FindField(ExtraType0, FieldType, NumberOfBaryonFields);
    PopIIIMetalNum      = FindField(ExtraType1, FieldType, NumberOfBaryonFields);
    PopIIIPISNeMetalNum = FindField(MetalPISNeDensity, FieldType, NumberOfBaryonFields);
    WindMetalNum        = FindField(MetalWindDensity, FieldType, NumberOfBaryonFields);
    WindMetalNum2       = FindField(MetalWindDensity2, FieldType, NumberOfBaryonFields);
    RProcMetalNum       = FindField(MetalRProcessDensity, FieldType, NumberOfBaryonFields);

    ExtraMetalNum0 = FindField(ExtraMetalField0, FieldType, NumberOfBaryonFields);
    ExtraMetalNum1 = FindField(ExtraMetalField1, FieldType, NumberOfBaryonFields);
    ExtraMetalNum2 = FindField(ExtraMetalField2, FieldType, NumberOfBaryonFields);

    if (IndividualStarTrackAGBMetalDensity){
      BoundaryMassFluxFieldNumbers[count++] = AGBMetalNum;
    }

    if (IndividualStarPopIIIFormation){
      BoundaryMassFluxFieldNumbers[count++] = PopIIIMetalNum;
      BoundaryMassFluxFieldNumbers[count++] = PopIIIPISNeMetalNum;

      if (IndividualStarPopIIISeparateYields){
        for (int yield = 0; yield < StellarYieldsNumberOfSpecies; yield++){
          int anum = StellarYieldsAtomicNumbers[yield];
          if (anum == 1 || anum == 2) continue; // ignore entirely here

          int field_num;
          this->IdentifyChemicalTracerSpeciesFieldsByNumber(field_num, anum, 0, 2);
          BoundaryMassFluxFieldNumbers[count++] = field_num;
        }
      }
    }



    if (IndividualStarTrackWindDensity){
      BoundaryMassFluxFieldNumbers[count++] = WindMetalNum;
      BoundaryMassFluxFieldNumbers[count++] = WindMetalNum2;
    }

    if (IndividualStarTrackSNMetalDensity && IndividualStarSNIaModel == 2){
      BoundaryMassFluxFieldNumbers[count++] = ExtraMetalNum0;
      BoundaryMassFluxFieldNumbers[count++] = ExtraMetalNum1;
      BoundaryMassFluxFieldNumbers[count++] = ExtraMetalNum2;
    }

    if (IndividualStarRProcessModel){
      BoundaryMassFluxFieldNumbers[count++] = RProcMetalNum;
    }
  }

  return SUCCESS;
}
