/***********************************************************************
/
/  GLOBAL DATA DECLARATIONS FOR THE STAR PARTICLES
/
/  written by: Greg Bryan
/  date:       February, 1997
/  modified1:
/
/  PURPOSE:
/
************************************************************************/
#ifndef STAR_PARTICLE_DATA_DEFINED__
#define STAR_PARTICLE_DATA_DEFINED__

#ifdef DEFINE_STORAGE
# define SPEXTERN
#else /* DEFINE_STORAGE */
# define SPEXTERN extern
#endif /* DEFINE_STORAGE */

#define STAR_PARTICLE_NUMBER_START 1000000000

/* #include "macros_and_parameters.h" */

struct ParticleEntry {
  FLOAT Position[3];
  float Mass;
  float Velocity[3];
  float Attribute[MAX_NUMBER_OF_PARTICLE_ATTRIBUTES];
  PINT Number;
  int Type;
};


/* Number of Star particles. */

SPEXTERN int NumberOfStarParticles;
SPEXTERN int NumberOfActiveParticles;
SPEXTERN int NumberOfDeletedParticles;
SPEXTERN PINT NumberOfOtherParticles; //all the particles other than type=2
SPEXTERN PINT NextActiveParticleID;
SPEXTERN int G_TotalNumberOfStars;

/* Star particle parameters. */

SPEXTERN int StarFeedbackType;
SPEXTERN int StarMakerTypeIaSNe;
SPEXTERN int StarMakerTypeIISNeMetalField;
SPEXTERN int StarMakerPlanetaryNebulae;
SPEXTERN int StarMakerTimeIndependentFormation;
SPEXTERN float StarMakerOverDensityThreshold;
SPEXTERN int StarMakerUseOverDensityThreshold;
SPEXTERN float StarMakerMaximumFractionCell;
SPEXTERN float StarMakerSHDensityThreshold;
SPEXTERN float StarMakerMassEfficiency;
SPEXTERN float StarMakerMinimumMass;
SPEXTERN float StarMakerMinimumDynamicalTime;
SPEXTERN float StarMassEjectionFraction;
SPEXTERN float StarMetalYield;
SPEXTERN float StarEnergyToThermalFeedback;
SPEXTERN float StarEnergyFeedbackRate;
SPEXTERN float StarEnergyToStellarUV;
SPEXTERN float StarEnergyToQuasarUV;
SPEXTERN int StarFeedbackDistRadius;
SPEXTERN int StarFeedbackDistCellStep;
SPEXTERN int StarFeedbackDistTotalCells;
SPEXTERN int StarFeedbackPreSN;
SPEXTERN float StarFeedbackKineticFraction;
SPEXTERN float StarMakerExplosionDelayTime;
SPEXTERN int StarMakerUseJeansMass;

SPEXTERN float PopIIIStarMass;
SPEXTERN int   PopIIIInitialMassFunction;
SPEXTERN int   PopIIIInitialMassFunctionSeed;
SPEXTERN int   PopIIIInitialMassFunctionCalls;
SPEXTERN float PopIIILowerMassCutoff;
SPEXTERN float PopIIIUpperMassCutoff;
SPEXTERN float PopIIIInitialMassFunctionSlope;
SPEXTERN int   PopIIIBlackHoles;
SPEXTERN float PopIIIBHLuminosityEfficiency;
SPEXTERN float PopIIIOverDensityThreshold;
SPEXTERN float PopIIIH2CriticalFraction;
SPEXTERN float PopIIIMetalCriticalFraction;
SPEXTERN int   PopIIIHeliumIonization;
SPEXTERN float PopIIISupernovaRadius;
SPEXTERN int   PopIIISupernovaUseColour;
SPEXTERN int   PopIIISupernovaMustRefine;
SPEXTERN int   PopIIISupernovaMustRefineResolution;
SPEXTERN float PopIIIColorDensityThreshold;
SPEXTERN float PopIIIColorMass;
SPEXTERN int   PopIIIUseHypernova;
SPEXTERN int   PopIIISupernovaExplosions;
SPEXTERN int   PopIIIOutputOnFeedback;
SPEXTERN int   PopIIIRadiationModel;

SPEXTERN int    StarClusterUseMetalField;
SPEXTERN int    StarClusterHeliumIonization;
SPEXTERN float  StarClusterMinDynamicalTime;
SPEXTERN double StarClusterIonizingLuminosity;
SPEXTERN double StarClusterSNEnergy;
SPEXTERN float  StarClusterSNRadius;
SPEXTERN float  StarClusterFormEfficiency;
SPEXTERN float  StarClusterMinimumMass;
SPEXTERN float  StarClusterCombineRadius;
SPEXTERN int    StarClusterUnresolvedModel;
SPEXTERN float  StarClusterRegionLeftEdge[3];
SPEXTERN float  StarClusterRegionRightEdge[3];

SPEXTERN float  MBHMinDynamicalTime;
SPEXTERN float  MBHMinimumMass;
SPEXTERN int    MBHAccretion;
SPEXTERN float  MBHAccretionRadius;
SPEXTERN float  MBHAccretingMassRatio;
SPEXTERN float  MBHAccretionFixedTemperature;
SPEXTERN float  MBHAccretionFixedRate;
SPEXTERN int    MBHTurnOffStarFormation;
SPEXTERN float  MBHCombineRadius;

SPEXTERN float UnfulfilledStarFormationMass;

SPEXTERN int    MBHFeedback;
SPEXTERN float  MBHFeedbackRadiativeEfficiency;
SPEXTERN float  MBHFeedbackEnergyCoupling;
SPEXTERN float  MBHFeedbackMassEjectionFraction;
SPEXTERN float  MBHFeedbackMetalYield;
SPEXTERN float  MBHFeedbackThermalRadius;
SPEXTERN float  MBHFeedbackJetsThresholdMass;

SPEXTERN float  H2StarMakerEfficiency;
SPEXTERN float  H2StarMakerNumberDensityThreshold;
SPEXTERN float  H2StarMakerMinimumMass;
SPEXTERN float  H2StarMakerMinimumH2FractionForStarFormation;
SPEXTERN int    H2StarMakerStochastic;
SPEXTERN int    H2StarMakerUseSobolevColumn;
SPEXTERN float  H2StarMakerSigmaOverR;
SPEXTERN int    H2StarMakerAssumeColdWarmPressureBalance;
SPEXTERN float  H2StarMakerH2DissociationFlux_MW;
SPEXTERN float  H2StarMakerH2FloorInColdGas;
SPEXTERN float  H2StarMakerColdGasTemperature;

SPEXTERN int   IndividualStarRefineToLevel;
SPEXTERN float IndividualStarRefineToPhysicalRadius;
SPEXTERN int   IndividualStarRefineForRadiation;
SPEXTERN int   IndividualStarRefineBufferSize;
SPEXTERN float IndividualStarTemperatureLimit;
SPEXTERN float IndividualStarTemperatureLimitFactor;
SPEXTERN float IndividualStarICSupernovaRate;
SPEXTERN float IndividualStarICSupernovaTime;
SPEXTERN float IndividualStarICSupernovaR;
SPEXTERN float IndividualStarICSupernovaZ;
SPEXTERN int   IndividualStarICSupernovaMethod;
SPEXTERN int   IndividualStarICSupernovaInjectionMethod;
SPEXTERN int   IndividualStarICSupernovaFromFile;
SPEXTERN float IndividualStarWDFixedLifetime;
SPEXTERN float IndividualStarRefineTime;
SPEXTERN int   ICSupernovaNumberOfPoints;
SPEXTERN float ICSupernovaTimeArray[500];
SPEXTERN float ICSupernovaSNRArray[500];
SPEXTERN float IndividualStarICSupernovaPos[MAX_DIMENSION];
SPEXTERN int   IndividualStarCheckVelocityDiv;
SPEXTERN int   IndividualStarICLifetimeMode;
SPEXTERN int   IndividualStarFeedbackOverlapSample;
SPEXTERN int   IndividualStarUseWindMixingModel;
SPEXTERN int   IndividualStarExtrapolateYields;
SPEXTERN float IndividualStarSFGasMassThreshold;
SPEXTERN float IndividualStarTemperatureThreshold;
SPEXTERN float IndividualStarSecondaryOverDensityThreshold;
SPEXTERN float IndividualStarIMFUpperMassCutoff;
SPEXTERN float IndividualStarIMFLowerMassCutoff;
SPEXTERN float IndividualStarIMFMassFloor;
SPEXTERN float IndividualStarVelocityDispersion;
SPEXTERN int   IndividualStarIMFSeed;
SPEXTERN int   IndividualStarIMF;
SPEXTERN int   IndividualStarIMFCalls;
SPEXTERN int   IndividualStarIgnoreNegativeMass;
SPEXTERN float IndividualStarSalpeterSlope;
SPEXTERN float IndividualStarKroupaAlpha1;
SPEXTERN float IndividualStarKroupaAlpha2;
SPEXTERN float IndividualStarKroupaAlpha3;
SPEXTERN float IndividualStarMassFraction;
SPEXTERN float IndividualStarSNIIMassCutoff;
SPEXTERN float IndividualStarAGBThreshold;
SPEXTERN int   IndividualStarTrackAGBMetalDensity;
SPEXTERN int   IndividualStarTrackSNMetalDensity;
SPEXTERN int   IndividualStarTrackWindDensity;
SPEXTERN int   IndividualStarRProcessModel;
SPEXTERN float IndividualStarRProcessMinMass;
SPEXTERN float IndividualStarRProcessMaxMass;
SPEXTERN float IndividualStarAGBWindVelocity;
SPEXTERN float IndividualStarDirectCollapseThreshold;
SPEXTERN float IndividualStarRadiationMinimumMass;
SPEXTERN int   IndividualStarStellarWinds;
SPEXTERN float IndividualStarWindTemperature;
SPEXTERN int   IndividualStarFollowStellarYields;
SPEXTERN int   IndividualStarSurfaceAbundances;
SPEXTERN int   LimongiAbundances;
SPEXTERN int   IndividualStarOutputChemicalTags;
SPEXTERN char* IndividualStarChemicalTagFilename;
SPEXTERN int   IndividualStarSaveTablePositions;
SPEXTERN float IndividualStarWDMinimumMass;
SPEXTERN float IndividualStarWDMaximumMass;
SPEXTERN float IndividualStarSNIaMinimumMass;
SPEXTERN float IndividualStarSNIaMaximumMass;
SPEXTERN float IndividualStarDTDSlope;
SPEXTERN float IndividualStarSNIaFraction;
SPEXTERN int   IndividualStarSNIaModel;
SPEXTERN int   IndividualStarBlackBodyOnly;
SPEXTERN float IndividualStarBlackBodyq0Factors[2];
SPEXTERN float IndividualStarBlackBodyq1Factors[2];
SPEXTERN float IndividualStarBlackBodyq2Factors[2];
SPEXTERN float IndividualStarBlackBodyFUVFactors[2];
SPEXTERN float IndividualStarBlackBodyLWFactors[2];
SPEXTERN float IndividualStarBlackBodyIRFactors[2];
SPEXTERN int   IndividualStarCreationStencilSize;
SPEXTERN float IndividualStarFeedbackStencilSize;
SPEXTERN float IndividualStarFeedbackRadius;
SPEXTERN float IndividualStarSupernovaEnergy;
SPEXTERN float IndividualStarStellarWindVelocity;
SPEXTERN float IndividualStarMaximumStellarWindVelocity;
SPEXTERN int   IndividualStarFUVHeating;
SPEXTERN int   IndividualStarLWRadiation;
SPEXTERN int   IndividualStarIRRadiation;
SPEXTERN int   IndividualStarPrintSNStats;

SPEXTERN float IndividualStarIonizingRadiationMinimumMass;
SPEXTERN float IndividualStarOTRadiationMass;
SPEXTERN float IndividualStarFUVTemperatureCutoff;

SPEXTERN int   IndividualStarPopIIIFormation;
SPEXTERN int   IndividualStarPopIIISeparateYields;

SPEXTERN int AccretingParticleRadiation;
SPEXTERN double AccretingParticleLuminosity;

SPEXTERN float minStarLifetime;
SPEXTERN FLOAT LastSupernovaTime;

SPEXTERN float *IMFData;
SPEXTERN float *SecondaryIMFData;
SPEXTERN float *EventDTD;

SPEXTERN float TypeIILowerMass, TypeIIUpperMass;
SPEXTERN float PISNLowerMass, PISNUpperMass;
SPEXTERN float PopIIIPISNEnergy;

/* for star particle minimum mass ramp */
SPEXTERN int StarMakerMinimumMassRamp;
SPEXTERN float StarMakerMinimumMassRampStartTime;
SPEXTERN float StarMakerMinimumMassRampStartMass;
SPEXTERN float StarMakerMinimumMassRampEndTime;
SPEXTERN float StarMakerMinimumMassRampEndMass;

#endif
