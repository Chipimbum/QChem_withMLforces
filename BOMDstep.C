#include "qcmath.h"
#include "qchem.h"
#include "aimdman.h"
#include "../setman/setman.h"
#include "BSetMgr.hh"
#include "convfac.h"
#include "rem_values.h"
#include "../drvman/drvman.h"
#include "Parallel.h"
#include <cassert>

#define ZETA1 0.1931833275037836
#define ZETA2 0.6136333449924328

//#include <vector> // -----JAGS------------
/******************************************************************
 *                                                                *
 * BOMDstep updates the nuclear coordinates and velocities per    *
 * on BOMD step.  					          *
 *                                                                *
 *  JMH (3/05)                                                    *
 *                                                                *
 ******************************************************************/

int BOMDstep(Nukes &nukes, double &NuclearKE)
{  
  static int NCalls = 0;


  int NAtoms = rem_read(REM_NATOMS);
  int IntMethod = rem_read(REM_AIMD_INTEGRATION);
  int IPrint = rem_read(REM_AIMD_PRINT);
  int NCycles;


  double dt = time_step_conversion * (double)rem_read(REM_TIME_STEP);
  double TimeAU = dt * (double)NCalls;
  double au2fs  = 1.0e15*ConvFac(AU_TIME_IN_SEC);
  double TimeFs = TimeAU*au2fs; 

  int NSampleNuc = rem_read(REM_AIMD_NUCL_SAMPLE_RATE);
  double NucSampleTime = dt * (double)NSampleNuc;


  double *NucCarts = nukes.Ptr2Carts();
  double *NucVeloc = nukes.Ptr2Veloc();
  double *NucMass  = nukes.Ptr2Mass();

  bool iscapaimd = false;  //JAGS

  static double *PrevNucGrad = QAllocDouble(3*NAtoms), 
                *NucGrad = QAllocDouble(3*NAtoms); 

  static FILE *File_NucForces = QOpen("AIMD/NucForces","w"); 

  if (NCalls == 0)
    fprintf(File_NucForces,"# Time/fs  Nuclear cartesian forces (a.u.)\n");
 
  if (rem_read(REM_CS_HF)>0) iscapaimd = true; //

  FileMan_Open_Read(FILE_NUCLEAR_GRADIENT);
  FileMan(FM_READ,FILE_NUCLEAR_GRADIENT,FM_DP,3*NAtoms,0,FM_BEG,NucGrad);
  FileMan_Close(FILE_NUCLEAR_GRADIENT); 
 
  //ERM added: I mimic what is done here above with the file NucForces, to create my New empty file for sGDML----
  //File number was added to a dictionary, its number is 11105
   static FILE *File_MLNucForces = QOpen("AIMD/MLNucForces","w"); //ERM Added

   if (NCalls == 0) //ERM Added
     fprintf(File_MLNucForces,"# Empty file created by Estefi during debugging \n"); //ERM Added
     fprintf(File_MLNucForces,"# some contain TBD at a later stage \n"); //ERM Added

  //end ERM -----------------------------------------------------------------------------------------------------

  /*   
  //ERM added: ----------------------------------------------------------------
  //I want a file with Energy and cartesian coordinates, pretty much like /AIMD/View.xyz but with E instead of tstep
  //File number was added to a dictionary, its number is 135
   static FILE *File_ECart = QOpen("AIMD/ECart","w");

   if (NCalls == 0) //ERM Added
     fprintf(File_ECart,"# Cartesian Coordinates and Energy for sGDML Force prediction.\n"); //ERM Added
  //end ERM -----------------------------------------------------------------------------------------------------
 */
  //ERM added: ----------------------------------------------------------------
  FileMan_Open_Read(FILE_ML_GRADIENT);
  //FileMan(FM_READ,FILE_NUCLEAR_GRADIENT,FM_DP,3*NAtoms,0,FM_BEG,NucGrad);
  FileMan_Close(FILE_ML_GRADIENT);
  //end ERM -----------------------------------------------------------------------------------------------------


  // ERM added 14March to control the printing to the NucForces file... ----------
  FILE *MLF_File;
    //MLF_File = fopen("/home/scr/fanirm/h20_dyns/FORCE_QUERY-2-good/integral/read_predF_from_file/6March/MLforces.dat", "r");
    MLF_File = fopen("/home/SCR/qcscr/test/AIMD/last_MLForces.dat", "r");
  double MLFor[9];
  for (int i = 0; i < 3*NAtoms; i++)
    {
        fscanf(MLF_File, "%lf,", &MLFor[i]);
        // ********** KEY CHANGE *************** 
        NucGrad[i] = MLFor[i]; //Replace the calculated NucGrad by the ML Force for the PRINTING IN THE SCRATCH FILE! 
        // *************************************
    }
 // end ERM added 14march ---------------------------------------------------------

   //ERM added: only a tmp printing while debugging ----------

   printf("-------------------------\n");
   printf("NucGrad [i] = MLFor[i] PRINTED to NucForces file by BOMDstep.C\n");
   printf("-------------------------\n");

   for (int i = 0; i < 3*NAtoms; i++) {
            printf(" %f \n",  NucGrad[i]);

        }
  
   for (int i = 0; i < NAtoms; i++) {
            printf(" %f \n",  NucMass[i]);

        }
   //end ERM added --------------------------------------------


  if (NCalls == 0) RemoveTransRot(NucGrad,nukes);

  int i;
  fprintf(File_NucForces,"%11.6f ",TimeFs);
  for (i = 0; i < NAtoms; i++)
    for (int j = 0; j < 3; j++)
      fprintf(File_NucForces,"%18.10e ",NucGrad[3*i+j]/NucMass[i]);
  fprintf(File_NucForces,"\n");

  if(rem_read(REM_AIMD_FRAG_PRINT) > 0) {
     rem_write(1, REM_STATUS_VELOC);
      set_veloc(NucVeloc,1);
  }

  if (IntMethod == OVV){
    double dtZeta = dt*ZETA1;
    for (i = 0; i < NAtoms; i++)
      for (int j = 0; j < 3; j++) NucCarts[3*i+j] += dtZeta*NucVeloc[3*i+j];

    set_carts(NucCarts,1); 
    CAPAimdFixer(NCalls,iscapaimd); //JAGS
    NCycles = BOMD_EnergyAndGrad(NucGrad);
    RemoveTransRot(NucGrad,nukes);

    dtZeta = dt*ZETA2;
    for (i = 0; i < NAtoms; i++)
      for (int j = 0; j < 3; j++) {
        int I = 3*i+j;
        NucVeloc[I] -= 0.5*dt*NucGrad[I]/NucMass[i];
        NucCarts[I] += dtZeta*NucVeloc[I];
      }

    set_carts(NucCarts,1); 
    if(rem_read(REM_AIMD_FRAG_PRINT) > 0)
    set_veloc(NucVeloc,1);
    CAPAimdFixer(NCalls,iscapaimd); // JAGS;
    NCycles += BOMD_EnergyAndGrad(NucGrad);
    RemoveTransRot(NucGrad,nukes);

    dtZeta = dt*ZETA1;
    for (i = 0; i < NAtoms; i++)
      for (int j = 0; j < 3; j++) {
        int I = 3*i+j;
        NucVeloc[I] -= 0.5*dt*NucGrad[I]/NucMass[i];
        NucCarts[I] += dtZeta*NucVeloc[I];
      }

    set_carts(NucCarts,1); 
    set_veloc(NucVeloc,1);
    if (IPrint >= 1) nukes.Print();

  }else{ 

    NuclearKE = NuclearUpdate(dt,nukes,NucGrad,PrevNucGrad);
    // ERM added  ------------
    printf ("NuclearKE in the BOMDstep.C step: %lf \n", NuclearKE); //ERM added
    // end ERM added --------------
    VRcopy(PrevNucGrad,NucGrad,3*NAtoms);
    set_carts(NucCarts,1);
    set_veloc(NucVeloc,1);
    CAPAimdFixer(NCalls,iscapaimd); // JAGS
    NCycles = BOMD_EnergyAndGrad(NucGrad);
    
    RemoveTransRot(NucGrad,nukes);
    if (IntMethod == VVERLET || IntMethod == BEEMAN || IntMethod == RATTLE) {
      NuclearKE = NuclearUpdatePart2(dt,nukes,NucGrad);
      set_veloc(NucVeloc,1);
    }
    if (IPrint >= 1) nukes.Print();
  }



  if (NSampleNuc > 0)
    if (NCalls%NSampleNuc == 0) {
      nukes.Disk((TimeAU+dt)*au2fs);
      //if (IPrint >= 3) nukes.Snapshot((TimeAU+dt)*au2fs);
    }

  // Velocity monitoring for quasiclassical calculations (DSL)
  if (rem_read(REM_AIMD_INIT_VELOC) == QUASICLASSICAL &&
      rem_read(REM_AIMD_STEPS) > 0) {
    MonitorKineticEnergy(nukes);
  }


  NCalls++;
  return NCycles;
}
