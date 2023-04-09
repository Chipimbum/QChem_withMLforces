#include <stdlib.h>
#include "qcio.h"
#include "qcio.h"
#include "qcio.h"
#include "qcmath.h"
#include "qchem.h"
#include "aimdman.h"
#include "convfac.h"
#include "Nukes.h"
#include "QString.h"
#include "MoleculeInput.h"
#include "rem_values.h"

/******************************************************************
 *                                                                *
 *  Member functions for class Nukes, which cotains nuclear       *
 *  cartesian coordinates and velocities.			  *
 *                                                                *
 *  JMH (3/05)                                                    *
 *                                                                *
 ******************************************************************/


//Nukes::Nukes(char *CartFile, char *VelocFile)  

//ERM COMMENTED THE FOLLOWING LINE by rps. Be careful!!)----------------------------------------
//Nukes::Nukes(const char *CartFile, const char *VelocFile, const char *ViewFile)   //rps modified
//!!!!!!!!!!!!!!!!!! ERM end of comment line-----------------------------------------------------

//ERM added this line to add ECart to the list that already rps had: -----------------------------
Nukes::Nukes(const char *CartFile, const char *VelocFile, const char *ViewFile, const char *ECartFile) //ERM modified
//end ERM added ECart to Nukes--------------------------------------------------------------------
{ 
   /* standard constructor -- inputs are file names for position, 
      velocity, and interatomic distance output */

   int NAtoms = rem_read(REM_NATOMS);
   int InitCoords = rem_read(REM_AIMD_INIT_COORDS);
   int IPrint = rem_read(REM_AIMD_PRINT);

   HasNHC = (rem_read(REM_AIMD_THERMOSTAT) == NOSE_HOOVER) ? true : false;
   if (HasNHC){
      AllocateNHC();
      //InitializeNHC();
   }else{
      NChain_NHC=0;
      Q=Xi=XiDot=NULL;
   }

   NucCarts   = QAllocDouble(3*NAtoms);
   NucVeloc   = QAllocDouble(3*NAtoms);  
   NucMass    = QAllocDouble(NAtoms);
   InvNucMass = QAllocDouble(3*NAtoms);
   
   double *Carts4 = NULL;
   double *Carts3 = NULL;
   get_carts(&Carts4,&Carts3,&AtNum,&NAtoms);

 
   // get atomic masses and convert from amu to a.u.
   defmass(&NAtoms,AtNum,NucMass);  
   double AMUtoAU = 1.0/ConvFac(ELECTRON_MASS_IN_AMU);
   VRscale(NucMass,NAtoms,AMUtoAU);
   for (int i = 0; i < NAtoms; i++) 
      for (int j = 0; j < 3; j++) InvNucMass[3*i+j] = 1.0/NucMass[i];
   if (IPrint > 3){
      printf("Nuclear masses in A.U. and A.M.U.\n");
      for (int i = 0; i < NAtoms; i++)
         printf("%4d %15.6e %12.6f\n",i+1,NucMass[i],NucMass[i]/AMUtoAU);
   }

   if (InitCoords <= 0) // the $molecule carts were just fine
      for (int i = 0; i < 3*NAtoms; i++) NucCarts[i] = Carts3[i];

   else if (InitCoords == RESTART){ // get carts from previous final coords
      printf("Reading nuclear coordinates from previous final geometry\n");
      FileMan_Open_Read(FILE_NUKES);
      FileMan(FM_READ,FILE_NUKES,FM_DP,3*NAtoms,0,FM_BEG,NucCarts);
      FileMan_Close(FILE_NUKES);
   }
   else QCrash("Illegal coordinate input option");



   // load velocities only if they exist on disk
   if (rem_read(REM_STATUS_VELOC) != 0) get_veloc(NucVeloc);


   File_NucCarts = QOpen(CartFile,"w");
   File_NucVeloc = QOpen(VelocFile,"w");
   File_NucView  = QOpen(ViewFile,"w");
 // ERM added----------------
   File_ECart = QOpen(ECartFile, "w");
 // end ERM added------------ 

   fprintf(File_NucCarts,"# Time/fs  Nuclear cartestian coordinates ");
   fprintf(File_NucVeloc,"# Time/fs  Nuclear cartestian velocities  ");

   fprintf(File_NucVeloc,"(atomic units)\n");

   if (rem_read(REM_INPUT_BOHR) > 0)
      fprintf(File_NucCarts,"(bohrs)\n");
   else
      fprintf(File_NucCarts,"(angstroms)\n");
  
   printf("Generating viewable xyz trajectory,\n");
   printf(" file (in Angstroms) at AIMD/View.xyz\n");

}



//Nukes::Nukes(Nukes &nukes, char *CartFile, char *VelocFile)
//ERM commented the following, in order to add ECart. BE CAREFUL!! --------------------------------
//Nukes::Nukes(Nukes &nukes, const char *CartFile, const char *VelocFile, const char *ViewFile) //rps
//!!!!!!!!!!!!!!!end ERM commented------------------------------------------------------------------

//ERM added----------------
Nukes::Nukes(Nukes &nukes, const char *CartFile, const char *VelocFile, const char *ViewFile, const char *ECartFile) //ERM 
//end ERM -----------------

{
   // copy constructor
   int NAtoms = rem_read(REM_NATOMS);
   this->NucCarts = QAllocDouble(3*NAtoms);  
   this->NucVeloc = QAllocDouble(3*NAtoms);  
   this->NucMass  = QAllocDouble(NAtoms);
 
   VRcopy(this->NucCarts,nukes.NucCarts,3*NAtoms);
   VRcopy(this->NucVeloc,nukes.NucVeloc,3*NAtoms);
   VRcopy(this->NucMass, nukes.NucMass,   NAtoms);
   
   File_NucCarts = QOpen(CartFile,"w");
   File_NucVeloc = QOpen(VelocFile,"w");
   File_NucView  = QOpen(ViewFile,"w");
 // ERM added----------------
   File_ECart = QOpen(ECartFile, "w");
 // end ERM added------------ 

   this->HasNHC = hasNHC();
   if (HasNHC){
      int N = NChain_NHC;
      this->NChain_NHC = N;
      this->Xi    = QAllocDouble(N);
      this->XiDot = QAllocDouble(N);
      this->Q     = QAllocDouble(N);
      VRcopy(this->Xi,   Xi,   N);
      VRcopy(this->XiDot,XiDot,N);
      VRcopy(this->Q,    Q,    N);
   }
}


void QFreeLoc(Nukes &xyz, const char *Module, int Line)
{   
  double *r = xyz.Ptr2Carts();
  double *v = xyz.Ptr2Veloc();
  double *M = xyz.Ptr2Mass();
  double *iM = xyz.Ptr2InvMass();

  // pass the info for which routine QFree'd the Nukes
  QFreeLoc(iM,Module,Line); 
  QFreeLoc(M,Module,Line); 
  QFreeLoc(v,Module,Line); 
  QFreeLoc(r,Module,Line);

  if (xyz.hasNHC()){
     double *ptr = xyz.Ptr2Xi_NHC();
     QFreeLoc(ptr,Module,Line);
     ptr = xyz.Ptr2XiDot_NHC();
     QFreeLoc(ptr,Module,Line);
     ptr = xyz.Ptr2Q_NHC();
     QFreeLoc(ptr,Module,Line);
  }
}




void Nukes::Print() 
{ 
  double Fac;
  int NAtoms = rem_read(REM_NATOMS);
  int AtomicUnits = rem_read(REM_INPUT_BOHR);

  
  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  if (AtomicUnits)
     printf("             Nuclear coordinates and velocities (a.u.)\n");
    
  else {
     printf("       Nuclear coordinates (Angst) and velocities (a.u.)\n");

     Fac = ConvFac(BOHRS_TO_ANGSTROMS);
     VRscale(NucCarts,3*NAtoms,Fac);

     // no unit conversion for the velocities -- they are always in a.u.
  
  }

  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  printf("             x        y        z            ");
  printf("v_x         v_y        v_z\n");
  printf("-----------------------------------");
  printf("-----------------------------------------\n");

   if(rem_read(REM_QM_MM_INTERFACE) > 0 && NAtoms > 10 && rem_read(REM_QMMM_PRINT) < 1){
      for(int i=0;i<5;i++){
        String AtSymb = AtomicSymbol(AtNum[i]);  
        char *atsymb = (char*)AtSymb;
        printf("%-4d %3s  %8.4f %8.4f %8.4f %14.4e %12.4e %12.4e\n",i+1,atsymb,
           NucCarts[3*i],NucCarts[3*i+1],NucCarts[3*i+2], 
           NucVeloc[3*i],NucVeloc[3*i+1],NucVeloc[3*i+2]);
      }
      printf("                                      .                             \n");
      printf("                                      .     (omitted for QM/MM)     \n");
      printf("                                      .                             \n");
      for(int i=NAtoms-5;i<NAtoms;i++){
        String AtSymb = AtomicSymbol(AtNum[i]);  
        char *atsymb = (char*)AtSymb;
        printf("%-4d %3s  %8.4f %8.4f %8.4f %14.4e %12.4e %12.4e\n",i+1,atsymb,
           NucCarts[3*i],NucCarts[3*i+1],NucCarts[3*i+2], 
           NucVeloc[3*i],NucVeloc[3*i+1],NucVeloc[3*i+2]);
      }
   }
   else{ 
     for (int i = 0; i < NAtoms; i++) {
        String AtSymb = AtomicSymbol(AtNum[i]);  
        char *atsymb = (char*)AtSymb;
        printf("%-4d %3s  %8.4f %8.4f %8.4f %14.4e %12.4e %12.4e\n",i+1,atsymb,
           NucCarts[3*i],NucCarts[3*i+1],NucCarts[3*i+2], 
           NucVeloc[3*i],NucVeloc[3*i+1],NucVeloc[3*i+2]);
     //ERM added -----
     //I want to print in a new file EandCart.xyz the energy and xyz at each time step for sGDML.
      /*  FILE *File_CartE;
        File_CartE = fopen("/home/scr/fanirm/", "ab+");
        
	fprintf(File_CartE, "%d", NAtoms); // I want the energy

        //fprintf(File_CartE,"%18.10e ",NucGrad[3*i+j]/NucMass[i]); //I want the energy
    
        fprintf(File_CartE,"%-4d %3s  %8.4f %8.4f %8.4f %14.4e %12.4e %12.4e\n",i+1,atsymb,
           NucCarts[3*i],NucCarts[3*i+1],NucCarts[3*i+2], 
           NucVeloc[3*i],NucVeloc[3*i+1],NucVeloc[3*i+2]);
    
*/
     //end ERM added ----------------------
     }  
  }
  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  if (!AtomicUnits){
     Fac = 1.0/Fac;
     VRscale(NucCarts,3*NAtoms,Fac); 
  }
}

//BRL FSSH begin
void Nukes::Print(FILE * pfile)
{
  double Fac;
  int NAtoms = rem_read(REM_NATOMS);
  int AtomicUnits = rem_read(REM_INPUT_BOHR);


  //printf("-----------------------------------");
  //printf("-----------------------------------------\n");

  if (AtomicUnits) {
    //printf("             Nuclear coordinates and velocities (a.u.)\n");
  }
  else {
    //printf("       Nuclear coordinates (Angst) and velocities (a.u.)\n");

     Fac = ConvFac(BOHRS_TO_ANGSTROMS);
     VRscale(NucCarts,3*NAtoms,Fac);

     // no unit conversion for the velocities -- they are always in a.u.

  }

  //printf("-----------------------------------");
  //printf("-----------------------------------------\n");

  //printf("             x        y        z            ");
  // printf("v_x         v_y        v_z\n");
  //printf("-----------------------------------");
  //printf("-----------------------------------------\n");

  fprintf(pfile,"%-4d",NAtoms);
  fprintf(pfile,"\n\n");

  for (int i = 0; i < NAtoms; i++) {
     String AtSymb = AtomicSymbol(AtNum[i]);
     char *atsymb = (char*)AtSymb;
     fprintf(pfile, "%1s  %17.15f %17.15f %17.15f %17.15e %17.15e %17.15e\n",atsymb,
        NucCarts[3*i],NucCarts[3*i+1],NucCarts[3*i+2],
            NucVeloc[3*i],NucVeloc[3*i+1],NucVeloc[3*i+2]);  //BRL more sigfigs for debugging
  }
  //printf("-----------------------------------");
  //printf("-----------------------------------------\n");

  if (!AtomicUnits){
     Fac = 1.0/Fac;
     VRscale(NucCarts,3*NAtoms,Fac);
  }
}


//BRL FSSH end

//Printing geometry in standard QChem format for WebMO and other software
void Nukes::PrintStandard() 
{ 
  double Fac= ConvFac(BOHRS_TO_ANGSTROMS);
  int NAtoms = rem_read(REM_NATOMS);
  int AtomicUnits = rem_read(REM_INPUT_BOHR);


  printf(" ----------------------------------------------------\n");
  if (AtomicUnits)
    printf("       Standard Nuclear Orientation (Bohr)\n");
  else
    {
      printf("       Standard Nuclear Orientation (Angstroms)\n");
      VRscale(NucCarts,3*NAtoms,Fac);
    }
  
  printf("    I     Atom         X            Y            Z\n");
  printf(" ----------------------------------------------------\n");
  for (int i=0; i < NAtoms; i++){
    String AtSymb = AtomicSymbol(AtNum[i]);
    char *atsymb = (char*)AtSymb;
    printf("%5d      %-2s %13.6f %12.6f %12.6f\n",i+1,atsymb,
	   NucCarts[3*i],NucCarts[3*i+1],NucCarts[3*i+2]);
  }
  printf(" ----------------------------------------------------\n");
  
  if (!AtomicUnits)
    {
      Fac = 1.0/Fac;
      VRscale(NucCarts,3*NAtoms,Fac); 
    }
}



void Nukes::PrintDiff(double *OldCarts, double *OldVeloc) 
{ 

  // print difference between two sets of coordinates 

  double Fac;
  int NAtoms = rem_read(REM_NATOMS);
  int AtomicUnits = rem_read(REM_INPUT_BOHR);

  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  if (AtomicUnits)
     printf("            Position and velocity displacements(a.u.)\n");
   
  else{
     printf("       Position (Angst) and velocities (a.u.) displacements\n");

     Fac = ConvFac(BOHRS_TO_ANGSTROMS);
     VRscale(NucCarts,3*NAtoms,Fac);
     VRscale(OldCarts,3*NAtoms,Fac);

     // no unit conversion for the velocities -- they are always in a.u.
  
  }

  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  printf("            dx       dy       dz            ");
  printf("dv_x        dv_y       dv_z\n");
  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  for (int i = 0; i < NAtoms; i++) {
     String AtSymb = AtomicSymbol(AtNum[i]);  
     char *atsymb = (char*)AtSymb;
     printf("%-4d %3s  %8.4f %8.4f %8.4f %14.4e %12.4e %12.4e\n",i+1,atsymb,
        NucCarts[3*i]-OldCarts[3*i],NucCarts[3*i+1]-OldCarts[3*i+1],
        NucCarts[3*i+2]-OldCarts[3*i+2], 
        NucVeloc[3*i]-OldVeloc[3*i],NucVeloc[3*i+1]-OldVeloc[3*i+1],
        NucVeloc[3*i+2]-OldVeloc[3*i+2]);
  }
  printf("-----------------------------------");
  printf("-----------------------------------------\n");

  if (!AtomicUnits){
     Fac = 1.0/Fac;
     //VRscale(NucVeloc,3*NAtoms,Fac);
     //Fac = 1.0/ConvFac(BOHRS_TO_ANGSTROMS);
     VRscale(NucCarts,3*NAtoms,Fac); 
     VRscale(OldCarts,3*NAtoms,Fac);
  }

}
 

void Nukes::Disk(double Time)
{ 
  int AtomicUnits = rem_read(REM_INPUT_BOHR); 
  int NAtoms  = rem_read(REM_NATOMS);
  int NSample = rem_read(REM_AIMD_NUCL_SAMPLE_RATE);

  static int NCalls = 0, Charge, Mult;
  static String Directory, FileRoot;


  double Fac;
  if (!AtomicUnits) {
    Fac = ConvFac(BOHRS_TO_ANGSTROMS);
    VRscale(NucCarts,3*NAtoms,Fac);
  }
      

  fprintf(File_NucCarts, "%.5f ", Time);
  fprintf(File_NucVeloc, "%.5f ", Time);

  // coordinates are written in whatever the input coords were
  // but velocities are always output in a.u.

  for (int i = 0; i < 3*NAtoms; i++)
     fprintf(File_NucCarts, "%16.6e ", NucCarts[i]); 
  fprintf(File_NucCarts,"\n");

  for (int i = 0; i < 3*NAtoms; i++)
     fprintf(File_NucVeloc, "%16.6e ", NucVeloc[i]); 
  fprintf(File_NucVeloc,"\n");

  //For viewable xyz file...

    //1st line is # of atoms
      fprintf(File_NucView,"%i\n",NAtoms);
    //2nd line is comment (here, is time in au)
      fprintf(File_NucView,"%8.3f\n",Time);
    //Rest is xyz format
      for (int i = 0; i < NAtoms; i++)
         fprintf(File_NucView, "%s %16.10f %16.10f %16.10f\n",
                 (char*)AtomicSymbol(AtNum[i]),
                 NucCarts[i*3+0],NucCarts[i*3+1],NucCarts[i*3+2]);

    //ERM copied the printing style used above for File_NucView that is the AIMD/View.xyz-------------------
    //For an xyz file that is used to estimate forces by sGDML

    //1st line is # of atoms
      fprintf(File_ECart,"%i\n",NAtoms);
    //2nd line is comment (here, SCF energy)
    double eSCF,eTot;
    //double Fac; // *** 
    //Fac = ConvFac(HARTREES_TO_EV); //*** this is the way they define in several files when they need to use it.
    FileMan(FM_READ,FILE_ENERGY,FM_DP,1,FILE_POS_CRNT_TOTAL_ENERGY,FM_BEG,&eTot);
    //VRscale(************) ; //*** example: VRscale(NucCarts,3*NAtoms,Fac);
    fprintf(File_ECart,"%.12f\n", eTot); //ERM: I need to later print this in eV because it is the units that sGDML uses
    
    //Rest is xyz format
      for (int i = 0; i < NAtoms; i++)
         fprintf(File_ECart, "%s %16.10f %16.10f %16.10f\n",
                 (char*)AtomicSymbol(AtNum[i]),
                 NucCarts[i*3+0],NucCarts[i*3+1],NucCarts[i*3+2]);
    //end ERM -------------------------------------------------------------------------------------------------

  if (!AtomicUnits){
     Fac = 1.0/Fac;
     VRscale(NucCarts,3*NAtoms,Fac); 
  }

  // these are too precious to lose in case of a crash
  fflush(File_NucCarts);
  fflush(File_NucVeloc);
  fflush(File_NucView);
    
  //ERM added this to have file with E and xyz from last frame alone
  fflush(File_ECart); //SUPER IMPORTANT! otherwise the file is lost and I read nothing! 
  system(". /home/scr/fanirm/interf_QC_SGDML/ECart_to_lastxyz.sh"); //work in python and see
  system(". /home/scr/fanirm/interf_QC_SGDML/force_query.sh"); 
  //end ERM added for printing of AIMD/ECart file and tailing it for sGMDL
  NCalls++;
}



void Nukes::Snapshot(double Time)
{ 
  int AtomicUnits = rem_read(REM_INPUT_BOHR); 
  int NAtoms = rem_read(REM_NATOMS);
  int NSample = rem_read(REM_AIMD_NUCL_SAMPLE_RATE);

  static int NCalls = 0, Charge, Mult;
  static String Directory, FileRoot;


  if (NCalls == 0){
     // open a directory for molecular snapshots

     Directory = QMkDir("AIMD/snapshots");
     FileRoot = Directory + String("/molecule.");

     MoleculeInput Mol(QFileName("molecule"));
     Charge = Mol.charge();
     Mult   = Mol.multiplicity();
  }


  double Fac;
  if (!AtomicUnits) {
    Fac = ConvFac(BOHRS_TO_ANGSTROMS);
    VRscale(NucCarts,3*NAtoms,Fac);
  }
      
 
  // write a separate file for each call 
 
  String FileName = FileRoot;
  FileName += NCalls*NSample;
  FILE *Snapshot = QOpen(FileName,"w");
   
  fprintf(Snapshot,"$comment\nMolecular geometry at t = %.5f\n$end\n\n",Time);
  fprintf(Snapshot,"$molecule\n%d %d\n",Charge,Mult);

  for (int i = 0; i < NAtoms; i++){
     String AtSymb = AtomicSymbol(AtNum[i]);  
     char *atsymb = (char*)AtSymb;
     fprintf(Snapshot,"%s  ",atsymb);
     for (int j = 0; j < 3; j++)
        fprintf(Snapshot,"%16.8f ",NucCarts[3*i+j]); 
     fprintf(Snapshot,"\n");
  } 
  fprintf(Snapshot,"$end\n");
  fclose(Snapshot);


  if (!AtomicUnits){  // convert back to a.u.
     Fac = 1.0/Fac;
     //VRscale(NucVeloc,3*NAtoms,Fac);
     //Fac = 1.0/ConvFac(BOHRS_TO_ANGSTROMS);
     VRscale(NucCarts,3*NAtoms,Fac); 
  }

  NCalls++;
}


  
double Nukes::KE()
{
   int NAtoms = rem_read(REM_NATOMS);
   double T = 0.0; 
   for (int i = 0; i < NAtoms; i++)
	   
      for (int j = 0; j < 3; j++)
         T += NucMass[i] * NucVeloc[3*i+j] * NucVeloc[3*i+j];
                 /*  //ERM added-----
      printf("print inside Nukes::KE, temperature %d \n",T);
      printf("print inside Nukes::KE, temperature %d \n",&T);
      printf("print inside Nukes::KE, temperature %p \n",T);
      printf("print inside Nukes::KE, temperature %p \n",&T);

   //printf("NucMass [i]= %.2f \n", &NucMass);
   //end of ERM added -----
*/
   return 0.5*T;
}



void Nukes::PrintNHC()
{
  printf("NHC thermost masses, coordinates, velocities\n");
  for (int i=0; i < nChain_NHC(); i++)
     printf("%4d %16.5e %16.5e %16.5e\n",i+1,Q[i],Xi[i],XiDot[i]);
}
