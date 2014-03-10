/* GO-MODEL BUILDER V1.1 */
// This program makes prm and rtf files for a single bead GO-liked model.
// For the input only XYZ coordinates of the C_alpha of single chain protein is required.
// The assumption is that there are no missing amino acids in the middle. 
// "N" the number of amino acids and "rcut" need to be modified.  

#include <iostream>
#include <fstream>
#include <math.h>
#include <cstring>
#include <cstdlib>   // for exit function
using namespace std;

void read_cor (double a[][4], const char *name ); // read coordinates
void read_emin (const char *name );  // read emin values
void check_cor (double a[][4], const char *name ); 
void check_emin(const char *name);
void Mktop (const char *name);
void Mkprm (const char *name);

double rcut=8.;          // cut-off distance for GO model
const int n=93;          // number of amino acids

double ca[n][4]={0};
double cb[n][4]={0};
double emin[21][21]={0};

int main ()
{

  read_cor (ca, "ch1_ca.xyz");  // ch1_ca;alpha carbones; format: x,y,z,id ; id is defined based on the emin tables
  read_cor (cb, "ch1_cb.xyz");  // ch1_cb;center of mass of each residue; format: x,y,z, id 
  read_emin ("emin_b.dat");  // emin.dat (kgs table; 21*21) ; emin_b.dat (BT table, 20*20); 
  // WARNING: default value for reading emin file to an array=20*20!

  check_cor (ca, "check_ca.txt");
  check_cor (cb, "check_cb.txt");
  check_emin ("check_emin.txt");
 
  Mktop ("ch1_ca_cb.rtf");  
  Mkprm ("ch1_ca_cb.prm"); 
  return 0;
}

/* Functions start from here */
/*---------------------------------------------------------------*/
void Mkprm (const char *name)
{
  double d_aa;
  double d_ab;
  double d_bb;

  double r_ab;
  double s;
  double e_min;
  ofstream prmfile;
  prmfile.open (name);

  cout << "Making CHARMM parameter file for single bead coarse-grained GO-like model with the cutoff of " << rcut << " A ... " << endl; 

// BOND
  prmfile << "BOND	! bonded potential" <<endl ;
  prmfile << "! K(CA-CA)=100*epsilon/a^2" <<endl ;
  prmfile << "! atom  atom  force-ct  r_min" <<endl <<endl ;

  for( int i=1; i<n; i++){
      prmfile << "A" <<i << "   A" << i+1 <<"   100.0    3.8" <<endl ;

      r_ab= pow((ca[i-1][0]-cb[i-1][0]),2.0) + pow((ca[i-1][1]-cb[i-1][1]),2.0) + pow((ca[i-1][2]-cb[i-1][2]),2.0);

      prmfile << "A" <<i << "   B" << i <<"   200.0    " << sqrt(r_ab) <<endl ;
  }
 
// ANGLE
  prmfile << endl << "ANGLe	! bond-angle  potential" <<endl ;
  prmfile << "! atom  atom  atom  kappa   theta0" <<endl <<endl ;
  for( int i=1; i<n-1; i++){
      prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   12.5     105.0" <<endl ;
      prmfile << "A" <<i << "   A" << i+1  << "   B"<<i+1 <<"   12.5     105.0" <<endl ;
      prmfile << "B" <<i << "   A" << i    << "   A"<<i+1 <<"   12.5     105.0" <<endl ;
  }

// DIHEDRAL
  prmfile << endl << "PHI	! dihedral angle  potential" <<endl ;
  prmfile << "! dihedral_pot = kappa * [1+col(n*phi-delta]" << endl ;
  prmfile << "! atom  atom  atom atom  kappa  n  delta" <<endl <<endl ;
  for( int i=1; i<n-2; i++){
//      prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   A" << i+3 << "  0.90 3  -66.5"  <<endl ;
//      prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   A" << i+3 << "  2.27 2  -68.1"  <<endl ;
//      prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   A" << i+3 << "  2.91 1  -37.3"  <<endl ;

   prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   A" << i+3 << "  1.5  1  0.0" <<endl ;
   prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   A" << i+3 << "  1.5  3  0.0"   <<endl ;

   prmfile << "B" <<i << "   A" << i    << "   A"<<i+1 <<"   A" << i+2 << "  1.5  1  0.0" <<endl ;
   prmfile << "B" <<i << "   A" << i    << "   A"<<i+1 <<"   A" << i+2 << "  1.5  3  0.0"   <<endl ;

   prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   B" << i+2 << "  1.5  1  0.0" <<endl ;
   prmfile << "A" <<i << "   A" << i+1  << "   A"<<i+2 <<"   B" << i+2 << "  1.5  3  0.0"   <<endl ;

   prmfile << "B" <<i << "   A" << i    << "   A"<<i+1 <<"   B" << i+1 << "  1.5  1  0.0" <<endl ;
   prmfile << "B" <<i << "   A" << i    << "   A"<<i+1 <<"   B" << i+1 << "  1.5  3  0.0"   <<endl ;

  }
  
//NONBOND
  prmfile << endl << "NBONDed CUTNB 16.0 CTOFNB 14.0 CTONNB 12.0 WMIN 1.5 E14Fac 1.0 EPS 1.0 -" << endl ;
  prmfile << "ATOM CDIE SWITCH VSWITCH VATOM BYGROUP ! non-bonded potential" <<endl ;
  prmfile << "! atom           e_min    r_min/2" << endl;
  for( int i=1; i<=n; i++){
      prmfile <<"A" << i << "   0.0    -1e-12     22.71" << endl;
      prmfile <<"B" << i << "   0.0    -1e-12     22.71" << endl; 

 }

//NBFIX 
  prmfile << endl << "NBFIX" <<endl ;
  prmfile << "! modify certain interactions, here rmin is given instead of rmin/2" <<endl;
  prmfile << "! atom atom    e_min   r_min   s_ij"<<endl <<endl ;
  for(int i=1; i<n-2;i++)
     for(int j=i+3; j<=n; j++) {

     d_aa = ( (ca[j-1][0]-ca[i-1][0]),2.0) +
            ( (ca[j-1][1]-ca[i-1][1]),2.0) + 
            ( (ca[j-1][2]-ca[i-1][2]),2.0);

     d_bb = ( (cb[j-1][0]-cb[i-1][0]),2.0) +
            ( (cb[j-1][1]-cb[i-1][1]),2.0) + 
            ( (cb[j-1][2]-cb[i-1][2]),2.0);
 
     d_ab = ( (cb[j-1][0]-ca[i-1][0]),2.0) +
            ( (cb[j-1][1]-ca[i-1][1]),2.0) + 
            ( (cb[j-1][2]-ca[i-1][2]),2.0);
 

	 if ( i>=16 or i<=10 ) {
	
	s = (emin[int(ca[i-1][3])][int(ca[j-1][3])]);

	e_min = (double(1.+0.6*fabs(emin[int(ca[i-1][3])][int(ca[j-1][3])]))*(-2.));


	if ( sqrt(d_aa) <= rcut ) {
		if (s > 0)
            		prmfile << "A" << i << "  A" << j << "  " << e_min  << "  " << sqrt(d_aa) << "  " << "1"  <<endl; 
          	else if (s <0) 
          	  	prmfile << "A" << i << "  A" << j << "  " << e_min * (-2./3.)  << "  " << sqrt(d_aa) << "  " << "-1"  <<endl;
          	else
            		prmfile << "A" << i << "  A" << j << "  " << e_min  << "  " << sqrt(d_aa) << "  " << "0"  <<endl;
			 }

        if ( sqrt(d_ab) <= rcut ) {
                if (s > 0)
                        prmfile << "A" << i << "  B" << j << "  " << e_min  << "  " << sqrt(d_ab) << "  " << "1"  <<endl; 
                else if (s <0) 
                        prmfile << "A" << i << "  B" << j << "  " << e_min * (-2./3.)  << "  " << sqrt(d_ab) << "  " << "-1"  <<endl;
                else
                        prmfile << "A" << i << "  B" << j << "  " << e_min  << "  " << sqrt(d_ab) << "  " << "0"  <<endl;
                         }


        if ( sqrt(d_bb) <= rcut ) {
                if (s > 0)
                        prmfile << "B" << i << "  B" << j << "  " << e_min  << "  " << sqrt(d_bb) << "  " << "1"  <<endl; 
                else if (s <0) 
                        prmfile << "B" << i << "  B" << j << "  " << e_min * (-2./3.)  << "  " << sqrt(d_bb) << "  " << "-1"  <<endl;
                else
                        prmfile << "B" << i << "  B" << j << "  " << e_min  << "  " << sqrt(d_bb) << "  " << "0"  <<endl;
                         }

    		}
 }
  prmfile.close();
// [int(a[j][4])] 
// (1+0.6*(a2[int(a[i-1][3])][int(a[j-1][3])]))*(-2)


 cout << "Parameter file wriiten into "<< name << endl;

}
/*-----------------------------------------------------------------*/
void Mktop(const char *name)
{

  ofstream topfile;
  topfile.open (name);
  
  cout << "Making CHARMM topology file for single bead coarse-grained GO-like model... " << endl;

  topfile << "* Topology file for a single bead GO model system" <<endl;
  topfile << "*" << endl << "    28    6" << endl << endl;  

  for( int i=1; i<=n; i++) {
      topfile << "MASS     " << i << "   A" << i << "   100" ;  
      topfile << endl ;  }

  for( int i=n+1; i<=2*n; i++) {
      topfile << "MASS     " << i << "   B" << i-n << "   200" ;  
      topfile << endl ;  }

  topfile << endl ;

  topfile << "DECL -CA    ! '-'=preceding atom"  << endl ; 
  topfile << "DECL +CA  ! '+'=next atom"  << endl ; 
  topfile << "DECL #CA  ! '#'=second next atom"  << endl<<endl ; 

  topfile << endl << "AUTOGENERATE ANGLES DIHEDRAL"  << endl <<endl ; 


  for( int i=1; i<=n; i++) {
      topfile << "RESI P" << i << " 0.00" << endl ;  
      topfile << "GROUP" << endl ; 
      topfile << "ATOM CA A" << i << " 0.00" <<  endl; 
      topfile << "ATOM CB B" << i << " 0.00" <<  endl; 
      if (i<n)  topfile << "BOND CA +CA   CA CB" << endl ; 
      if (i==n)  topfile << "BOND CA CB" << endl ; 
      topfile << endl;  }

  topfile.close();

  cout << "Topology file written into " << name << endl; 

}

/*----------------------------------------------------------------*/
void read_cor(double a[][4], const char *name) {

 ifstream in(name);

  if(!in) {
     cout << "Error: file could not be opened" << endl;
     return;
  }

  // filling the arrays
  for( int i=0; i<n; i++)
    for( int j=0; j<4; j++) {
         in >> a[i][j]; 
       }
  in.close();

}

/*-----------------------------------------------------------*/
void read_emin(const char *name) {

 ifstream in(name);

if(!in) {
 cout << "Error: file could not be opened";
 return;
 }

for( int i=0; i<20; i++) {
  for( int j=0; j<20; j++) {
 in >> emin[i][j];
// emin[i][j]= emin[i][j];
     }
   }
 in.close();

 }

 /*------------------------------------------------------------*/
void check_cor (double a[][4], const char *name) 
 {
// checking the read data by writing into a file
 ofstream myfile;
 myfile.open (name);

 for( int i=0; i<n; i++) {
  for( int j=0; j<4; j++) {
    myfile << a[i][j] << " "  ;  } 
    myfile << endl ;  }

  myfile.close();

}
 /*---------------------------------------------------*/
void check_emin (const char *name) 
 {
// checking the read data by writing into a file
 ofstream myfile;
 myfile.open (name);

 for( int i=0; i<20; i++) {
  for( int j=0; j<20; j++) {
    myfile << emin[i][j] << " "  ;  } 
    myfile << endl ;  }

  myfile.close();

}

