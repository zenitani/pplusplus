// ------------------------------------------------------ -*- C++ -*-
//     RK.h         by S.Zenitani     last updated : 1999/05/24 
// ------------------------------------------------------------------
//   invariant matrix for 4/6th order Runge-Kutta Method


#ifndef _Z_RK_H_
#define _Z_RK_H_

#include <vector3.h>


const double st44[4][4] = {
  { 0.5, 0.0, 0.0,  0.5 },
  { 0.0, 0.5, 0.0,  0.5 },
  { 0.0, 0.0, 1.0,  1.0 },
  { 0.166666666666667,
    0.333333333333333,
    0.333333333333333,
    0.166666666666667
  }
};

const double st76[7][7] = {
  {
    0.33333333333333333333, 
    0.0, 
    0.0,
    0.0,
    0.0,
    0.0,
    0.33333333333333333333
  },
  {
    0.0,
    0.66666666666666666667, 
    0.0,
    0.0,
    0.0,
    0.0,
    0.66666666666666666667
  },
  {
    0.083333333333333333333, 
    0.33333333333333333333, 
    -0.083333333333333333333, 
    0.0, 
    0.0, 
    0.0,
    0.333333333333333333333
  },
  {
    -0.0625, 
    1.125, 
    -0.1875, 
    -0.375, 
    0.0, 
    0.0,
    0.5
  },
  {
    0.0, 
    1.125, 
    -0.375, 
    -0.75, 
    0.5, 
    0.0,
    0.5
  },
  {
    0.204545454545454545455, 
    -0.81818181818181818182, 
    1.43181818181818181818, 
    1.63636363636363636364, 
    0.0, 
    -1.45454545454545454545,
    1.0
  },
  {
    0.091666666666666666667,
    0.0,                    
    0.675,                  
    0.675,                  
    -0.266666666666666666667,
    -0.266666666666666666667,
    0.091666666666666666667
  }
};


vector3 f( const vector3&, const double& );
double  f( const double& , const double& );


void RK4( vector3& y, vector3 (*f)( const vector3&, const double& ),
	  double& x, const double& h )
{
  static int i,j;
  static vector3 k[4];
  static vector3 tmp;

  // k1
  k[0] = f( y, x );
  
  // k2 ... k4
  for( i=0; i<3; i++ ){
    tmp = st44[i][0] * k[0];
    for( j=1; j<(i+1); j++ ){
      tmp += st44[i][j] * k[j];
    }
    k[i+1] = f( y+tmp*h, x+st44[i][3]*h );
  }
  
  tmp = st44[3][0] * k[0];
  for( j=1; j<4; j++ ){
    tmp += st44[3][j] * k[j];
  }

  x += h;
  y += tmp*h;
  return;

}

void RK4( double& y, double (*f)( const double&, const double& ),
	  double& x, const double& h )
{
  static int i,j;
  static double k[4];
  static double tmp;

  // k1
  k[0] = f( y, x );
  
  // k2 ... k4
  for( i=0; i<3; i++ ){
    tmp = st44[i][0] * k[0];
    for( j=1; j<(i+1); j++ ){
      tmp += st44[i][j] * k[j];
    }
    k[i+1] = f( y+tmp*h, x+st44[i][3]*h );
  }
  
  tmp = st44[3][0] * k[0];
  for( j=1; j<4; j++ ){
    tmp += st44[3][j] * k[j];
  }

  x += h;
  y += tmp*h;
  return;

}


void RK6( vector3& y, vector3 (*f)( const vector3, const double& ),
	  double& x, const double& h )
{
  static int i,j;
  static vector3 k[7];
  static vector3 tmp;
  
  // k1
  k[0] = f( y, x );

  // k2 ... k7
  for( i=0; i<6; i++ ){
    tmp = st76[i][0] * k[0];
    for( j=1; j<(i+1); j++ ){
      tmp += st76[i][j] * k[j];
      }
    k[i+1] = f( y+tmp*h, x+st76[i][6]*h );
  }
  
  tmp = st76[6][0] * k[0];
  for( j=1; j<7; j++ ){
    tmp += st76[6][j] * k[j];
  }

  x += h;
  y += tmp*h;
  return;

}


void RK6( double& y, double (*f)( const double&, const double& ),
	  double& x, const double& h )
{
  static int i,j;
  static double k[7];
  static double tmp;
  
  // k1
  k[0] = f( y, x );

  // k2 ... k7
  for( i=0; i<6; i++ ){
    tmp = st76[i][0] * k[0];
    for( j=1; j<(i+1); j++ ){
      tmp += st76[i][j] * k[j];
      }
    k[i+1] = f( y+tmp*h, x+st76[i][6]*h );
  }
  
  tmp = st76[6][0] * k[0];
  for( j=1; j<7; j++ ){
    tmp += st76[6][j] * k[j];
  }

  x += h;
  y += tmp*h;
  return;

}

# endif

// end
