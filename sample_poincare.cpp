#include <particle.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* *********************************************************************
 Poincare map problem in a thin current sheet with a normal magnetic field.

 [1] S. Zenitani, I. Shinohara, T. Nagai, and T. Wada, Phys. Plasmas 20, 092120 (2013)
 [2] J. Chen and P. J. Palmadesso, J. Geophys. Res. 91, 1499 (1986)
 [3] J. Buchner and L. M. Zelenyi, J. Geophys. Res. 94, 11821 (1989)

 For general background, please consult Refs. [2,3].
For this specific implementation (noramlization etc.), please see
the Appendix chapter in Ref. [1] and references therein.
 ********************************************************************* */

// ************* initial parameters ************************************
// curvature parameter
const double kappa = 0.36178;
// number of particles
const int np = 256;
// ************* initial parameters ************************************

// electric field
vector3 E( const vector3& _r ){
  vector3 _E( 0.0, 0.0, 0.0 );
  return _E;
}
// magnetic field model
vector3 B( const vector3& _r ){
  vector3 _B( _r.z, 0.0, kappa );
  return _B;
}
// force
vector3 F( const vector3& _r, const vector3& _v,
           const  double& _t, const  double& _q ){
  return _q * ( E(_r) + _v * B(_r) );
}

int main()
{
  double dt = 0.01;
  particle p, pp, po;
  srand((unsigned) time(NULL)); // srand() may not be random enough on OSX/gcc

  // particle loop
  for( int ip=1; ip<=np; ip++ ){

    fprintf( stderr, "# running %d/%d th particle...\n", ip, np);

    // init
    p.sett(0);  p.setm(1);  p.setq(1);
    p.setr(0.0,0.0,0.0);
    p.setv(0.0,0.0,0.0);

    double PI2 = atan(1.0) * 8.0;
    double r1, r2;
    // ********** random scattering ***************
    r1 = (0.0+rand())/RAND_MAX;
    r2 = (0.0+rand())/RAND_MAX;
    p.v.x = 2*r1-1;
    r1 = sqrt( 1 - (p.v.x*p.v.x) );
    p.v.y = r1 * cos( PI2*r2 );
    p.v.z = r1 * sin( PI2*r2 );
    // ********** manual scattering ***************
//     r1 = 42.103; // very close to the fixed-point (parabolic)
//     //r1 = 45.0;
//     //r1 = 70.0;
//     r2 = PI2/360;
//     p.v.x =  0.0;
//     p.v.y = -sin(r1*r2);
//     p.v.z =  cos(r1*r2);
    // ***** adjusting the initial position *******
    p.r.x = -(1./kappa)*p.v.y;
    p.r.y = +(1./kappa)*p.v.x;

    // main loop
    for( int i=0; p.gett()<1000; i++ ){

//       if( i%20 == 0 ){
//      printf( "%f %f %f %f %f %f %f %f %d\n",
//              p.gett(),
//              p.r.x, p.r.y, p.r.z,
//              p.v.x, p.v.y, p.v.z,
//              p.v.abs2(), ip );
//       }

      // midplane crossing
      if( ( p.r.z * pp.r.z ) < 0.0 ){
        // linear interpolation
        // make sure that dt is small enough
        po.r = ( p.r.z * pp.r - pp.r.z * p.r ) / ( p.r.z - pp.r.z);
        po.v = ( p.r.z * pp.v - pp.r.z * p.v ) / ( p.r.z - pp.r.z);
        printf( "%f %f %f %f %f %f %d\n",
                po.r.x, po.r.y, po.r.z,
                po.v.x, po.v.y, po.v.z,
                ip );
      }
      pp = p;

      // nonrelativistic motion
      p.rk4(dt);

      // check the timestep
      if( ( B(p.r).abs() ) * dt > 0.3 ){
        fprintf( stderr, "# Exiting ... t = %lf", p.gett() );
        return -1;
      }
    }
  }

  return 0;
}
