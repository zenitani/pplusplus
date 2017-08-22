//  -*- C++ -*-
//  electromagnetic test particle code         last updated : 2001/09/27

//
//  Copyright (C) 1998-2001
//             Seiji Zenitani <zenitani@space.eps.s.u-tokyo.ac.jp>
//
//  Solar Terrestrial Physics Group, Department of Geophysics,
//  Graduate School of Science, University of Tokyo
//  7-3-1, Hongo, Bunkyo-ku, Tokyo, 113-0033, JAPAN
//
//  You may copy, use, modify and redistribute this code
//  for ANY PURPOSE, without significant change, as long as
//  all copyright notice are retained.
//  The author provides this code `as is', and declares that
//  there is no warranty for it.
//
//
// *** History ***
//
// 1998/10/10  Ver 0.1   started
// 1999/03/10  Ver 1.0   stable release
// 2000/09/28      1.1   ready for relativistic motion
// 2001/09/27  Ver 1.5   integrated with relativistic version
// 

// *** Notice ***
//
//   vector3.h is required for vector operation.
//   RK.h      is required for advancing particles.
//
//   external force function,
//      vector3 F( const vector3& _r, const vector3& _v,
//                 const  double& _t, const  double&  q  );
//      must be defined in your program.


#ifndef _Z_PARTICLE_H_
#define _Z_PARTICLE_H_

#include <math.h>
#include <vector3.h>
#include <RK.h>


//
// force function F() prototype
//

vector3 F( const vector3&, const vector3&,
           const  double&, const  double& );


//
// particle class
//
// rk4(), rk6() ==> non-relativistic motion
// RK4(), RK6() ==> relativistic motion  ( c = 1.0 )
//

class particle
{

// m,q,t is PROTECTED variable.
// use functions "set/get(m,q,t)".
protected:
  double m, m_inv, q, t;

public:
  vector3 r, v;

  // constructor
  particle( void );
  
  // operator
  particle& operator  = ( const particle& );

  // get functions
  double  getm( void ) const;
  double  getq( void ) const;
  double  gett( void ) const;
  vector3 getr( void ) const;
  vector3 getv( void ) const;

  // set functions
  void setm( const double& );
  void setq( const double& );
  void sett( const double& );

  void setr( const vector3& );
  void setv( const vector3& );
  void setr( void );
  void setv( void );
  void setr( const double&, const double&, const double& );
  void setv( const double&, const double&, const double& );
  void setr( const int&, const int&, const int& );
  void setv( const int&, const int&, const int& );
  
  void reset( void );

  // non-relativistic
  void rk4( const double& );
  void rk6( const double& );

  // relativistic
  void RK4( const double& );
  void RK6( const double& );

};


// ---- constructor -----

particle::particle( void )
  : m(1.0), m_inv(1.0), q(0.0), t(0.0)
{
  r.set(); v.set();
}

// ---- operators ----
inline particle& particle::operator = ( const particle& p )
{
  m = p.m ; m_inv = p.m_inv;
  q = p.q ; t = p.t ;
  r = p.r ; v = p.v ;
  return *this;
}

// ---- member functions -----

inline double  particle::getm( void ) const{ return m; }
inline double  particle::getq( void ) const{ return q; }
inline double  particle::gett( void ) const{ return t; }
inline vector3 particle::getr( void ) const{ return r; }
inline vector3 particle::getv( void ) const{ return v; }

inline void particle::setm( const double& _m ){
  m = _m;  m_inv = 1.0 / _m;
}
inline void particle::setq( const double& _q = 0.0 ){ q = _q; }
inline void particle::sett( const double& _t = 0.0 ){ t = _t; }

inline void particle::setr( const vector3& _r ){ r = _r; }
inline void particle::setv( const vector3& _v ){ v = _v; }
inline void particle::setr( void ){ r.set(); }
inline void particle::setv( void ){ v.set(); }
inline void particle::setr( const double& _x,
                            const double& _y, const double& _z )
{
  r.set( _x, _y, _z );
}
inline void particle::setv( const double& _x,
                            const double& _y, const double& _z )
{
  v.set( _x, _y, _z );
}
inline void particle::setr( const int& _x, const int& _y, const int& _z )
{
  r.set( _x, _y, _z );
}
inline void particle::setv( const int& _x, const int& _y, const int& _z )
{
  v.set( _x, _y, _z );
}

inline void particle::reset( void ){ v.set(); r.set(); t=0.0; }


// proceed by Runge-Kutta methods
inline void particle::rk4( const double& h )
{
  register int i,j;
  static vector3 kr[4],kv[4];
  static vector3 tmpr, tmpv;

  // k1
  kr[0] = v;
  kv[0] = m_inv * F( r,v,t, q );

  // k2 ... k4
  for( i=0; i<3; i++ ){
    tmpr = st44[i][0] * kr[0];
    tmpv = st44[i][0] * kv[0];
    for( j=1; j<(i+1); j++ ){
      tmpr += st44[i][j] * kr[j];
      tmpv += st44[i][j] * kv[j];
    }
    kr[i+1] = v + tmpv*h;
    kv[i+1] = m_inv * F( r+tmpr*h,v+tmpv*h,t+st44[i][3]*h, q );
  }

  tmpr = st44[3][0] * kr[0];
  tmpv = st44[3][0] * kv[0];
  for( j=1; j<4; j++ ){
    tmpr += st44[3][j] * kr[j];
    tmpv += st44[3][j] * kv[j];
  }
  
  t += h;
  r += tmpr*h;
  v += tmpv*h;
  return;

}

// 6th order
inline void particle::rk6( const double& h )
{
  register int i,j;
  static vector3 kr[7],kv[7];
  static vector3 tmpr, tmpv;

  // k1
  kr[0] = v;
  kv[0] = m_inv * F( r,v,t, q );

  // k2 ... k7
  for( i=0; i<6; i++ ){
    tmpr = st76[i][0] * kr[0];
    tmpv = st76[i][0] * kv[0];
    for( j=1; j<(i+1); j++ ){
      tmpr += st76[i][j] * kr[j];
      tmpv += st76[i][j] * kv[j];
    }
    kr[i+1] = v + tmpv*h;
    kv[i+1] = m_inv * F( r+tmpr*h,v+tmpv*h,t+st76[i][6]*h, q );
  }
  
  tmpr = st76[6][0] * kr[0];
  tmpv = st76[6][0] * kv[0];
  for( j=1; j<7; j++ ){
    tmpr += st76[6][j] * kr[j];
    tmpv += st76[6][j] * kv[j];
  }

  t += h;
  r += tmpr*h;
  v += tmpv*h;
  return;

}

// relativistic motion
// proceed by Runge-Kutta methods
void particle::RK4( const double& h )
{
  register int i,j;
  static vector3 kr[4],kv[4];
  static vector3 tmpr, tmpv;

  // k1
  kr[0] = v.uv2v();
  kv[0] = m_inv * F( r,v,t, q );

  // k2 ... k4
  for( i=0; i<3; i++ ){
    tmpr = st44[i][0] * kr[0];
    tmpv = st44[i][0] * kv[0];
    for( j=1; j<(i+1); j++ ){
      tmpr += st44[i][j] * kr[j];
      tmpv += st44[i][j] * kv[j];
    }
    kr[i+1] = ( v+tmpv*h ).uv2v();
    kv[i+1] = m_inv * F( r+tmpr*h,v+tmpv*h,t+st44[i][3]*h, q );
  }

  tmpr = st44[3][0] * kr[0];
  tmpv = st44[3][0] * kv[0];
  for( j=1; j<4; j++ ){
    tmpr += st44[3][j] * kr[j];
    tmpv += st44[3][j] * kv[j];
  }
  
  t += h;
  r += tmpr*h;
  v += tmpv*h;
  return;

}

void particle::RK6( const double& h )
{
  register int i,j;
  static vector3 kr[7],kv[7];
  static vector3 tmpr, tmpv;

  // k1
  kr[0] = v.uv2v();
  kv[0] = m_inv * F( r,v,t, q );

  // k2 ... k7
  for( i=0; i<6; i++ ){
    tmpr = st76[i][0] * kr[0];
    tmpv = st76[i][0] * kv[0];
    for( j=1; j<(i+1); j++ ){
      tmpr += st76[i][j] * kr[j];
      tmpv += st76[i][j] * kv[j];
    }
    kr[i+1] = ( v+tmpv*h ).uv2v();
    kv[i+1] = m_inv * F( r+tmpr*h,v+tmpv*h,t+st76[i][6]*h, q );
  }
  
  tmpr = st76[6][0] * kr[0];
  tmpv = st76[6][0] * kv[0];
  for( j=1; j<7; j++ ){
    tmpr += st76[6][j] * kr[j];
    tmpv += st76[6][j] * kv[j];
  }

  t += h;
  r += tmpr*h;
  v += tmpv*h;
  return;

}

# endif

// end
