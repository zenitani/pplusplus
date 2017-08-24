//  -*- C++ -*-
//  3-dimensional vector class                last updated : 2017/08/24


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
// 1998/09/05  Ver 0.1   project started
// 1999/03/10  Ver 1.0   stable release
// 2000/03/06      1.1   relativistic functions
// 


#ifndef _Z_VECTOR3_H_
#define _Z_VECTOR3_H_

#include <math.h>


class vector3
{
public:
  double x, y, z;

  // constructor
  vector3( void );
  vector3( const double&, const double&, const double& );
  vector3( const int&, const int&, const int& );

  // operators
  vector3& operator  = ( const vector3& );
  vector3& operator += ( const vector3& );
  vector3& operator -= ( const vector3& );
  vector3& operator *= ( const double& );
  vector3& operator *= ( const int& );
  vector3& operator /= ( const double& );
  vector3& operator /= ( const int& );

  friend vector3 operator + ( const vector3& );
  friend vector3 operator - ( const vector3& );
  friend vector3 operator + ( const vector3&, const vector3& );
  friend vector3 operator - ( const vector3&, const vector3& );

  friend vector3 operator * ( const vector3&, const vector3& );
  friend double  operator % ( const vector3&, const vector3& );

  // friend operators
  friend vector3 operator * ( const double&, const vector3& );
  friend vector3 operator * ( const int&, const vector3& );
  friend vector3 operator * ( const vector3&, const double& );
  friend vector3 operator * ( const vector3&, const int& );
  friend vector3 operator / ( const vector3&, const double& );
  friend vector3 operator / ( const vector3&, const int& );

  // logical operators
  // friend int operator == ( const vector3&,  const vector3&  );
  // friend int operator != ( const vector3&,  const vector3&  );

// #ifdef _IOSTREAM_H
//   // stream operators
//   friend ostream &operator<< ( ostream &s, vector3 v );
//   friend istream &operator>> ( istream &s, vector3 v );
// #endif

  // member functions
  void set( const double&, const double&, const double& );
  void set( const int&, const int&, const int& );
  void set( void );
  void reset( void );

  double abs( void ) const;
  double abs2( void ) const;
  double gamma( void ) const;
  double ugamma( void ) const;
  vector3 v2uv( void ) const;
  vector3 uv2v( void ) const;

  // friend functions
  friend vector3 cross( const vector3&, const vector3& );
  friend double dot( const vector3&, const vector3& );
  friend double abs( const vector3& );

};


// ---- constructor ----

vector3::vector3( void )
  : x( 0.0 ), y( 0.0 ), z( 0.0 ) {}
vector3::vector3( const double& _x, const double& _y, const double& _z )
  : x( _x ), y( _y ), z( _z ) {}
vector3::vector3( const int& _x, const int& _y, const int& _z )
  : x(double(_x)), y(double(_y)), z(double(_z)) {}


// ----  operators ------

// unary operators
inline vector3& vector3::operator = ( const vector3& v )
{
  x = v.x ; y = v.y ; z = v.z ;
  return *this;
}

// assign operators
inline vector3& vector3::operator += ( const vector3& v )
{
  x += v.x ; y += v.y ; z += v.z ;
  return *this;
}
inline vector3& vector3::operator -= ( const vector3& v )
{
  x -= v.x ; y -= v.y ; z -= v.z ;
  return *this;
}

inline vector3& vector3::operator *= ( const double& d )
{
  x *= d ; y *= d ; z *= d ;
  return *this;
}
inline vector3& vector3::operator *= ( const int& i )
{
  double d = (double) i ;
  x *= d ; y *= d ; z *= d ;
  return *this;
}
inline vector3& vector3::operator /= ( const double& d )
{
  x /= d ; y /= d ; z /= d ;
  return *this;
}
inline vector3& vector3::operator /= ( const int& i )
{
  double d = (double) i ;
  x /= d ; y /= d ; z /= d ;
  return *this;
}

inline vector3 operator + ( const vector3& v )
{
  return v;
}
inline vector3 operator - ( const vector3& v )
{
  return vector3( -v.x, -v.y, -v.z );
}

// binary operators
inline vector3 operator + ( const vector3& a, const vector3& b )
{
  return vector3( a.x+b.x, a.y+b.y, a.z+b.z );
}
inline vector3 operator - ( const vector3& a, const vector3& b )
{
  return vector3( a.x-b.x, a.y-b.y, a.z-b.z );
}
inline vector3 operator * ( const vector3& a, const vector3& b )
{
  return vector3( a.y*b.z - a.z*b.y ,
                  a.z*b.x - a.x*b.z ,
                  a.x*b.y - a.y*b.x );
}
inline double  operator % ( const vector3& a, const vector3& b )
{
  return ( a.x*b.x + a.y*b.y + a.z*b.z );
}

// vector and scholar operations
inline vector3 operator * ( const double& d, const vector3& v )
{
  return vector3( d * v.x, d * v.y, d * v.z );
}
inline vector3 operator * ( const int& i, const vector3& v )
{
  double d = (double)i ;
  return vector3( d * v.x, d * v.y, d * v.z );
}
inline vector3 operator * ( const vector3& v, const double& d )
{
  return vector3( d * v.x, d * v.y, d * v.z );
}
inline vector3 operator * ( const vector3& v, const int& i )
{
  double d = (double)i ;
  return vector3( d * v.x, d * v.y, d * v.z );
}
inline vector3 operator / ( const vector3& v, const double& d )
{
  return vector3( v.x/d , v.y/d, v.z/d );
}
inline vector3 operator / ( const vector3& v, const int& i )
{
  double d = (double)i ;
  return vector3( v.x/d , v.y/d, v.z/d );
}

// logical operators
// int operator == ( const vector3& a, const vector3& b )
// {
//   return (( a.x == b.x ) && ( a.y == b.y ) && ( a.z == b.z ));
// }
// int operator != ( const vector3& a, const vector3& b )
// {
//   return (( a.x != b.x ) || ( a.y != b.y ) || ( a.z != b.z ));
// }

// #ifdef _IOSTREAM_H
// // stream operators
// ostream &operator<< ( ostream &s, vector3 v )
// {
//   s << v.x << ", " << v.y << ", " << v.z;
//   return s;
// }

// istream &operator>> ( istream &s, vector3 &v )
// {
//   s >> v.x >> v.y >> v.z;
//   return s;
// }
// #endif


// ---- member functions ------

void vector3::set( const double& _x, const double& _y, const double& _z )
{
  x = _x ; y = _y ; z = _z ;
}
void vector3::set( const int& _x, const int& _y, const int& _z )
{
  x = (double)_x ; y = (double)_y ; z = (double)_z ;
}
void vector3::set( void ){   x = 0.0 ; y = 0.0 ; z = 0.0 ; }
void vector3::reset( void ){ x = 0.0 ; y = 0.0 ; z = 0.0 ; }

inline double vector3::abs2( void ) const
{
  return( x*x + y*y + z*z );
}
inline double vector3::abs( void ) const
{
  return sqrt( x*x + y*y + z*z );
}
inline double vector3::gamma( void ) const
{
  return 1.0 / sqrt( 1.0 - x*x - y*y - z*z );
}
inline double vector3::ugamma( void ) const
{
  return sqrt( 1.0 + x*x + y*y + z*z );
}
inline vector3 vector3::v2uv( void ) const
{
  double f = 1.0 / sqrt( 1.0 - x*x - y*y - z*z );
  return vector3( f*x, f*y, f*z );
}
inline vector3 vector3::uv2v( void ) const
{
  double f = 1.0 / sqrt( 1.0 + x*x + y*y + z*z );
  return vector3( f*x, f*y, f*z );
}


// ---- useful functions ------

vector3 cross( const vector3& a, const vector3& b )
{
  return( a * b );
}
double dot( const vector3& a, const vector3& b )
{
  return( a % b );
}
double abs( const vector3& v )
{
  return( v.abs() );
}

# endif

// end
