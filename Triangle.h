#ifndef _TRIANGLE_H
#define _TRIANGLE_H

#include "Vector3.h"
#include "RigidTrans3.h"

// Triangle.h: interface for the triangle class.
//
//////////////////////////////////////////////////////////////////////

/* 
CLASS 
  Triangle

  Defines a triangle class for easy triangle handling
  (Includes transformation operator (|))2

KEYWORDS
  linear, algebra, triangle, transformation, baricenter

AUTHORS
  Yevgeny Menaker. (genmen@post.tau.ac.il)

  Copyright: SAMBA group, Tel-Aviv Univ. Israel, 2001.

CHANGES LOG
<UL>
<LI> Small changes, functionality is not changed. M.Shatsky 2002.
<LI> operator| computed trans incorrectly. Replaced. M.Shatsky 2003. 
<LI> area and sideLengths methods added Dina 2007.
</UL>

GOALS
  This class simplifies the work with triangle form
  allowing easy and fast operation with them. The core
  of this class is transformation operator, which performs
  3D transformation, transforming one triangle to fit the 
  second one effectively.
  
USAGE
  There are provided three operators: assignment(=), shifting the triangle
  by Vector3 (+) and transforming operator(|).
  In addition the Baricenter method is provided to calculate the 
  baricenter of triangle.

  A new triangle may be constructed with 3 3D vectors - one
  for each vertix of triangle
  EXAMPLE
    Vector3 v_A(1.0, 2.0, 3.0);
    Vector3 v_B(1.0, 5.0, 1.0);
    Vector3 v_C(0.0, 2.5, 2.0);

    triangle tr(v_A, v_B, v_C);

  END

  You may calculate the baricenter of given triangle:

  EXAMPLE
    Vector3 baricenter = tr.Baricenter();
  END

  Shifting of triangle:

  EXAMPLE
	triangle tr2 = tr + Vector3(1.0,1.0,1.0);
  END

  And finally the calculation of 3D transformation
  of two triangles, assuming they can be fit

  EXAMPLE
    RigidTrans3 rt3 = tr | tr2;

    cout << "The transformation is: " << rt3;

  END

*/
class Triangle  
{
public:
  //// Default constructor
  Triangle();

  //// Constructor with vertices initialization
  Triangle(const Vector3& a,const Vector3& b,const Vector3& c);

  //// Copy constructor
  Triangle(const Triangle& tr);

  virtual ~Triangle();
  
  //// Assignment operator
  Triangle& operator=(const Triangle& tr);
  
  //// Calculates the Baricenter of the triangle
  // We assume that every vertix of triangle has
  // the same weight (1)
  Vector3 baricenter() const;
  
  //// Calculates a shift of Baricenter of triangle (add vector)
  // (doesn't change the object itself)
  Triangle operator+(const Vector3& vec) const;
  
  //// Returns a transformation to fit triangle tr onto 
  // this triangle
  RigidTrans3 operator|(const Triangle& tr) const{	
    RigidTrans3 tr1(m_A,m_B,m_C);
    RigidTrans3 tr2(tr.m_A,tr.m_B,tr.m_C);
    
    return (tr1*(!tr2));
  }

  //// returns area of triangle formed by a b c
  static float area(const Vector3& a,const Vector3& b,const Vector3& c);

  //// returns side lengths of triangle formed by a b c
  static Vector3 sideLengths(const Vector3& a,const Vector3& b,const Vector3& c);

protected:
  // Triangle vertices
  Vector3 m_A;
  Vector3 m_B;
  Vector3 m_C;

};

#endif
