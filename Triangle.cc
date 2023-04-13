#include "Triangle.h"
#include "Matrix3.h"
#include <math.h>

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Triangle::Triangle()
{

}

Triangle::Triangle(const Vector3& a,const Vector3& b,const Vector3& c)
{
	m_A = a;
	m_B = b;
	m_C = c;
}

// Copy constructor
Triangle::Triangle(const Triangle& tr)
{
	m_A = tr.m_A;
	m_B = tr.m_B;
	m_C = tr.m_C;
}

// Destructor
Triangle::~Triangle()
{

}

// Assignment operator
Triangle& Triangle::operator =(const Triangle& tr)
{
	m_A = tr.m_A;
	m_B = tr.m_B;
	m_C = tr.m_C;

	return(*this);
}


// Calculate the Baricenter of triangle
Vector3 Triangle::baricenter()const
{
	Vector3 v;

	v = (m_A + m_B + m_C)*((Vector3::real)1/(Vector3::real)3);

	return(v);
}

// Shift the Baricenter (+)
Triangle Triangle::operator +(const Vector3& vec) const
{
	Triangle tr;

	tr.m_A = m_A + vec;
	tr.m_B = m_B + vec;
	tr.m_C = m_C + vec;

	return(tr);
}

float Triangle::area(const Vector3& a,const Vector3& b,const Vector3& c) {
  Vector3 ab = a-b;
  Vector3 cb = c-b;
  
  // Check for flat triangle case:
  // To calculate the area of the triangle, we use the following
  // property of vector multiplication, namely, the norm of the
  // orthogonal vector is the area of the parallelogram defined by the
  // two vectors. Thus, the area of the basis triangle is half of the
  // norm of the cross product of two of its sides.    
  return (ab&cb).norm() / 2;
}

// returns side lengths of triangle formed by a b c
Vector3 Triangle::sideLengths(const Vector3& a,const Vector3& b,const Vector3& c) {
  Vector3 ab = a-b;
  Vector3 cb = c-b;
  Vector3 ac = a-c;
  return Vector3(ab.norm(), cb.norm(), ac.norm());
}
