///////////////////////////////////////////////////////////////////////////////
//
//	C3Vector.h
//
//	Header file for the C3Vector class
//	updated to use templates
//
//	Alan Brunton 2009
//
///////////////////////////////////////////////////////////////////////////////


#ifndef C3VECTOR_H
#define C3VECTOR_H


///////////////////////////////////////////////////////////////////////////////
//C3Vector
///////////////////////////////////////////////////////////////////////////////

template<typename scalar>
class C3Vector
{
public:

	scalar			x, y, z;

	C3Vector(): x(0.0), y(0.0), z(0.0) { }
	C3Vector(scalar sx, scalar sy, scalar sz): x(sx), y(sy), z(sz) { }
	C3Vector(scalar sx, scalar sy): x(sx), y(sy), z(0.0) { }
	C3Vector(scalar s): x(s), y(s), z(s) { }
	C3Vector(const C3Vector& v): x(v.x), y(v.y), z(v.z) { }
	C3Vector(scalar* pv): x(pv[0]), y(pv[1]), z(pv[2]) { }
	
	C3Vector& set(scalar sx, scalar sy, scalar sz)
	{
		x = sx;
		y = sy;
		z = sz;
		return *this;
	}
	
	C3Vector& set(scalar* pv)
	{
		x = pv[0];
		y = pv[1];
		z = pv[2];
		return *this;
	}

	
	C3Vector& operator=(const C3Vector& v)
	{
		x = v.x;
		y = v.y;
		z = v.z;
		return *this;
	}
	
	C3Vector& operator*=(const scalar s)
	{
		x *= s;
		y *= s;
		z *= s;
		return *this;
	}
	
	C3Vector& operator/=(const scalar s)
	{
		scalar rcpS = 0.0;
//		if (fabs(s) > g_epsilon)
		if (fabs(s) > (scalar)ABUTIL_EPSILON)
			rcpS = 1.0 / s;
		x *= rcpS;
		y *= rcpS;
		z *= rcpS;
		return *this;
	}
	
	C3Vector& operator+=(const C3Vector& v)
	{
		x += v.x;
		y += v.y;
		z += v.z;
		return *this;
	}
	
	C3Vector& operator-=(const C3Vector& v)
	{
		x -= v.x;
		y -= v.y;
		z -= v.z;
		return *this;
	}

	C3Vector operator*(const scalar s) const
	{
		C3Vector v(x, y, z);
		v *= s;
		return v;
	}
	
	C3Vector operator/(const scalar s) const
	{
		C3Vector v(x, y, z);
		v /= s;
		return v;
	}
	
	C3Vector operator+(const C3Vector& v) const
	{
		C3Vector result(x, y, z);
		result += v;
		return result;
	}
	
	C3Vector operator-(const C3Vector& v) const
	{
		C3Vector result(x, y, z);
		result -= v;
		return result;
	}
	
	C3Vector operator-() const
	{
		return C3Vector(-x, -y, -z);
	}

	scalar dot(const C3Vector& v) const
	{
		return x*v.x + y*v.y + z*v.z;
	}
	
	scalar lengthSquared() const
	{
		return dot(*this);
	}
	
	scalar length() const
	{
		return (scalar) ::sqrt(dot(*this));
	}

	
	C3Vector& normalize()
	{
		scalar rcpMag = 1.0 / length();
		x *= rcpMag;
		y *= rcpMag;
		z *= rcpMag;
		return *this;
	}

	C3Vector cross(const C3Vector& v) const
	{
		C3Vector result;
		result.x = y * v.z - z * v.y;
		result.y = z * v.x - x * v.z;
		result.z = x * v.y - y * v.x;
		return result;
	}

	void crossFast(const C3Vector& v, C3Vector& result) const
	{
		result.x = y * v.z - z * v.y;
		result.y = z * v.x - x * v.z;
		result.z = x * v.y - y * v.x;
	}

	bool operator==(const C3Vector& v)
	{
		C3Vector vdiff(x, y, z);
		vdiff -= v;
		return ::fabs(vdiff.x) < g_epsilon && ::fabs(vdiff.y) < g_epsilon && ::fabs(vdiff.z) < g_epsilon;
	}

	bool operator!=(const C3Vector& v)
	{
		C3Vector vdiff(x, y, z);
		vdiff -= v;
		return fabs(vdiff.x) >= g_epsilon || fabs(vdiff.y) >= g_epsilon || fabs(vdiff.z) >= g_epsilon;
	}
};

typedef C3Vector<float>		C3Vectorf;
typedef C3Vector<double>	C3Vectord;

#endif //C3VECTOR_H


