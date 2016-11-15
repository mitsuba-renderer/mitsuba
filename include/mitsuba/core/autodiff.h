/**
    Automatic differentiation data type for C++, depends on the Eigen
	linear algebra library.

    Copyright (c) 2012 by Wenzel Jakob. Based on code by Jon Kaldor
    and Eitan Grinspun.

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
*/

#if !defined(__AUTODIFF_H)
#define __AUTODIFF_H
#define EIGEN_DONT_PARALLELIZE
#define EIGEN_NO_DEBUG

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cmath>
#include <stdexcept>

#if defined(__WINDOWS__)
  #define __tls_decl __declspec(thread)
#else
  #define __tls_decl __thread
#endif

/**
 * \brief Base class of all automatic differentiation types
 *
 * This class records the number of independent variables with respect
 * to which derivatives are computed.
 */
struct DiffScalarBase {
	// ======================================================================
	/// @{ \name Configuration
	// ======================================================================

	/**
	 * \brief Set the independent variable count used by the automatic
	 * differentiation layer
	 *
	 * This function must be called before doing any computations with
	 * \ref DScalar1 or \ref DScalar2. The value will be recorded in
	 * thread-local storage.
	 */
	static inline void setVariableCount(size_t value) {
		m_variableCount = value;
	}

	/// Get the variable count used by the automatic differentiation layer
	static inline size_t getVariableCount() {
		return m_variableCount;
	}

	/// @}
	// ======================================================================

	static __tls_decl size_t m_variableCount;
};

#define DECLARE_DIFFSCALAR_BASE() \
	__tls_decl size_t DiffScalarBase::m_variableCount = 0

/**
 * \brief Automatic differentiation scalar with first-order derivatives
 *
 * This class provides an instrumented "scalar" value, which may be dependent on
 * a number of independent variables. The implementation keeps tracks of
 * first -order drivatives with respect to these variables using a set
 * of overloaded operations and implementations of special functions (sin,
 * tan, exp, ..).
 *
 * This is extremely useful for numerical zero-finding, particularly when
 * analytic derivatives from programs like Maple or Mathematica suffer from
 * excessively complicated expressions.
 *
 * The class relies on templates, which makes it possible to fix the
 * number of independent variables at compile-time so that instances can
 * be allocated on the stack. Otherwise, they will be placed on the heap.
 *
 * This is an extended C++ port of Jon Kaldor's implementation, which is
 * based on a C version by Eitan Grinspun at Caltech)
 *
 * \sa DScalar2
 * \author Wenzel Jakob
 */
template <typename _Scalar, typename _Gradient = Eigen::Matrix<_Scalar, Eigen::Dynamic, 1> >
	struct DScalar1 : public DiffScalarBase {
public:
	typedef _Scalar                                         Scalar;
	typedef _Gradient                                       Gradient;
	typedef Eigen::Matrix<DScalar1, 2, 1>                   DVector2;
	typedef Eigen::Matrix<DScalar1, 3, 1>                   DVector3;

	// ======================================================================
	/// @{ \name Constructors and accessors
	// ======================================================================

	/// Create a new constant automatic differentiation scalar
	explicit DScalar1(Scalar value = (Scalar) 0) : value(value) {
		size_t variableCount = getVariableCount();
		grad.resize(Eigen::NoChange_t(), variableCount);
		grad.setZero();
	}

	/// Construct a new scalar with the specified value and one first derivative set to 1
	DScalar1(size_t index, const Scalar &value)
	 : value(value) {
		size_t variableCount = getVariableCount();
		grad.resize(Eigen::NoChange_t(), variableCount);
		grad.setZero();
		grad(index) = 1;
	}

	/// Construct a scalar associated with the given gradient
	DScalar1(Scalar value, const Gradient &grad)
	 : value(value), grad(grad) { }

	/// Copy constructor
	DScalar1(const DScalar1 &s)
	 : value(s.value), grad(s.grad) { }

	inline const Scalar &getValue() const { return value; }
	inline const Gradient &getGradient() const { return grad; }

	// ======================================================================
	/// @{ \name Addition
	// ======================================================================
	friend DScalar1 operator+(const DScalar1 &lhs, const DScalar1 &rhs) {
		return DScalar1(lhs.value+rhs.value, lhs.grad+rhs.grad);
	}

	friend DScalar1 operator+(const DScalar1 &lhs, const Scalar &rhs) {
		return DScalar1(lhs.value+rhs, lhs.grad);
	}

	friend DScalar1 operator+(const Scalar &lhs, const DScalar1 &rhs) {
		return DScalar1(rhs.value+lhs, rhs.grad);
	}

	inline DScalar1& operator+=(const DScalar1 &s) {
		value += s.value;
		grad += s.grad;
		return *this;
	}

	inline DScalar1& operator+=(const Scalar &v) {
		value += v;
		return *this;
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Subtraction
	// ======================================================================

	friend DScalar1 operator-(const DScalar1 &lhs, const DScalar1 &rhs) {
		return DScalar1(lhs.value-rhs.value, lhs.grad-rhs.grad);
	}

	friend DScalar1 operator-(const DScalar1 &lhs, const Scalar &rhs) {
		return DScalar1(lhs.value-rhs, lhs.grad);
	}

	friend DScalar1 operator-(const Scalar &lhs, const DScalar1 &rhs) {
		return DScalar1(lhs-rhs.value, -rhs.grad);
	}

	friend DScalar1 operator-(const DScalar1 &s) {
		return DScalar1(-s.value, -s.grad);
	}

	inline DScalar1& operator-=(const DScalar1 &s) {
		value -= s.value;
		grad -= s.grad;
		return *this;
	}

	inline DScalar1& operator-=(const Scalar &v) {
		value -= v;
		return *this;
	}
	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Division
	// ======================================================================
	friend DScalar1 operator/(const DScalar1 &lhs, const Scalar &rhs) {
		if (rhs == 0)
			throw std::runtime_error("DScalar1: Division by zero!");
		Scalar inv = 1.0f / rhs;
		return DScalar1(lhs.value*inv, lhs.grad*inv);
	}

	friend DScalar1 operator/(const Scalar &lhs, const DScalar1 &rhs) {
		return lhs * inverse(rhs);
	}

	friend DScalar1 operator/(const DScalar1 &lhs, const DScalar1 &rhs) {
		return lhs * inverse(rhs);
	}

	friend DScalar1 inverse(const DScalar1 &s) {
		Scalar valueSqr = s.value*s.value,
			invValueSqr = (Scalar) 1 / valueSqr;

		// vn = 1/v, Dvn = -1/(v^2) Dv
		return DScalar1((Scalar) 1 / s.value, s.grad * -invValueSqr);
	}

	inline DScalar1& operator/=(const Scalar &v) {
		value /= v;
		grad /= v;
		return *this;
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Multiplication
	// ======================================================================
	inline friend DScalar1 operator*(const DScalar1 &lhs, const Scalar &rhs) {
		return DScalar1(lhs.value*rhs, lhs.grad*rhs);
	}

	inline friend DScalar1 operator*(const Scalar &lhs, const DScalar1 &rhs) {
		return DScalar1(rhs.value*lhs, rhs.grad*lhs);
	}

	inline friend DScalar1 operator*(const DScalar1 &lhs, const DScalar1 &rhs) {
		// Product rule
		return DScalar1(lhs.value*rhs.value,
			rhs.grad * lhs.value + lhs.grad * rhs.value);
	}

	inline DScalar1& operator*=(const Scalar &v) {
		value *= v;
		grad *= v;
		return *this;
	}
	
	inline DScalar1& operator*=(const DScalar1 &v) {
		grad = v.grad * value + grad * v.value;
		value *= v.value;
		return *this;
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Miscellaneous functions
	// ======================================================================

	friend DScalar1 sqrt(const DScalar1 &s) {
		Scalar sqrtVal = std::sqrt(s.value),
		       temp    = (Scalar) 1 / ((Scalar) 2 * sqrtVal);

		// vn = sqrt(v)
		// Dvn = 1/(2 sqrt(v)) Dv
		return DScalar1(sqrtVal, s.grad * temp);
	}

	friend DScalar1 pow(const DScalar1 &s, const Scalar &a) {
		Scalar powVal = std::pow(s.value, a),
		       temp   = a * std::pow(s.value, a-1);
		// vn = v ^ a, Dvn = a*v^(a-1) * Dv
		return DScalar1(powVal, s.grad * temp);
	}

	friend DScalar1 exp(const DScalar1 &s) {
		Scalar expVal = std::exp(s.value);

		// vn = exp(v), Dvn = exp(v) * Dv
		return DScalar1(expVal, s.grad * expVal);
	}

	friend DScalar1 log(const DScalar1 &s) {
		Scalar logVal = std::log(s.value);

		// vn = log(v), Dvn = Dv / v
		return DScalar1(logVal, s.grad / s.value);
	}

	friend DScalar1 sin(const DScalar1 &s) {
		// vn = sin(v), Dvn = cos(v) * Dv
		return DScalar1(std::sin(s.value), s.grad * std::cos(s.value));
	}

	friend DScalar1 cos(const DScalar1 &s) {
		// vn = cos(v), Dvn = -sin(v) * Dv
		return DScalar1(std::cos(s.value), s.grad * -std::sin(s.value));
	}

	friend DScalar1 acos(const DScalar1 &s) {
		if (std::abs(s.value) >= 1)
			throw std::runtime_error("acos: Expected a value in (-1, 1)");

		Scalar temp = -std::sqrt((Scalar) 1 - s.value*s.value);

		// vn = acos(v), Dvn = -1/sqrt(1-v^2) * Dv
		return DScalar1(std::acos(s.value),
			s.grad * ((Scalar) 1 / temp));
	}

	friend DScalar1 asin(const DScalar1 &s) {
		if (std::abs(s.value) >= 1)
			throw std::runtime_error("asin: Expected a value in (-1, 1)");

		Scalar temp = std::sqrt((Scalar) 1 - s.value*s.value);

		// vn = asin(v), Dvn = 1/sqrt(1-v^2) * Dv
		return DScalar1(std::asin(s.value),
			s.grad * ((Scalar) 1 / temp));
	}

	friend DScalar1 atan2(const DScalar1 &y, const DScalar1 &x) {
		Scalar denom = x.value*x.value + y.value*y.value;

		// vn = atan2(y, x), Dvn = (x*Dy - y*Dx) / (x^2 + y^2)
		return DScalar1(std::atan2(y.value, x.value),
			y.grad * (x.value / denom) - x.grad * (y.value / denom));
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Comparison and assignment
	// ======================================================================

	inline void operator=(const DScalar1& s) { value = s.value; grad = s.grad; }
	inline void operator=(const Scalar &v) { value = v; grad.setZero(); }
	inline bool operator<(const DScalar1& s) const { return value < s.value; }
	inline bool operator<=(const DScalar1& s) const { return value <= s.value; }
	inline bool operator>(const DScalar1& s) const { return value > s.value; }
	inline bool operator>=(const DScalar1& s) const { return value >= s.value; }
	inline bool operator<(const Scalar& s) const { return value < s; }
	inline bool operator<=(const Scalar& s) const { return value <= s; }
	inline bool operator>(const Scalar& s) const { return value > s; }
	inline bool operator>=(const Scalar& s) const { return value >= s; }
	inline bool operator==(const Scalar& s) const { return value == s; }
	inline bool operator!=(const Scalar& s) const { return value != s; }

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Comparison and assignment
	// ======================================================================

	#if defined(__MITSUBA_MITSUBA_H_) /* Mitsuba-specific */
		/// Initialize a constant two-dimensional vector
		static inline DVector2 vector(const mitsuba::TVector2<Scalar> &v) {
			return DVector2(DScalar1(v.x), DScalar1(v.y));
		}

		/// Initialize a constant two-dimensional vector
		static inline DVector2 vector(const mitsuba::TPoint2<Scalar> &p) {
			return DVector2(DScalar1(p.x), DScalar1(p.y));
		}

		/// Create a constant three-dimensional vector
		static inline DVector3 vector(const mitsuba::TVector3<Scalar> &v) {
			return DVector3(DScalar1(v.x), DScalar1(v.y), DScalar1(v.z));
		}

		/// Create a constant three-dimensional vector
		static inline DVector3 vector(const mitsuba::TPoint3<Scalar> &p) {
			return DVector3(DScalar1(p.x), DScalar1(p.y), DScalar1(p.z));
		}
	#endif

	/// @}
	// ======================================================================
protected:
	Scalar value;
	Gradient grad;
};

template <typename Scalar, typename VecType>
		std::ostream &operator<<(std::ostream &out, const DScalar1<Scalar, VecType> &s) {
	out << "[" << s.getValue()
		<< ", grad=" << s.getGradient().format(Eigen::IOFormat(4, 1, ", ", "; ", "", "", "[", "]"))
		<< "]";
	return out;
}

/**
 * \brief Automatic differentiation scalar with first- and second-order derivatives
 *
 * This class provides an instrumented "scalar" value, which may be dependent on
 * a number of independent variables. The implementation keeps tracks of first
 * and second-order drivatives with respect to these variables using a set
 * of overloaded operations and implementations of special functions (sin,
 * tan, exp, ..).
 *
 * This is extremely useful for numerical optimization, particularly when
 * analytic derivatives from programs like Maple or Mathematica suffer from
 * excessively complicated expressions.
 *
 * The class relies on templates, which makes it possible to fix the
 * number of independent variables at compile-time so that instances can
 * be allocated on the stack. Otherwise, they will be placed on the heap.
 *
 * This is an extended C++ port of Jon Kaldor's implementation, which is
 * based on a C version by Eitan Grinspun at Caltech)
 *
 * \sa DScalar1
 * \author Wenzel Jakob
 */
template <typename _Scalar, typename _Gradient = Eigen::Matrix<_Scalar, Eigen::Dynamic, 1>,
		  typename _Hessian = Eigen::Matrix<_Scalar, Eigen::Dynamic, Eigen::Dynamic> >
		  struct DScalar2 : public DiffScalarBase {
public:
	typedef _Scalar                                         Scalar;
	typedef _Gradient                                       Gradient;
	typedef _Hessian                                        Hessian;
	typedef Eigen::Matrix<DScalar2, 2, 1>                   DVector2;
	typedef Eigen::Matrix<DScalar2, 3, 1>                   DVector3;

	// ======================================================================
	/// @{ \name Constructors and accessors
	// ======================================================================

	/// Create a new constant automatic differentiation scalar
	explicit DScalar2(Scalar value = (Scalar) 0) : value(value) {
		size_t variableCount = getVariableCount();

		grad.resize(Eigen::NoChange_t(), variableCount);
		grad.setZero();
		hess.resize(variableCount, variableCount);
		hess.setZero();
	}

	/// Construct a new scalar with the specified value and one first derivative set to 1
	DScalar2(size_t index, const Scalar &value)
	 : value(value) {
		size_t variableCount = getVariableCount();

		grad.resize(Eigen::NoChange_t(), variableCount);
		grad.setZero();
		grad(index) = 1;
		hess.resize(variableCount, variableCount);
		hess.setZero();
	}

	/// Construct a scalar associated with the given gradient and Hessian
	DScalar2(Scalar value, const Gradient &grad, const Hessian &hess)
	 : value(value), grad(grad), hess(hess) { }

	/// Copy constructor
	DScalar2(const DScalar2 &s)
	 : value(s.value), grad(s.grad), hess(s.hess) { }

	inline const Scalar &getValue() const { return value; }
	inline const Gradient &getGradient() const { return grad; }
	inline const Hessian &getHessian() const { return hess; }

	// ======================================================================
	/// @{ \name Addition
	// ======================================================================
	friend DScalar2 operator+(const DScalar2 &lhs, const DScalar2 &rhs) {
		return DScalar2(lhs.value+rhs.value,
			lhs.grad+rhs.grad, lhs.hess+rhs.hess);
	}

	friend DScalar2 operator+(const DScalar2 &lhs, const Scalar &rhs) {
		return DScalar2(lhs.value+rhs, lhs.grad, lhs.hess);
	}

	friend DScalar2 operator+(const Scalar &lhs, const DScalar2 &rhs) {
		return DScalar2(rhs.value+lhs, rhs.grad, rhs.hess);
	}

	inline DScalar2& operator+=(const DScalar2 &s) {
		value += s.value;
		grad += s.grad;
		hess += s.hess;
		return *this;
	}

	inline DScalar2& operator+=(const Scalar &v) {
		value += v;
		return *this;
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Subtraction
	// ======================================================================

	friend DScalar2 operator-(const DScalar2 &lhs, const DScalar2 &rhs) {
		return DScalar2(lhs.value-rhs.value, lhs.grad-rhs.grad, lhs.hess-rhs.hess);
	}

	friend DScalar2 operator-(const DScalar2 &lhs, const Scalar &rhs) {
		return DScalar2(lhs.value-rhs, lhs.grad, lhs.hess);
	}

	friend DScalar2 operator-(const Scalar &lhs, const DScalar2 &rhs) {
		return DScalar2(lhs-rhs.value, -rhs.grad, -rhs.hess);
	}

	friend DScalar2 operator-(const DScalar2 &s) {
		return DScalar2(-s.value, -s.grad, -s.hess);
	}

	inline DScalar2& operator-=(const DScalar2 &s) {
		value -= s.value;
		grad -= s.grad;
		hess -= s.hess;
		return *this;
	}

	inline DScalar2& operator-=(const Scalar &v) {
		value -= v;
		return *this;
	}
	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Division
	// ======================================================================
	friend DScalar2 operator/(const DScalar2 &lhs, const Scalar &rhs) {
		if (rhs == 0)
			throw std::runtime_error("DScalar2: Division by zero!");
		Scalar inv = 1.0f / rhs;
		return DScalar2(lhs.value*inv, lhs.grad*inv, lhs.hess*inv);
	}

	friend DScalar2 operator/(const Scalar &lhs, const DScalar2 &rhs) {
		return lhs * inverse(rhs);
	}

	friend DScalar2 operator/(const DScalar2 &lhs, const DScalar2 &rhs) {
		return lhs * inverse(rhs);
	}

	friend DScalar2 inverse(const DScalar2 &s) {
		Scalar valueSqr = s.value*s.value,
			valueCub = valueSqr * s.value,
			invValueSqr = (Scalar) 1 / valueSqr;

		// vn = 1/v
		DScalar2 result((Scalar) 1 / s.value);

		// Dvn = -1/(v^2) Dv
		result.grad = s.grad * -invValueSqr;

		// D^2vn = -1/(v^2) D^2v + 2/(v^3) Dv Dv^T
		result.hess = s.hess * -invValueSqr;
		result.hess += s.grad * s.grad.transpose()
			* ((Scalar) 2 / valueCub);

		return result;
	}

	inline DScalar2& operator/=(const Scalar &v) {
		value /= v;
		grad /= v;
		hess /= v;
		return *this;
	}
	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Multiplication
	// ======================================================================
	friend DScalar2 operator*(const DScalar2 &lhs, const Scalar &rhs) {
		return DScalar2(lhs.value*rhs, lhs.grad*rhs, lhs.hess*rhs);
	}

	friend DScalar2 operator*(const Scalar &lhs, const DScalar2 &rhs) {
		return DScalar2(rhs.value*lhs, rhs.grad*lhs, rhs.hess*lhs);
	}

	friend DScalar2 operator*(const DScalar2 &lhs, const DScalar2 &rhs) {
		DScalar2 result(lhs.value*rhs.value);

		/// Product rule
		result.grad = rhs.grad * lhs.value + lhs.grad * rhs.value;

		// (i,j) = g*F_xixj + g*G_xixj + F_xi*G_xj + F_xj*G_xi
		result.hess = rhs.hess * lhs.value;
		result.hess += lhs.hess * rhs.value;
		result.hess += lhs.grad * rhs.grad.transpose();
		result.hess += rhs.grad * lhs.grad.transpose();

		return result;
	}

	inline DScalar2& operator*=(const Scalar &v) {
		value *= v;
		grad *= v;
		hess *= v;
		return *this;
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Miscellaneous functions
	// ======================================================================

	friend DScalar2 sqrt(const DScalar2 &s) {
		Scalar sqrtVal = std::sqrt(s.value),
		       temp    = (Scalar) 1 / ((Scalar) 2 * sqrtVal);

		// vn = sqrt(v)
		DScalar2 result(sqrtVal);

		// Dvn = 1/(2 sqrt(v)) Dv
		result.grad = s.grad * temp;

		// D^2vn = 1/(2 sqrt(v)) D^2v - 1/(4 v*sqrt(v)) Dv Dv^T
		result.hess = s.hess * temp;
		result.hess += s.grad * s.grad.transpose()
			* (-(Scalar) 1 / ((Scalar) 4 * s.value * sqrtVal));

		return result;
	}

	friend DScalar2 pow(const DScalar2 &s, const Scalar &a) {
		Scalar powVal = std::pow(s.value, a),
		       temp   = a * std::pow(s.value, a-1);
		// vn = v ^ a
		DScalar2 result(powVal);

		// Dvn = a*v^(a-1) * Dv
		result.grad = s.grad * temp;

		// D^2vn = a*v^(a-1) D^2v - 1/(4 v*sqrt(v)) Dv Dv^T
		result.hess = s.hess * temp;
		result.hess += s.grad * s.grad.transpose()
			* (a * (a-1) * std::pow(s.value, a-2));

		return result;
	}

	friend DScalar2 exp(const DScalar2 &s) {
		Scalar expVal = std::exp(s.value);

		// vn = exp(v)
		DScalar2 result(expVal);

		// Dvn = exp(v) * Dv
		result.grad = s.grad * expVal;

		// D^2vn = exp(v) * Dv*Dv^T + exp(v) * D^2v
		result.hess = (s.grad * s.grad.transpose()
			+ s.hess) * expVal;

		return result;
	}

	friend DScalar2 log(const DScalar2 &s) {
		Scalar logVal = std::log(s.value);

		// vn = log(v)
		DScalar2 result(logVal);

		// Dvn = Dv / v
		result.grad = s.grad / s.value;

		// D^2vn = (v*D^2v - Dv*Dv^T)/(v^2)
		result.hess = s.hess / s.value -
			(s.grad * s.grad.transpose() / (s.value*s.value));

		return result;
	}

	friend DScalar2 sin(const DScalar2 &s) {
		Scalar sinVal = std::sin(s.value),
		       cosVal = std::cos(s.value);

		// vn = sin(v)
		DScalar2 result(sinVal);

		// Dvn = cos(v) * Dv
		result.grad = s.grad * cosVal;

		// D^2vn = -sin(v) * Dv*Dv^T + cos(v) * Dv^2
		result.hess = s.hess * cosVal;
		result.hess += s.grad * s.grad.transpose() * -sinVal;

		return result;
	}

	friend DScalar2 cos(const DScalar2 &s) {
		Scalar sinVal = std::sin(s.value),
		       cosVal = std::cos(s.value);
		// vn = cos(v)
		DScalar2 result(cosVal);

		// Dvn = -sin(v) * Dv
		result.grad = s.grad * -sinVal;

		// D^2vn = -cos(v) * Dv*Dv^T - sin(v) * Dv^2
		result.hess = s.hess * -sinVal;
		result.hess += s.grad * s.grad.transpose() * -cosVal;

		return result;
	}

	friend DScalar2 acos(const DScalar2 &s) {
		if (std::abs(s.value) >= 1)
			throw std::runtime_error("acos: Expected a value in (-1, 1)");

		Scalar temp = -std::sqrt((Scalar) 1 - s.value*s.value);

		// vn = acos(v)
		DScalar2 result(std::acos(s.value));

		// Dvn = -1/sqrt(1-v^2) * Dv
		result.grad = s.grad * ((Scalar) 1 / temp);

		// D^2vn = -1/sqrt(1-v^2) * D^2v - v/[(1-v^2)^(3/2)] * Dv*Dv^T
		result.hess = s.hess * ((Scalar) 1 / temp);
		result.hess += s.grad * s.grad.transpose()
			* s.value / (temp*temp*temp);

		return result;
	}

	friend DScalar2 asin(const DScalar2 &s) {
		if (std::abs(s.value) >= 1)
			throw std::runtime_error("asin: Expected a value in (-1, 1)");

		Scalar temp = std::sqrt((Scalar) 1 - s.value*s.value);

		// vn = asin(v)
		DScalar2 result(std::asin(s.value));

		// Dvn = 1/sqrt(1-v^2) * Dv
		result.grad = s.grad * ((Scalar) 1 / temp);

		// D^2vn = 1/sqrt(1-v*v) * D^2v + v/[(1-v^2)^(3/2)] * Dv*Dv^T
		result.hess = s.hess * ((Scalar) 1 / temp);
		result.hess += s.grad * s.grad.transpose()
			* s.value / (temp*temp*temp);

		return result;
	}

	friend DScalar2 atan2(const DScalar2 &y, const DScalar2 &x) {
		// vn = atan2(y, x)
		DScalar2 result(std::atan2(y.value, x.value));

		// Dvn = (x*Dy - y*Dx) / (x^2 + y^2)
		Scalar denom = x.value*x.value + y.value*y.value,
			denomSqr = denom*denom;
		result.grad = y.grad * (x.value / denom)
			- x.grad * (y.value / denom);

		// D^2vn = (Dy*Dx^T + xD^2y - Dx*Dy^T - yD^2x) / (x^2+y^2)
		//    - [(x*Dy - y*Dx) * (2*x*Dx + 2*y*Dy)^T] / (x^2+y^2)^2
		result.hess = (y.hess*x.value
			+ y.grad * x.grad.transpose()
			- x.hess*y.value
			- x.grad*y.grad.transpose()
		) / denom;

		result.hess -=
			(y.grad*(x.value/denomSqr) - x.grad*(y.value/denomSqr)) *
			(x.grad*((Scalar) 2 * x.value) + y.grad*((Scalar) 2 * y.value)).transpose();

		return result;
	}

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Comparison and assignment
	// ======================================================================

	inline void operator=(const DScalar2& s) { value = s.value; grad = s.grad; hess = s.hess; }
	inline void operator=(const Scalar &v) { value = v; grad.setZero(); hess.setZero(); }
	inline bool operator<(const DScalar2& s) const { return value < s.value; }
	inline bool operator<=(const DScalar2& s) const { return value <= s.value; }
	inline bool operator>(const DScalar2& s) const { return value > s.value; }
	inline bool operator>=(const DScalar2& s) const { return value >= s.value; }
	inline bool operator<(const Scalar& s) const { return value < s; }
	inline bool operator<=(const Scalar& s) const { return value <= s; }
	inline bool operator>(const Scalar& s) const { return value > s; }
	inline bool operator>=(const Scalar& s) const { return value >= s; }
	inline bool operator==(const Scalar& s) const { return value == s; }
	inline bool operator!=(const Scalar& s) const { return value != s; }

	/// @}
	// ======================================================================

	// ======================================================================
	/// @{ \name Comparison and assignment
	// ======================================================================

	#if defined(__MITSUBA_MITSUBA_H_) /* Mitsuba-specific */
		/// Initialize a constant two-dimensional vector
		static inline DVector2 vector(const mitsuba::TVector2<Scalar> &v) {
			return DVector2(DScalar2(v.x), DScalar2(v.y));
		}

		/// Initialize a constant two-dimensional vector
		static inline DVector2 vector(const mitsuba::TPoint2<Scalar> &p) {
			return DVector2(DScalar2(p.x), DScalar2(p.y));
		}

		/// Create a constant three-dimensional vector
		static inline DVector3 vector(const mitsuba::TVector3<Scalar> &v) {
			return DVector3(DScalar2(v.x), DScalar2(v.y), DScalar2(v.z));
		}

		/// Create a constant three-dimensional vector
		static inline DVector3 vector(const mitsuba::TPoint3<Scalar> &p) {
			return DVector3(DScalar2(p.x), DScalar2(p.y), DScalar2(p.z));
		}
	#endif

	/// @}
	// ======================================================================
protected:
	Scalar value;
	Gradient grad;
	Hessian hess;
};

template <typename Scalar, typename VecType, typename MatType>
		std::ostream &operator<<(std::ostream &out, const DScalar2<Scalar, VecType, MatType> &s) {
	out << "[" << s.getValue()
		<< ", grad=" << s.getGradient().format(Eigen::IOFormat(4, 1, ", ", "; ", "", "", "[", "]"))
		<< ", hess=" << s.getHessian().format(Eigen::IOFormat(4, 0, ", ", "; ", "", "", "[", "]"))
		<< "]";
	return out;
}

#endif /* __AUTODIFF_H */
