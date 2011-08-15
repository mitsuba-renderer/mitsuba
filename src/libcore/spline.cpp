/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2011 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/core/spline.h>

MTS_NAMESPACE_BEGIN

CubicSpline::CubicSpline(Stream *stream, InstanceManager *manager) {
	size_t size = stream->readSize();
	m_x.resize(size); m_y.resize(size); m_deriv.resize(size);
	stream->readFloatArray(&m_x[0], size);
	stream->readFloatArray(&m_y[0], size);
	stream->readFloatArray(&m_deriv[0], size);
}

void CubicSpline::build() {
	if (m_x.size() != m_y.size())
		Log(EError, "build(): encountered a different number "
			" of X and Y components!");
	size_t size = m_x.size();
	if (size < 3) 
		Log(EError, "build(): need at least three points!");

	Float *temp = new Float[size];
	m_deriv.resize(size);

	/* Solve a tridiagonal system (based on Numerical Recipes) */
	m_deriv[0] = temp[0] = 0.0f;
	for(size_t i=1; i<size-1; i++){
		Float sig = (m_x[i] - m_x[i-1]) / (m_x[i+1] - m_x[i-1]);
		Float invP = 1.0f / (sig * m_deriv[i-1] + 2);
		m_deriv[i]=(sig-1) * invP;
		temp[i]=(m_y[i+1]-m_y[i])  / (m_x[i+1]-m_x[i]) 
		       - (m_y[i]-m_y[i-1]) / (m_x[i]-m_x[i-1]);
		temp[i]=(6*temp[i] / (m_x[i+1]-m_x[i-1]) - sig*temp[i-1]) * invP;
	}
	m_deriv[size-1] = 0.0f;

	/* Backsubstitute */
	for(ptrdiff_t k=size-2; k>=0; k--)
		m_deriv[k] = m_deriv[k] * m_deriv[k+1] + temp[k];

	delete[] temp;
}

void CubicSpline::clear() {
	m_x.clear();
	m_y.clear();
	m_deriv.clear();
}

Float CubicSpline::eval(Float x) const {
	size_t lo = 0, hi = m_x.size()-1;

	/* Binary search */
	while (hi-lo > 1){
		size_t k = (hi+lo) >> 1;
		if (m_x[k] > x)
			hi = k;
		else
			lo = k;
	}

	const Float h = m_x[hi]-m_x[lo],
		        invH = 1.0f / h,
		        a = (m_x[hi]-x) * invH,
				b = (x- m_x[lo]) * invH;

	return a*m_y[lo] + b*m_y[hi] + (1.0f/6.0f) *
		((a*a*a-a)*m_deriv[lo] + (b*b*b-b)*m_deriv[hi]) * (h*h);
}

void CubicSpline::serialize(Stream *stream, InstanceManager *manager) const {
	Assert(m_x.size() == m_y.size() && m_x.size() == m_deriv.size());
	stream->writeSize(m_x.size());
	stream->writeFloatArray(&m_x[0], m_x.size());
	stream->writeFloatArray(&m_y[0], m_x.size());
	stream->writeFloatArray(&m_deriv[0], m_x.size());
}

std::string CubicSpline::toString() const {
	std::ostringstream oss;
	Assert(m_x.size() == m_y.size());
	oss << "CubicSpline[" << endl
		<< "  nodeCount = " << m_x.size() << "," << endl
		<< "  nodes = {" << endl;
	for (size_t i=0; i<m_x.size(); ++i) {
		oss << "    " << m_x[i] << " => " << m_y[i];
		if (i+1 < m_x.size())
			oss << ",";
		oss << endl;
	}
	oss << "  }" << endl
		<< "]";
	return oss.str();
}

MTS_IMPLEMENT_CLASS_S(CubicSpline, false, SerializableObject)
MTS_NAMESPACE_END
