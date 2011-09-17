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

#include <mitsuba/render/util.h>
#include <mitsuba/core/quad.h>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/bind.hpp>
#include <Eigen/SVD>
#include <iomanip>
#include "/Users/wenzel/mitsuba/lookup.h"

MTS_NAMESPACE_BEGIN

class UFlakeFit : public Utility {
public:
	struct Reference {
		Reference(size_t res, size_t nSines) : m_res(res), m_sines(nSines) {
			m_directions = new Vector[m_res];
			m_result = new Float[m_res];
			m_error = new Float[m_res];
			for (size_t i=0; i<m_res; ++i) {
				Float theta = i * M_PI / ((Float) m_res - 1);
				m_directions[i] = sphericalDirection(theta, 0);
				m_result[i] = m_error[i] = 0;
			}
		}

		~Reference() {
			delete[] m_directions;
			delete[] m_result;
			delete[] m_error;
		}

		void f(size_t nPts, const Float *in, Float *out) {
			Vector wf(0.0f, 0.0f, 1.0f);
			#pragma omp parallel for
			for (int i=0; i<(int) nPts; ++i) {
				Vector d = sphericalDirection(in[2*i], in[2*i+1]);
				Float sinTheta = std::sin(in[2*i]);
				Float D = sinTheta * std::fastexp(-std::pow(d.z, 2)/(2*m_gamma*m_gamma)) / 
					(std::pow(2*M_PI, (Float) 3 / (Float) 2) * m_gamma * erf(1/(std::sqrt(2)*m_gamma)));
				for (size_t j=0; j<m_res; ++j)
					out[j*nPts+i] = D * absDot(d, m_directions[j]);
			}
		}

		void run() {
			std::ofstream os("reference.m");
			std::ofstream os2("coeffs.m");
			std::ofstream os3("coeffs.h");
			os << std::fixed << std::setprecision(16) << endl;
			os2 << std::fixed << std::setprecision(16) << endl;
			os3 << std::scientific << std::setprecision(16) << endl;
			size_t nElements = 1000;
			os << "ref = {" << endl;
			os2 << "coeffs = {" << endl;
			os3 << "const double coeffs[" << nElements << "][" << m_sines << "] = {" << endl;
			Float range = 4;
			Float factor = std::pow(range, (Float) 1 / (Float) 4)/nElements;
			Float globalAbsErr = 0, globalRelErr = 0, globalAvgAbsErr = 0, globalAvgRelErr = 0;
			Float refAvgAbsErr = 0, refAvgRelErr = 0;

			for (size_t j=1; j<=nElements; ++j) {
				Float tmp = factor * j, tmp2 = tmp*tmp;
				m_gamma = tmp2*tmp2;
				size_t attempt = 0, evals = 0;
				Float min[2] = { 0, 0 } , max[2] = { M_PI, 2*M_PI };
				ref<Timer> timer = new Timer();
				Float base = refAvgAbsErr;
				do {
					min[0] = M_PI/2 * (1.0 - std::pow((Float) 2, - (Float) attempt)),
					max[0] = M_PI/2 * (1.0 + std::pow((Float) 2, - (Float) attempt));
					NDIntegrator quad(m_res, 2, 200000, 0, 1e-8f);
					quad.integrateVectorized(boost::bind(
						&Reference::f, this, _1, _2, _3), min, max, m_result, m_error, evals);
					Float maxError = 0;
					for (size_t i=0; i<m_res; ++i)
						maxError = std::max(m_error[i], maxError);
					Log(EInfo, "gamma=%e (attempt %i): used " SIZE_T_FMT " function evaluations, error = %f, took %i ms",
							m_gamma, attempt+1, evals, maxError, timer->getMilliseconds());
					++attempt;
				} while (evals < 20000);

				Eigen::MatrixXd A(m_res, m_sines);
				Eigen::VectorXd b(m_res);
				for (size_t i=0; i<m_res; ++i) {
					Float theta = i * M_PI / ((Float) m_res - 1),
						  sinTheta = std::sin(theta),
						  value = 1;
					for (size_t k=0; k<m_sines; ++k) {
						A(i, k) = value;
						value *= sinTheta;
					}
					b(i) = m_result[i];
					
					refAvgAbsErr += std::abs(m_result[i]-sigmaT(m_gamma, sinTheta)) / (Float) b.rows();
					refAvgRelErr += std::abs(m_result[i]-sigmaT(m_gamma, sinTheta)) / (m_result[i] * b.rows());
				}
				Eigen::VectorXd coeffs = A.jacobiSvd(Eigen::ComputeThinU
					| Eigen::ComputeThinV).solve(b);

				Float absErr = (A*coeffs-b).array().abs().maxCoeff();
				Float avgAbsErr = (A*coeffs-b).array().abs().sum() / b.rows();
				Log(EInfo, "Average abs. error = %f, reference = %f", avgAbsErr, refAvgAbsErr-base);

				Float relErr = ((A*coeffs-b).array() / b.array()).abs().maxCoeff();
				Float avgRelErr = ((A*coeffs-b).array() / b.array()).abs().sum() / b.rows();
				globalAbsErr = std::max(globalAbsErr, absErr);
				globalRelErr = std::max(globalRelErr, relErr);
				globalAvgAbsErr += avgAbsErr;
				globalAvgRelErr += avgRelErr;

				os2 << "\t{ ";
				os3 << "\t{ ";
				for (size_t i=0; i<m_sines; ++i) {
					if (i != 0) {
						os2.width(24);
						os3.width(24);
					}
					os2 << coeffs[i];
					os3 << coeffs[i];
					if (i+1 < m_sines) {
						os2 << ", ";
						os3 << ", ";
					}
				}

				os << "\t{ ";
				for (size_t i=0; i<m_res; ++i) {
					os << m_result[i];
					if (i+1 < m_res)
						os << ", ";
				}
				os << " }";
				os2 << " }";
				os3 << " }";
				if (j+1 <= nElements) {
					os << "," << endl;
					os2 << "," << endl;
					os3 << "," << endl;
				}
			}
			os << endl << "};" << endl;
			os2 << endl << "};" << endl;
			os3 << endl << "};" << endl;
			os.close();
			os2.close();
			os3.close();
			globalAvgAbsErr /= nElements;
			globalAvgRelErr /= nElements;
			refAvgAbsErr /= nElements;
			refAvgRelErr /= nElements;
			cout << "Maximum errors:" << endl;
			cout << "  Global max. abs. err = " << globalAbsErr << endl;
			cout << "  Global max. rel. err = " << globalRelErr << endl << endl;
			cout << "Average errors:" << endl;
			cout << "  Global avg. abs. err = " << globalAvgAbsErr << endl;
			cout << "  Reference avg. abs. err = " << refAvgAbsErr << endl;
			cout << "  Global avg. rel. err = " << globalAvgRelErr << endl;
			cout << "  Reference avg. rel. err = " << refAvgRelErr << endl;
		}
	private:
		size_t m_res, m_sines;
		Vector *m_directions;
		Float *m_result, *m_error;
		Float m_gamma;
	};

	int run(int argc, char **argv) {
		for (int i=10; i<=10; i+=2) {
			Reference ref(200, i);
			ref.run();
			cout << "(that was with " << i << " terms)" << endl;
			cout << endl;
			cout << endl;
		}
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(UFlakeFit, "Precompute data for the microflake model")
MTS_NAMESPACE_END
