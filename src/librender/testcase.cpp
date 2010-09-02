/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2010 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/scene.h>
#include <mitsuba/render/testcase.h>
#include <boost/math/distributions/students_t.hpp>
#include <fstream>

MTS_NAMESPACE_BEGIN

void TestCase::init() { }
void TestCase::shutdown() { }

struct Sample {
	Float value;
	Float variance;
	int nSamples;
};

static std::vector<Float> parseRefFile(std::ifstream &is) {
	std::string line;
	std::vector<Float> result;
	while (!is.eof() && !is.fail()) {
		std::getline(is, line);
		std::vector<std::string> tokens = tokenize(line, " \t;,[]");

		for (size_t i=0; i<tokens.size(); ++i) {
			char *end_ptr = NULL;
			Float val = (Float) strtod(tokens[i].c_str(), &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Error while parsing a testcase output file");
			result.push_back(val);
		}
	}
	return result;
}

static std::vector<Sample> parseMFile(std::ifstream &is, int testType) {
	std::string line;
	std::vector<Sample> result;
	while (!is.eof() && !is.fail()) {
		std::getline(is, line);
		std::vector<std::string> tokens = tokenize(line, " \t;,[]");
		SAssert(testType == Scene::ERelativeError || (tokens.size() % 3) == 0);

		for (size_t i=0; i<tokens.size(); ) {
			Sample sample;
			char *end_ptr = NULL;
			sample.value = (Float) strtod(tokens[i++].c_str(), &end_ptr);
			if (*end_ptr != '\0')
				SLog(EError, "Error while parsing a testcase output file");
			if (testType == Scene::ETTest) {
				sample.variance = (Float) strtod(tokens[i++].c_str(), &end_ptr);
				if (*end_ptr != '\0')
					SLog(EError, "Error while parsing a testcase output file");
				sample.nSamples = (int) strtol(tokens[i++].c_str(), &end_ptr, 10);
				if (*end_ptr != '\0')
					SLog(EError, "Error while parsing a testcase output file");
				}
			result.push_back(sample);
		}
	}
	return result;
}

TestSupervisor::TestSupervisor(size_t total) 
 : m_total(total), m_numFailed(0), m_numSucceeded(0) {
	m_mutex = new Mutex();
}

void TestSupervisor::analyze(const Scene *scene) {
	using namespace boost::math;

	TestResult result;

	result.input = scene->getSourceFile();
	result.output = scene->getDestinationFile() + ".m";
	result.success = false;
	std::string refFilename = scene->getDestinationFile() + ".ref";

	std::ifstream is(result.output.c_str());
	std::ifstream is_ref(refFilename.c_str());
	if (is.fail()) {
		result.message = formatString("Could not open '%s'!", result.output.c_str());
		m_mutex->lock();
		m_numFailed++; m_results.push_back(result);
		m_mutex->unlock();
		return;
	}
	if (is_ref.fail()) {
		result.message = formatString("Could not open '%s'!", refFilename.c_str());
		m_mutex->lock();
		m_numFailed++; m_results.push_back(result);
		m_mutex->unlock();
		return;
	}

	std::vector<Sample> actual = parseMFile(is, scene->getTestType());
	std::vector<Float> ref = parseRefFile(is_ref);

	is.close();
	is_ref.close();

	if (actual.size() != ref.size()) {
		result.message = formatString("Output format does not match the reference (%i vs %i pixels)!", 
			(int) actual.size(), (int) ref.size());
		m_mutex->lock();
		m_numFailed++; m_results.push_back(result);
		m_mutex->unlock();
		return;
	}

	if (scene->getTestType() == Scene::ETTest) {
		for (size_t i=0; i<actual.size(); ++i) {
			int df = actual[i].nSamples-1;
			Float var = std::max(actual[i].variance, Epsilon);
			Float T = (actual[i].value - ref[i]) * std::sqrt(actual[i].nSamples / var);
			students_t dist(df);
			Float pval = (Float) (2*cdf(complement(dist, std::abs(T))));
			Log(EDebug, "Performing a t-test: result=%f (ref=%f), diff=%e, var=%f, T-stat=%f, df=%i, p-value=%f", actual[i].value, ref[i], actual[i].value - ref[i], var, T, df, pval);
			if (pval <= scene->getTestThreshold()) {
				Log(EWarn, "t-test REJECTS!");
				result.message = formatString("t-test REJECTS: result=%f (ref=%f), diff=%e, var=%f T-stat=%f, df=%i, p-value=%f", actual[i].value, ref[i], actual[i].value - ref[i], var, T, df, pval);
				m_mutex->lock();
				m_numFailed++; m_results.push_back(result);
				m_mutex->unlock();
				return;
			} else {
				Log(EDebug, "t-test accepts.");
			}
		}
	} else if (scene->getTestType() == Scene::ERelativeError) {
		for (size_t i=0; i<actual.size(); ++i) {
			Float relerr = std::abs((actual[i].value - ref[i]) / std::max(Epsilon, ref[i]));
			Log(EDebug, "Testing the relative error: result=%f (ref=%f), diff=%e, relerr=%f", actual[i].value, ref[i], actual[i].value - ref[i], relerr);
			if (relerr > scene->getTestThreshold()) {
				Log(EWarn, "Relativ error threshold EXCEEDED!");
				result.message = formatString("Relative error threshold EXCEEDED: result=%f (ref=%f), diff=%e, relerr=%f", actual[i].value, ref[i], actual[i].value - ref[i], relerr);
				m_mutex->lock();
				m_numFailed++; m_results.push_back(result);
				m_mutex->unlock();
				return;
			} else {
				Log(EDebug, "Relative error accepted.");
			}
		}
	} else if (scene->getTestType() == Scene::ENone) {
		Log(EError, "No test type specified, don't know what to do!");
	} else {
		Log(EError, "Unknown test type!");
	}

	result.success = true;
	m_mutex->lock();
	m_numSucceeded++;
	m_results.push_back(result);
	m_mutex->unlock();
}

void TestSupervisor::printSummary() const {
	m_mutex->lock();
	Log(EInfo, "Ran %i/%i testcases, %i succeeded, %i failed.", (int) (m_numFailed+m_numSucceeded), (int) m_total,
			(int) m_numSucceeded, (int) m_numFailed);

	for (size_t i=0; i<m_results.size(); ++i) {
		const TestResult &result = m_results[i];
		if (result.success)
			continue;
		Log(EWarn, "============================================================");
		Log(EWarn, " Failure: Test case %zi (\"%s\")", i+1, result.input.c_str());
		Log(EWarn, " Message: \"%s\"", result.message.c_str());
		Log(EWarn, "============================================================");
	}
	m_mutex->unlock();
}

MTS_IMPLEMENT_CLASS(TestSupervisor, false, Object)
MTS_IMPLEMENT_CLASS(TestCase, false, Utility)
MTS_NAMESPACE_END
