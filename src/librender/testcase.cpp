/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

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

#include <mitsuba/render/scene.h>
#include <mitsuba/render/testcase.h>
#include <boost/math/distributions/students_t.hpp>
#include <boost/filesystem/fstream.hpp>

MTS_NAMESPACE_BEGIN

void TestCase::init() { }
void TestCase::shutdown() { }

void TestCase::assertTrueImpl(bool value, const char *expr, const char *file, int line) {
	if (!value)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion '%s == true' failed!", expr);
}

void TestCase::assertFalseImpl(bool value, const char *expr, const char *file, int line) {
	if (value)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion '%s == false' failed!", expr);
}

void TestCase::assertEqualsImpl(int actual, int expected, Float epsilon, const char *file, int line) {
	if (std::abs(actual-expected)>epsilon)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected integer value %i, got %i.", expected, actual);
}

void TestCase::assertEqualsImpl(Float actual, Float expected, Float epsilon, const char *file, int line) {
	if (std::abs(actual-expected) > epsilon)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected floating point value %f, got %f.", expected, actual);
}

void TestCase::assertEqualsImpl(const Spectrum &actual, const Spectrum &expected, Float epsilon, const char *file, int line) {
	bool match = true;
	for (int i=0; i<SPECTRUM_SAMPLES; ++i)
		if (std::abs(actual[i]-expected[i]) > epsilon)
			match = false;
	if (!match)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected vector %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
}

void TestCase::assertEqualsImpl(const Vector2 &actual, const Vector2 &expected, Float epsilon, const char *file, int line) {
	bool match = true;
	for (int i=0; i<2; ++i)
		if (std::abs(actual[i]-expected[i]) > epsilon)
			match = false;
	if (!match)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected vector %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
}


void TestCase::assertEqualsImpl(const Point2 &actual, const Point2 &expected, Float epsilon, const char *file, int line) {
	bool match = true;
	for (int i=0; i<2; ++i)
		if (std::abs(actual[i]-expected[i]) > epsilon)
			match = false;
	if (!match)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected point %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
}

void TestCase::assertEqualsImpl(const Vector &actual, const Vector &expected, Float epsilon, const char *file, int line) {
	bool match = true;
	for (int i=0; i<3; ++i)
		if (std::abs(actual[i]-expected[i]) > epsilon)
			match = false;
	if (!match)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected vector %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
}

void TestCase::assertEqualsImpl(const Point &actual, const Point &expected, Float epsilon, const char *file, int line) {
	bool match = true;
	for (int i=0; i<3; ++i)
		if (std::abs(actual[i]-expected[i]) > epsilon)
			match = false;
	if (!match)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected point %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
}

void TestCase::assertEqualsImpl(const Vector4 &actual, const Vector4 &expected, Float epsilon, const char *file, int line) {
	bool match = true;
	for (int i=0; i<4; ++i)
		if (std::abs(actual[i]-expected[i]) > epsilon)
			match = false;
	if (!match)
		Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
			"expected vector %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
}

void TestCase::failAndContinueImpl(const std::string &msg, const char *file, int line) {
	Thread::getThread()->getLogger()->log(EWarn, NULL, file, line, "Failure: %s", msg.c_str());
	m_executed++;
}

void TestCase::succeed() {
	m_executed++;
	m_succeeded++;
}

MTS_IMPLEMENT_CLASS(TestCase, false, Utility)
MTS_NAMESPACE_END
