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

#if !defined(__TESTCASE_H)
#define __TESTCASE_H

#include <mitsuba/render/util.h>

MTS_NAMESPACE_BEGIN

#if defined(MTS_TESTCASE)
/**
 * When a testcase is being compiled, define the following preprocessor macros for convenience
 */
#define assertEquals(expected, actual) assertEqualsImpl(expected, actual, 0, __FILE__, __LINE__)
#define assertEqualsEpsilon(expected, actual, epsilon) assertEqualsImpl(expected, actual, epsilon, __FILE__, __LINE__)
#define assertTrue(expr) assertTrueImpl(expr, #expr, __FILE__, __LINE__)
#define assertFalse(expr) assertFalseImpl(expr, #expr, __FILE__, __LINE__)
#define failAndContinue(msg) failAndContinueImpl(msg, __FILE__, __LINE__)
#endif

/** \brief Base class of all testcases. Implementations of this
 * interface can be executed using the 'mtsutil' command. The execution
 * order is as follows: after initializaiton using init(), any tests
 * declared using the MTS_DECLARE_TEST() macro are executed. Finally,
 * the shutdown() method is called. See the files in 'mitsuba/src/tests'
 * for examples.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER TestCase : public Utility {
public:
	/**
	 * Perform any required initializations. The default
	 * implementation simply returns
	 */
	virtual void init();

	/**
	 * Execute any required shutdown code. The default
	 * implementation simply returns
	 */
	virtual void shutdown();

	/// Return the number of executed testcases
	inline int getExecuted() const { return m_executed; }
	
	/// Return the number of successfully executed testcases
	inline int getSucceeded() const { return m_succeeded; }

	MTS_DECLARE_CLASS()
protected:
	/// Virtual destructor
	virtual ~TestCase() { }

	/// Asserts that the two integer values are equal
	void assertEqualsImpl(int expected, int actual, const char *file, int line);

	/// Asserts that the two floating point values are equal
	void assertEqualsImpl(Float expected, Float actual, Float epsilon, const char *file, int line);
	
	/// Asserts that the two spectral power distributions are equal
	void assertEqualsImpl(const Spectrum &expected, const Spectrum &actual, Float epsilon, const char *file, int line);

	/// Asserts that the two 2D vectors are equal
	void assertEqualsImpl(const Vector2 &expected, const Vector2 &actual, Float epsilon, const char *file, int line);

	/// Asserts that the two 3D vectors are equal
	void assertEqualsImpl(const Vector &expected, const Vector &actual, Float epsilon, const char *file, int line);

	/// Asserts that the two 4D vectors are equal
	void assertEqualsImpl(const Vector4 &expected, const Vector4 &actual, Float epsilon, const char *file, int line);

	/// Asserts that the two 2D points are equal
	void assertEqualsImpl(const Point2 &expected, const Point2 &actual, Float epsilon, const char *file, int line);

	/// Asserts that the two 3D points are equal
	void assertEqualsImpl(const Point &expected, const Point &actual, Float epsilon, const char *file, int line);

	/// Asserts that the two 4x4 matrices are equal
	template<int M, int N> void assertEqualsImpl(const Matrix<M, N, Float> &expected, const Matrix<M, N, Float> &actual, Float epsilon, const char *file, int line) {
		bool match = true;
		for (int i=0; i<M; ++i)
			for (int j=0; j<N; ++j)
				if (std::abs(expected.m[i][j]-actual.m[i][j]) > epsilon)
					match = false;
		if (!match)
			Thread::getThread()->getLogger()->log(EError, NULL, file, line, "Assertion failure: "
				"expected matrix %s, got %s.", expected.toString().c_str(), actual.toString().c_str());
	}

	/// Asserts that a condition is true
	void assertTrueImpl(bool condition, const char *expr, const char *file, int line);

	/// Asserts that a condition is false
	void assertFalseImpl(bool condition, const char *expr, const char *file, int line);

	/// Note a failure and continue
	void failAndContinueImpl(const std::string &msg, const char *file, int line);

	/// Increase the number of succeeded tests
	void succeed();
protected:
	int m_executed, m_succeeded;
};

/**
 * The test supervisor is used when rendering a collection of testcases in the
 * form of scenes with analytic solutions. Reference output is compared against
 * the actual generated results and a user-specified type of test is executed
 * to decide between equality or inequality. Any problems are kept in a log, 
 * which can later be printed using <tt>printSummary()</tt>.
 * \ingroup librender
 */
class MTS_EXPORT_RENDER TestSupervisor : public Object {
public:
	/// Initialize a test supervisor for a specified number of testcases
	TestSupervisor(size_t total);

	/// Analyze the output of a rendered scene
	void analyze(const Scene *scene);

	/// Summarize the executed testcases
	void printSummary() const;

	MTS_DECLARE_CLASS()
protected:
	virtual ~TestSupervisor() { }
private:
	struct TestResult {
		bool success;
		fs::path input, output;
		std::string message;
	};

	size_t m_total, m_numFailed, m_numSucceeded;
	std::vector<TestResult> m_results;
	mutable ref<Mutex> m_mutex;
};

MTS_NAMESPACE_END

#define EXECUTE_GUARDED(name) \
	try { \
		Log(EInfo, "Executing test \"%s\" ..", #name); \
		m_executed++;\
		name();\
		m_succeeded++;\
	} catch (std::exception &e) {\
		Log(EInfo, "Testcase failed with error: %s", e.what());\
	}

#define MTS_BEGIN_TESTCASE() \
	MTS_DECLARE_CLASS() \
	int run(int argc, char **argv) {\
		init(); \
		Log(EInfo, "Executing testcase \"%s\" ..", getClass()->getName().c_str()); \
		m_executed = m_succeeded = 0;

#define MTS_DECLARE_TEST(name) \
		EXECUTE_GUARDED(name)

#define MTS_END_TESTCASE()\
		shutdown();\
		return 0;\
	}

#define MTS_EXPORT_TESTCASE(name, descr) \
	MTS_IMPLEMENT_CLASS(name, false, TestCase) \
	extern "C" { \
		void MTS_EXPORT *CreateUtility() { \
			return new name(); \
		} \
		const char MTS_EXPORT *GetDescription() { \
			return descr; \
		} \
	}

#endif /* __TESTCASE_H */
