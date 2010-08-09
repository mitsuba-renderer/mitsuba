#if !defined(__TESTCASE_H)
#define __TESTCASE_H

#include <mitsuba/mitsuba.h>

MTS_NAMESPACE_BEGIN

/**
 * The test supervisor is used when rendering a collection of testcases in the
 * form of scenes with analytic solutions. Reference output is compared against
 * the actual generated results and a user-specified type of test is executed
 * to decide between equality or inequality. Any problems are kept in a log, 
 * which can later be printed using <tt>printSummary()</tt>.
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
		std::string input, output;
		std::string message;
	};

	size_t m_total, m_numFailed, m_numSucceeded;
	std::vector<TestResult> m_results;
	mutable ref<Mutex> m_mutex;
};

MTS_NAMESPACE_END

#endif /* __TESTCASE_H */
