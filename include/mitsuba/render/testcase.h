#if !defined(__TESTCASE_H)
#define __TESTCASE_H

#include <mitsuba/render/util.h>

MTS_NAMESPACE_BEGIN

/** \brief Base class of all testcases. Implementations of this
 * interface can be executed using the 'mtsutil' command. The execution
 * order is as follows: after initializaiton using init(), any tests
 * declared using the MTS_DECLARE_TEST() macro are executed. Finally,
 * the shutdown() method is called. See the files in 'mitsuba/src/tests'
 * for examples.
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
protected:
	int m_executed, m_succeeded;
};

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
