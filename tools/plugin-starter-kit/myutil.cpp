#include <mitsuba/render/util.h>

MTS_NAMESPACE_BEGIN

class MyUtility : public Utility {
public:
	MyUtility(UtilityServices *us) : Utility(us) { }

	int run(int argc, char **argv) {
		cout << "Hello world!" << endl;
		return 0;
	}

	MTS_DECLARE_CLASS()
};

MTS_IMPLEMENT_CLASS(MyUtility, false, Utility)
MTS_EXPORT_UTILITY(MyUtility, "Example utility")
MTS_NAMESPACE_END
