#include <mitsuba/render/common.h>

MTS_NAMESPACE_BEGIN

std::ostream &operator<<(std::ostream &os, const ETransportQuantity &quantity) {
	switch (quantity) {
		case EImportance: os << "importance"; break;
		case ERadiance:   os << "radiance"; break;
		default: os << "invalid"; break;
	};
	return os;
}

MTS_NAMESPACE_END
