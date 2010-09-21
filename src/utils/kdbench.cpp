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

#include <mitsuba/render/util.h>
#include <mitsuba/core/timer.h>

MTS_NAMESPACE_BEGIN

class KDBench : public Utility {
public:
	int run(int argc, char **argv) {
		if (argc != 2) {
			cout << "Test the kd-Tree performance using a specified scene" << endl;
			cout << "Syntax: mtsutil kdbench <scene.xml>" << endl;
			return -1;
		}

		ref<Scene> scene = loadScene(argv[1]);
		Thread::getThread()->getLogger()->setLogLevel(EDebug);
		scene->initialize();
		BSphere bsphere(scene->getBSphere());

		for (int j=0; j<3; ++j) {
			ref<Random> random = new Random();
			ref<Timer> timer = new Timer();
			size_t nRays = 5000000, nIntersections = 0;

			Log(EInfo, "Bounding sphere: %s", bsphere.toString().c_str());
			Log(EInfo, "Shooting " SIZE_T_FMT " rays ..", nRays);

			for (size_t i=0; i<nRays; ++i) {
				Point2 sample1(random->nextFloat(), random->nextFloat()),
					sample2(random->nextFloat(), random->nextFloat());
				Point p1 = bsphere.center + squareToSphere(sample1) * bsphere.radius;
				Point p2 = bsphere.center + squareToSphere(sample2) * bsphere.radius;
				Ray r(p1, normalize(p2-p1));

				Intersection its;
				if (scene->rayIntersect(r, its))
					nIntersections++;
			}

			Log(EInfo, "Found " SIZE_T_FMT " intersections in %i ms",
				nIntersections, timer->getMilliseconds());
			Log(EInfo, "%.3f MRays/s", 
				nRays / (timer->getMilliseconds() * (Float) 1000));
		}

		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(KDBench, "kd-tree performance benchmark")
MTS_NAMESPACE_END
