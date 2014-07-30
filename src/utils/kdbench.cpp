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

#include <mitsuba/render/util.h>
#include <mitsuba/core/timer.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/plugin.h>
#include <boost/algorithm/string.hpp>
#if defined(WIN32)
#include <mitsuba/core/getopt.h>
#endif

MTS_NAMESPACE_BEGIN

class KDBench : public Utility {
public:
	void help() {
		cout << endl;
		cout << "Synopsis: kd-tree performance benchmark. Traces uniformly distributed rays" << endl;
		cout << "though the bounding sphere of a scene and reports the resulting number of" << endl;
		cout << "rays per second. The main intent of this utility is to optimize the kd-tree" << endl;
		cout << "construction parameters for particular scenes and machines." << endl;
		cout << endl;
		cout << "Usage: mtsutil kdbench [options] <Scene XML file or PLY file>" << endl;
		cout << "Options/Arguments:" << endl;
		cout << "   -h             Display this help text" << endl << endl;
		cout << "   -t value       Specify the SAH traversal cost" << endl << endl;
		cout << "   -i value       Specify the SAH intersection cost" << endl << endl;
		cout << "   -e value       Specify the SAH empty space bonus" << endl << endl;
		cout << "   -b value       Specify the number of min-max bins" << endl << endl;
		cout << "   -c true/false  Enable/disable primitive clipping (aka. \"perfect splits\")" << endl << endl;
		cout << "   -p true/false  Enable/disable parallel tree construction" << endl << endl;
		cout << "   -r true/false  Enable/disable retraction of bad splits" << endl << endl;
		cout << "   -l value       Specify the primitive count, below which a leaf node" << endl;
		cout << "                  will always be created" << endl << endl;
		cout << "   -d depth       Specify the maximum tree depth" << endl << endl;
	 	cout << "   -x value       Specify the number of primitives, at which the " << endl;
		cout << "                  builder will switch from (approximate) Min-Max " << endl;
		cout << "                  binning to the more accurate O(n log n) SAH-based " << endl;
		cout << "                  optimization method." << endl << endl;
		cout << "   -f             Try to empirically find the best SAH cost values by" << endl;
		cout << "                  fitting the cost model to collected performance data" << endl << endl;
		cout << "Examples:" << endl;
		cout << "  E.g. to build a tree for the Stanford bunny having a low SAH cost, type " << endl << endl;
		cout << "  $ mtsutil kdbench -e .9 -l1 -d48 -x100000 data/tests/bunny.ply" << endl << endl;
		cout << "  To get SAH costs comparable to [Wald and Havran 06], also specify -t15 -i20" << endl << endl;
		cout << "  The high -x paramer effectively disables Min-Max binning, which " << endl;
		cout << "  leads to a slower and more memory-intensive build, so don't try" << endl;
		cout << "  this on a huge model." << endl << endl;
	}

	int run(int argc, char **argv) {
		ref<FileResolver> fileResolver = Thread::getThread()->getFileResolver();
		int optchar;
		char *end_ptr = NULL;
		Float intersectionCost = -1, traversalCost = -1, emptySpaceBonus = -1;
		int stopPrims = -1, maxDepth = -1, exactPrims = -1, minMaxBins = -1;
		bool clip = true, parallel = true, retract = true, fitParameters = false;
		optind = 1;

		/* Parse command-line arguments */
		while ((optchar = getopt(argc, argv, "i:t:e:c:p:r:l:x:b:d:hf")) != -1) {
			switch (optchar) {
				case 'h': {
						help();
						return 0;
					}
					break;
				case 'f':
					fitParameters = true;
					break;
				case 'i':
					intersectionCost = (Float) strtod(optarg, &end_ptr);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the intersection cost!");
					break;
				case 't':
					traversalCost = (Float) strtod(optarg, &end_ptr);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the traversal cost!");
					break;
				case 'e':
					emptySpaceBonus = (Float) strtod(optarg, &end_ptr);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the empty space bonus!");
					break;
				case 'l':
					stopPrims = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the stopping primitive count!");
					break;
				case 'd':
					maxDepth = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the maximum depth!");
					break;
				case 'x':
					exactPrims = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the -e parameter!");
					break;
				case 'b':
					minMaxBins = strtol(optarg, &end_ptr, 10);
					if (*end_ptr != '\0')
						SLog(EError, "Could not parse the min-max bins parameter!");
					break;
				case 'c':
					if (strcmp(optarg, "true") == 0)
						clip = true;
					else if (strcmp(optarg, "false") == 0)
						clip = false;
					else
						SLog(EError, "Could not parse the clipping parameter!");
					break;
				case 'p':
					if (strcmp(optarg, "true") == 0)
						parallel = true;
					else if (strcmp(optarg, "false") == 0)
						parallel = false;
					else
						SLog(EError, "Could not parse the parallel build parameter!");
					break;
				case 'r':
					if (strcmp(optarg, "true") == 0)
						retract = true;
					else if (strcmp(optarg, "false") == 0)
						retract = false;
					else
						SLog(EError, "Could not parse the retraction parameter!");
					break;
			};
		}

		if (optind == argc || optind+1 < argc) {
			help();
			return 0;
		}

		ref<Scene> scene;
		ref<ShapeKDTree> kdtree;

		std::string lowercase = boost::to_lower_copy(std::string(argv[optind]));
		if (boost::ends_with(lowercase, ".xml")) {
			fs::path
				filename = fileResolver->resolve(argv[optind]),
				filePath = fs::absolute(filename).parent_path(),
				baseName = filename.stem();
			ref<FileResolver> frClone = fileResolver->clone();
			frClone->prependPath(filePath);
			Thread::getThread()->setFileResolver(frClone);
			scene = loadScene(argv[optind]);
			kdtree = scene->getKDTree();
		} else if (boost::ends_with(lowercase, ".ply")) {
			Properties props("ply");
			props.setString("filename", argv[optind]);
			ref<TriMesh> mesh;
			mesh = static_cast<TriMesh *> (PluginManager::getInstance()->
					createObject(MTS_CLASS(TriMesh), props));
			mesh->configure();
			kdtree = new ShapeKDTree();
			kdtree->addShape(mesh);
		} else {
			Log(EError, "The supplied scene filename must end in either PLY or XML!");
		}

		if (intersectionCost != -1)
			kdtree->setQueryCost(intersectionCost);
		if (traversalCost != -1)
			kdtree->setTraversalCost(traversalCost);
		if (emptySpaceBonus != -1)
			kdtree->setEmptySpaceBonus(emptySpaceBonus);
		if (stopPrims != -1)
			kdtree->setStopPrims(stopPrims);
		if (maxDepth != -1)
			kdtree->setMaxDepth(maxDepth);
		if (exactPrims != -1)
			kdtree->setExactPrimitiveThreshold(exactPrims);
		if (minMaxBins != -1)
			kdtree->setMinMaxBins(minMaxBins);
		kdtree->setClip(clip);
		kdtree->setRetract(retract);
		kdtree->setParallelBuild(parallel);

		/* Show some statistics, and make sure it roughly fits in 80cols */
		Logger *logger = Thread::getThread()->getLogger();
		DefaultFormatter *formatter = ((DefaultFormatter *) logger->getFormatter());
		logger->setLogLevel(EDebug);
		formatter->setHaveDate(false);

		if (scene)
			scene->initialize();
		else
			kdtree->build();

		BSphere bsphere(kdtree->getAABB().getBSphere());
		const size_t nRays = 5000000;

		if (!fitParameters) {
			Log(EInfo, "Bounding sphere: %s", bsphere.toString().c_str());
			Float best = 0;
			for (int j=0; j<3; ++j) {
				ref<Random> random = new Random();
				ref<Timer> timer = new Timer();
				size_t nIntersections = 0;

				Log(EInfo, "Shooting " SIZE_T_FMT " rays (1 thread, incoherent) ..", nRays);

				for (size_t i=0; i<nRays; ++i) {
					Point2 sample1(random->nextFloat(), random->nextFloat()),
						sample2(random->nextFloat(), random->nextFloat());
					Point p1 = bsphere.center + warp::squareToUniformSphere(sample1) * bsphere.radius;
					Point p2 = bsphere.center + warp::squareToUniformSphere(sample2) * bsphere.radius;
					Ray r(p1, normalize(p2-p1), 0.0f);

					Intersection its;
					if (kdtree->rayIntersect(r, its))
						nIntersections++;
				}

				Log(EInfo, "Found " SIZE_T_FMT " intersections in %i ms",
					nIntersections, timer->getMilliseconds());
				Float mrays = nRays / (timer->getMilliseconds() * (Float) 1000);
				Log(EInfo, "-> %.3f MRays/s", mrays);
				Log(EInfo, "");
				best = std::max(best, mrays);
			}
			Log(EInfo, "Best of three: %.3f MRays/s", best);
		} else {
			Float intersectionCost, traversalCost;
			kdtree->findCosts(intersectionCost, traversalCost);
		}

		Thread::getThread()->getLogger()->setLogLevel(EInfo);
		return 0;
	}

	MTS_DECLARE_UTILITY()
};

MTS_EXPORT_UTILITY(KDBench, "kd-tree performance benchmark")
MTS_NAMESPACE_END
