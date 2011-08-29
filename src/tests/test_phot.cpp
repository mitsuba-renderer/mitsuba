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

#include <mitsuba/core/plugin.h>
#include <mitsuba/core/kdtree.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/render/testcase.h>
#include <mitsuba/render/gatherproc.h>

MTS_NAMESPACE_BEGIN

class TestPhotonMap : public TestCase {
public:
	MTS_BEGIN_TESTCASE()
	MTS_DECLARE_TEST(test01_photonmap)
	MTS_END_TESTCASE()

	void test01_photonmap() {
		ref<FileResolver> fres = Thread::getThread()->getFileResolver();
		fres->addPath("scenes/sponza2");
		ref<Scene> scene = loadScene("scenes/sponza2/sponza.xml");
		scene->configure();
		scene->initialize();

		size_t nPhotons = 1000000;

		ref<Scheduler> sched = Scheduler::getInstance();
		ref<GatherPhotonProcess> proc = new GatherPhotonProcess(
			GatherPhotonProcess::ESurfacePhotons, nPhotons,
			1000, 10, 10, false, true, NULL);

		ref<Sampler> sampler = static_cast<Sampler *>(
			PluginManager::getInstance()->createObject(Properties("halton")));

		int sceneResID = sched->registerResource(scene);
		int cameraResID = sched->registerResource(scene->getCamera());
		int qmcSamplerID = sched->registerResource(sampler);
		proc->bindResource("scene", sceneResID);
		proc->bindResource("camera", cameraResID);
		proc->bindResource("sampler", qmcSamplerID);

		sched->schedule(proc);
		sched->wait(proc);
			
		ref<PhotonMap> pmap = proc->getPhotonMap();
		pmap->setScaleFactor(1 / (Float) proc->getShotParticles());

		pmap->build();

		ref<Random> random = new Random();
		PhotonMap::SearchResult results[101];
		size_t nQueries = 100000;
		Point *queryPositions = new Point[nQueries];
		for (size_t i=0; i<nQueries; ++i)
			queryPositions[i] = (*pmap)[random->nextUInt(pmap->size())].getPosition();

		ref<Timer> timer = new Timer();
		Log(EInfo, "Testing query performance (photon map) ..");
		size_t nResults1 = 0;
		for (size_t i=0; i<nQueries; ++i) {
			Float searchRadius = std::numeric_limits<Float>::infinity();
			nResults1 += pmap->nnSearch(queryPositions[i], searchRadius, 100, results);
		}
		Log(EInfo, "Done. Took %i ms, got " SIZE_T_FMT " results",  timer->getMilliseconds(), nResults1);
		Log(EInfo, "Photon size: %i bytes", (int) sizeof(Photon));
	}
};

MTS_EXPORT_TESTCASE(TestPhotonMap, "Testcase for kd-tree related code")
MTS_NAMESPACE_END
