#include <mitsuba/render/imageproc_wu.h>
#include "preview_proc.h"

PreviewProcess::PreviewProcess(const Scene *scene, int sceneResID, int blockSize) 
	: m_vpl(NULL) {
	m_blockSize = blockSize;
	m_logLevel = ETrace;
	m_mutex = new Mutex();
	m_scene = scene;
	m_film = m_scene->getFilm();
	bindResource("scene", sceneResID);
}

bool PreviewProcess::isLocal() const {
	return true;
}

void PreviewProcess::processResult(const WorkResult *result, bool cancelled) {
	const ImageBlock *block = static_cast<const ImageBlock *>(result);
	const int sx = block->getOffset().x - m_film->getCropOffset().x, 
		sy = block->getOffset().y - m_film->getCropOffset().y;
	const int ex = sx + block->getSize().x, ey = sy + block->getSize().y;
	int pos = 0;
	Float r=0, g=0, b=0;

	if (m_source) {
		for (int y=sy; y<ey; ++y) {
			const float *source = m_source->getFloatData() + ((m_target->getHeight() - 1 - y) * m_source->getWidth() + sx) * 3;
			float *target = m_target->getFloatData() + ((m_target->getHeight() - 1 - y) 
				* m_target->getWidth() + sx) * 3;
			for (int x=sx; x<ex; ++x) {
				block->getPixel(pos++).toLinearRGB(r, g, b);
				*target++ = (float) (r + *source++); 
				*target++ = (float) (g + *source++); 
				*target++ = (float) (b + *source++); 
			}
		}
	} else {
		for (int y=sy; y<ey; ++y) {
			float *target = m_target->getFloatData() + ((m_target->getHeight() - 1 - y) 
				* m_target->getWidth() + sx) * 3;
			for (int x=sx; x<ex; ++x) {
				block->getPixel(pos++).toLinearRGB(r, g, b);
				*target++ = (float) r;
				*target++ = (float) g;
				*target++ = (float) b;
			}
		}
	}

	m_mutex->lock();
	m_numRays += block->getExtra();
	m_mutex->unlock();
}

void PreviewProcess::configure(const VPL &vpl, Float minDist, const Point2 &jitter, 
		const Bitmap *source, Bitmap *target, bool coherent) {
	BlockedImageProcess::init(m_film->getCropOffset(), m_film->getCropSize(), m_blockSize);
	m_source = source;
	m_target = target;
	m_numRays = 0;
	m_vpl = &vpl;
	m_minDist = minDist;
	m_coherent = coherent;

	/* It is not necessary to shoot normalized rays. Instead, interpolate: 
	   here, we generate the upper left corner ray as well as the 
	   per-pixel increments */
	Ray topLeftRay, rightRay, bottomRay;
	const Point2 topLeft(jitter);
	const Point2 right(topLeft.x + m_film->getSize().x, topLeft.y);
	const Point2 bottom(topLeft.x, topLeft.y + m_film->getSize().y);
	const Point2 lens(0, 0);

	const Camera *camera = m_scene->getCamera();
	camera->generateRay(topLeft, lens, topLeftRay);
	camera->generateRay(right, lens, rightRay);
	camera->generateRay(bottom, lens, bottomRay);
	m_cameraTL = Vector(topLeftRay.d);
	m_cameraO = camera->getPosition();
	m_cameraDx = (rightRay.d - topLeftRay.d)
		/ (Float) m_film->getSize().x;
	m_cameraDy = (bottomRay.d - topLeftRay.d)
		/ (Float) m_film->getSize().y;
}

ref<WorkProcessor> PreviewProcess::createWorkProcessor() const {
	return new PreviewWorker(m_blockSize, m_cameraO, m_cameraTL,
		m_cameraDx, m_cameraDy, *m_vpl, m_minDist, m_coherent);
}


MTS_IMPLEMENT_CLASS(PreviewProcess, false, BlockedImageProcess)
