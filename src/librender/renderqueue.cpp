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

#include <mitsuba/core/timer.h>
#include <mitsuba/render/renderjob.h>

MTS_NAMESPACE_BEGIN

RenderQueue::RenderQueue() {
	m_mutex = new Mutex();
	m_joinMutex = new Mutex();
	m_cond = new ConditionVariable(m_mutex);
	m_timer = new Timer();
}
	
RenderQueue::~RenderQueue() {
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->decRef();
}

void RenderQueue::addJob(RenderJob *job) {
	m_mutex->lock();
	m_jobs[job] = JobRecord(m_timer->getMilliseconds());
	job->incRef();
	m_mutex->unlock();
}

void RenderQueue::registerListener(RenderListener *listener) {
	listener->incRef();
	m_mutex->lock();
	m_listeners.push_back(listener);
	m_mutex->unlock();
}

void RenderQueue::unregisterListener(RenderListener *listener) {
	m_mutex->lock();
	m_listeners.erase(std::remove(m_listeners.begin(), m_listeners.end(), listener));
	m_mutex->unlock();
	listener->decRef();
}
	
void RenderQueue::flush() {
	m_mutex->lock();
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.begin();
	for (; it != m_jobs.end(); ++it) {
		(*it).first->flush();
	}
	m_mutex->unlock();
}

void RenderQueue::removeJob(RenderJob *job, bool cancelled) {
	m_mutex->lock();
	std::map<RenderJob *, JobRecord>::iterator it = m_jobs.find(job);
	if (it == m_jobs.end()) {
		Log(EError, "RenderQueue::removeRenderJob() - job not found!");
		m_mutex->unlock();
	}
	JobRecord &rec = (*it).second;
	unsigned int ms = m_timer->getMilliseconds() - rec.startTime;
	Log(EInfo, "Render time: %s", timeString(ms/1000.0f, true).c_str());
	m_jobs.erase(job);
	m_cond->broadcast();
	m_joinMutex->lock();
	m_joinList.push_back(job);
	m_joinMutex->unlock();
	signalFinishJob(job, cancelled);
	m_mutex->unlock();
}
	
void RenderQueue::waitLeft(size_t njobs) const {
	m_mutex->lock();
	while (m_jobs.size() > njobs) 
		m_cond->wait();
	m_mutex->unlock();
	join();
}

void RenderQueue::join() const {
	m_joinMutex->lock();
	/* Wait for the proper termination of all stopping threads */
	for (size_t i=0; i<m_joinList.size(); ++i) {
		RenderJob *job = m_joinList[i];
		job->join();
		job->decRef();
	}
	m_joinList.clear();
	m_joinMutex->unlock();
}

void RenderQueue::signalWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workBeginEvent(job, wu, worker);
	m_mutex->unlock();
}

void RenderQueue::signalWorkEnd(const RenderJob *job, const ImageBlock *wr) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->workEndEvent(job, wr);
	m_mutex->unlock();
}

void RenderQueue::signalFinishJob(const RenderJob *job, bool cancelled) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->finishJobEvent(job, cancelled);
	m_mutex->unlock();
}

void RenderQueue::signalRefresh(const RenderJob *job, const Bitmap *bitmap) {
	m_mutex->lock();
	for (size_t i=0; i<m_listeners.size(); ++i)
		m_listeners[i]->refreshEvent(job, bitmap);
	m_mutex->unlock();
}

MTS_IMPLEMENT_CLASS(RenderQueue, false, Object)
MTS_IMPLEMENT_CLASS(RenderListener, true, Object)
MTS_NAMESPACE_END
