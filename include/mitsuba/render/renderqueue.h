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

#if !defined(__RENDERQUEUE_H)
#define __RENDERQUEUE_H

#include <mitsuba/mitsuba.h>
#include <mitsuba/render/imageproc_wu.h>

MTS_NAMESPACE_BEGIN

class RenderJob;

/** 
 * Abstract render listener - can be used to react to 
 * progress messages (e.g. in a GUI)
 */
class MTS_EXPORT_RENDER RenderListener : public Object {
public:
	/// Called when work has begun in a rectangular image region
	virtual void workBeginEvent(const RenderJob *job, const RectangularWorkUnit *wu, int worker) = 0;
	
	/// Called when work has finished in a rectangular image region
	virtual void workEndEvent(const RenderJob *job, const ImageBlock *wr) = 0;
	
	/// Called when the whole target image has been altered in some way
	virtual void refreshEvent(const RenderJob *job) = 0;
	
	/// Called when a render job has completed successfully or unsuccessfully
	virtual void finishJobEvent(const RenderJob *job, bool cancelled) = 0;

	MTS_DECLARE_CLASS()
protected:
	virtual ~RenderListener() { }
};

/** 
 * Render queue - used to keep track of a number of scenes
 * that are simultaneously being rendered. Also distributes
 * events regarding these scenes to registered listeners.
 */
class MTS_EXPORT_RENDER RenderQueue : public Object {
public:
	/// Create a new render queue
	RenderQueue();

	/// Return the current number of jobs in the queue
	inline int getJobCount() const { return m_jobs.size(); }

	/// Add a render job to the queue
	void addJob(RenderJob *thr);

	/// Remove a (finished) render job from the queue
	void removeJob(RenderJob *thr, bool wasCancelled);

	/// Register a render listener
	void registerListener(RenderListener *listener);

	/// Unregister a render listener
	void unregisterListener(RenderListener *listener);

	/** 
	 * Wait until the queue contains a certain number
	 * of scenes (or less).
	 */
	void waitLeft(size_t njobs) const;

	/// Releases resources held by recently finished jobs
	void join() const;

	/// Cause all render jobs to write out the current image
	void flush();

	/* Event distribution */
	void signalWorkBegin(const RenderJob *job, const RectangularWorkUnit *wu, int worker);
	void signalWorkEnd(const RenderJob *job, const ImageBlock *block);
	void signalFinishJob(const RenderJob *job, bool cancelled);
	void signalRefresh(const RenderJob *job);

	MTS_DECLARE_CLASS()
private:
	/// Virtual destructor
	virtual ~RenderQueue();
protected:
	struct JobRecord {
		/* Only starting time for now */
		unsigned int startTime;

		inline JobRecord() { }
		inline JobRecord(unsigned int startTime)
			: startTime(startTime) {
		}
	};

	std::map<RenderJob *, JobRecord> m_jobs;
	mutable std::vector<RenderJob *> m_joinList;
	mutable ref<Mutex> m_mutex, m_joinMutex;
	mutable ref<ConditionVariable> m_cond;
	ref<Timer> m_timer;
	std::vector<RenderListener *> m_listeners;
};

MTS_NAMESPACE_END

#endif /* __RENDERQUEUE_H */
