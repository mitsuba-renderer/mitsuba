#if !defined(__TLS_H)
#define __TLS_H

#include <mitsuba/mitsuba.h>
#include <pthread.h>

/**
 * Thin wrapper around posix thread local storage
 */
template <typename T> class ThreadLocal {
public:
	inline ThreadLocal() {
		pthread_key_create(&m_key, ThreadLocal::deletePtr);
	}

	inline ~ThreadLocal() {
		pthread_key_delete(m_key);
	}

	inline void set(T *value) {
		T *previous = get();
		if (previous != NULL)
			previous->decRef();
		if (value != NULL)
			value->incRef();
		pthread_setspecific(m_key, value);
	}

	inline T *get() {
		return static_cast<T *>(pthread_getspecific(m_key));
	}
protected:
	static void deletePtr(void *ptr) {
		static_cast<T *>(ptr)->decRef();
	}
protected:
	pthread_key_t m_key;
};

/**
 * As above, but meant for stack-allocated objects
 */
template <typename T> class PrimitiveThreadLocal {
public:
	inline PrimitiveThreadLocal() {
		pthread_key_create(&m_key, PrimitiveThreadLocal::deletePtr);
	}

	/**
	 * Beware: all threads that ever used this TLS (other than
	 * the caller) must either either died or called cleanup()
	 * before this destructor is called -- otherwise memory leaks
	 * will ensue.
	 */
	inline ~PrimitiveThreadLocal() {
		cleanup();
		pthread_key_delete(m_key);
	}

	inline void set(T &value) {
		get() = value;
	}

	inline T &get() {
		T *value = static_cast<T *>(pthread_getspecific(m_key));
		if (value == NULL) {
			value = new T();
			pthread_setspecific(m_key, value);
		}

		return *value;
	}

	inline const T &get() const {
		T *value = static_cast<T *>(pthread_getspecific(m_key));
		if (value == NULL) {
			value = new T();
			pthread_setspecific(m_key, value);
		}

		return *value;
	}

	void cleanup() {
		T *value = static_cast<T *>(pthread_getspecific(m_key));
		if (value != NULL)
			delete value;
		pthread_setspecific(m_key, NULL);
	}
protected:
	static void deletePtr(void *ptr) {
		delete static_cast<T *>(ptr);
	}
protected:
	pthread_key_t m_key;
};


#endif /* __TLS_H */
