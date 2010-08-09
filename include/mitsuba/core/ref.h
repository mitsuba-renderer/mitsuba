#if !defined(__REFERENCE_H)
#define __REFERENCE_H

MTS_NAMESPACE_BEGIN

/** \brief A simple wrapper class which takes care of
 * referencing and unreferencing objects
 */
template <typename T> class ref {
public:
	/// Create a NULL reference
	ref() : m_ptr(NULL) { }

	/// Construct a reference from a pointer
	ref(T *ptr) : m_ptr(ptr) { if (m_ptr) m_ptr->incRef(); }
	
	/// Copy-constructor
	ref(const ref &pRef) : m_ptr(pRef.m_ptr) { if (m_ptr) m_ptr->incRef(); }

	/// Destroy this reference
	~ref() { if (m_ptr) m_ptr->decRef(); }

	/// Overwrite this reference with another reference
	inline ref& operator= (const ref& pref) {
		if (m_ptr == pref.m_ptr)
			return *this;
		T* tmp = m_ptr;
		m_ptr = pref.m_ptr;
		if (m_ptr)
			m_ptr->incRef();
		if (tmp)
			tmp->decRef();
		return *this;
	}
	
	/// Overwrite this reference with a pointer to another object
	inline ref& operator= (T *ptr) {
		if (m_ptr == ptr)
			return *this;
		T* tmp = m_ptr;
		m_ptr = ptr;
		if (m_ptr)
			m_ptr->incRef();
		if (tmp)
			tmp->decRef();
		return *this;
	}

	/// Compare this reference with another reference
	inline bool operator== (const ref &pref) const { return (m_ptr == pref.m_ptr); }
	
	/// Compare this reference with another reference
	inline bool operator!= (const ref &pref) const { return (m_ptr != pref.m_ptr); }
	
	/// Compare this reference with a pointer
	inline bool operator== (const T* ptr) const { return (m_ptr == ptr); }
	
	/// Compare this reference with a pointer
	inline bool operator!= (const T* ptr) const { return (m_ptr != ptr); }

	/// Check whether this is a NULL reference
	inline bool operator!() const { return (m_ptr == NULL); }

	/// Access the object referenced by this reference
	inline T* operator-> () { return m_ptr; }
	
	/// Access the object referenced by this reference
	inline const T* operator-> () const { return m_ptr; }

	/// Return a C++ reference to the referenced object
	inline T& operator*() { return *m_ptr; }
	
	/// Return a C++ reference to the referenced object
	inline const T& operator*() const { return *m_ptr; }

	/// Return a pointer to the referenced object
	inline operator T* () { return m_ptr; }

	/// Return a pointer to the referenced object
	inline T* get() { return m_ptr; }
	
	/// Return a pointer to the referenced object
	inline const T* get() const { return m_ptr; }

	/**
	 * Return a string representation of this reference
	 */
	inline std::string toString() const {
		return formatString("ref<%s>[ref=%i, ptr=%s]",
				m_ptr == NULL ? "null" : m_ptr->getClass()->getName().c_str(),
				m_ptr == NULL ? -1 : m_ptr->getRefCount(),
				m_ptr == NULL ? "null" : m_ptr->toString().c_str());
	}
private:
	T *m_ptr;
};

MTS_NAMESPACE_END

#endif /* __REFERENCE_H */
