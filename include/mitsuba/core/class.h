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

#if !defined(__CLASS_H)
#define __CLASS_H

MTS_NAMESPACE_BEGIN

/* Forward declarations */
class Stream;
class Object;
class InstanceManager;

/** \brief Universal class descriptor.
 * @see ref, Object
 */
class MTS_EXPORT_CORE Class {
public:
	/// Construct a new class descriptor
	Class(const std::string &name, bool abstract, const std::string &superClassName, 
		void *instPtr = NULL, void *unSerPtr = NULL);

	/// Return the class' name
	inline const std::string &getName() const { return m_name; }

	/** \brief Return whether the class represented
	 * by this Class object is abstract
	 */
	inline bool isAbstract() const { return m_abstract; }

	/** \brief Does the class support instantiation over RTTI?
	 */
	inline bool isInstantiable() const { return m_instPtr != NULL; }

	/** \brief Does the class support serialization?
	 */
	inline bool isSerializable() const { return m_unSerPtr != NULL; }

	/** \brief Return this class' super class or NULL
	 * if it has no super class
	 */
	inline const Class *getSuperClass() const { return m_superClass; }

	/// Check whether this class derives from pClass
	bool derivesFrom(const Class *pClass) const;

	/// Look up a class by its name
	static const Class *forName(const std::string &name);

	/** \brief Look up a class by its name. Avoids allocating
	 * heap space by taking a character array as parameter
	 */
	static const Class *forName(const char *name);

	/** \brief Unserialize an instance of the class (if this is
	 * supported by the class).
	 */
	Object *unserialize(Stream *stream = NULL, InstanceManager *manager = NULL) const;

	/// Generate an instance of this class (if supported)
	Object *instantiate() const;

	/** \brief Initializes the built-in RTTI and creates
	 * a list of all compiled classes
	 */
	static void staticInitialization();

	/// Free the memory taken by staticInitialization()
	static void staticShutdown();
private:
	/** \brief Initialize a class - called by
	 * staticInitialization()
	 */
	static void initializeOnce(Class *theClass);
private:
	std::string m_name;
	bool m_abstract;
	Class *m_superClass;
	std::string m_superClassName;
	void *m_instPtr, *m_unSerPtr;
};

/* Class instantiation macros */
#define MTS_DECLARE_CLASS() \
	virtual const Class *getClass() const; \
public: \
	static Class *m_theClass; 

// basic RTTI support for a class
#define MTS_IMPLEMENT_CLASS(name, abstract, super) \
	Class *name::m_theClass = new Class(#name, abstract, #super); \
	const Class *name::getClass() const { \
		return m_theClass;\
	}

// Extended version, records that the class supports instantiation by name
#define MTS_IMPLEMENT_CLASS_I(name, abstract, super) \
	Object *__##name ##_inst() { \
		return new name(); \
	} \
	Class *name::m_theClass = new Class(#name, abstract, #super, (void *) &__##name ##_inst, NULL); \
	const Class *name::getClass() const { \
		return m_theClass;\
	}

// Extended version, records that the class can be unserialized from a binary data stream
#define MTS_IMPLEMENT_CLASS_S(name, abstract, super) \
	Object *__##name ##_unSer(Stream *stream, InstanceManager *manager) { \
		return new name(stream, manager); \
	} \
	Class *name::m_theClass = new Class(#name, abstract, #super, NULL, (void *) &__##name ##_unSer); \
	const Class *name::getClass() const { \
		return m_theClass;\
	}

// Extended version, records that the class can be unserialized from a binary data stream as well as instantiated by name
#define MTS_IMPLEMENT_CLASS_IS(name, abstract, super) \
	Object *__##name ##_unSer(Stream *stream, InstanceManager *manager) { \
		return new name(stream, manager); \
	} \
	Object *__##name ##_inst() { \
		return new name(); \
	} \
	Class *name::m_theClass = new Class(#name, abstract, #super, (void *) &__##name ##_inst, (void *) &__##name ##_unSer); \
	const Class *name::getClass() const { \
		return m_theClass;\
	}

MTS_NAMESPACE_END

#endif /* __CLASS_H */

