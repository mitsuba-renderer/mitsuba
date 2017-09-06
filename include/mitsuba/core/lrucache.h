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

#if !defined(__LRUCACHE_H)
#define __LRUCACHE_H

#include <mitsuba/mitsuba.h>
#include <boost/bimap.hpp>
#include <boost/bimap/list_of.hpp>
#include <boost/bimap/set_of.hpp>
#include <boost/function.hpp>

MTS_NAMESPACE_BEGIN

/**
 * \brief Generic LRU cache implementation
 *
 * Based on the bimap implementation by Tim Day
 * (http://timday.bitbucket.org/lru.html).
 *
 * This cache does not support multithreading out of the box -- it
 * will need to be protected using some form of locking mechanism.
 *
 * The original code is under the following license:
 *
 * <pre>
 * Copyright (c) 2010, Tim Day <timday@timday.com>
 * Permission to use, copy, modify, and/or distribute this software for any
 * purpose with or without fee is hereby granted, provided that the above
 * copyright notice and this permission notice appear in all copies.
 *
 * THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
 * WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
 * MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
 * ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
 * WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
 * ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
 * OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
 * </pre>
 *
 * \tparam K Key data type
 * \tparam KComp Key comparator
 * \tparam V Value data type
 * \ingroup libcore
 */
template <typename K, typename KComp, typename V> struct LRUCache : public Object {
public:
    typedef int dummy_type;

    // Bimap with key access on left view, key access
    // history on right view, and associated value.
    typedef boost::bimaps::bimap<
            boost::bimaps::set_of<K, KComp>,
            boost::bimaps::list_of<dummy_type>,
            boost::bimaps::with_info<V> > cache_type;

    LRUCache() { }

    // Constuctor specifies the cached function and
    // the maximum number of records to be stored.
    LRUCache(size_t capacity,
        const boost::function<V(const K&)>& generatorFunction,
        const boost::function<void (const V&)>& cleanupFunction = NULL)
        : m_capacity(capacity), m_generatorFunction(generatorFunction),
          m_cleanupFunction(cleanupFunction) {
        SAssert(m_capacity != 0);
    }

    virtual ~LRUCache() {
        typename cache_type::right_iterator
            src = m_cache.right.begin();
        if (m_cleanupFunction) {
            while (src != m_cache.right.end())
                m_cleanupFunction((*src++).info);
        }
    }

    bool isFull() const {
        return m_cache.size() == m_capacity;
    }

    // Obtain value of the cached function for k
    V get(const K& k, bool &hit) {
        // Attempt to find existing record
        const typename cache_type::left_iterator it
            = m_cache.left.find(k);

        if (it == m_cache.left.end()) {
            // We don't have it:
            // Evaluate function and create new record

            const V v = m_generatorFunction(k);
            insert(k,v);
            hit = false;
            return v;
        } else {
            // We do have it:
            // Update the access record view.

            m_cache.right.relocate(
                m_cache.right.end(),
                m_cache.project_right(it)
            );
            hit = true;
        }
        return it->info;
    }

    // Obtain the cached keys, most recently used element
    // at head, least recently used at tail.
    // This method is provided purely to support testing.
    template <typename IT> void get_keys(IT dst) const {
        typename cache_type::right_const_reverse_iterator
            src = m_cache.right.rbegin();
        while (src != m_cache.right.rend())
            *dst++=(*src++).second;
    }
protected:
    void insert(const K& k,const V& v) {
        SAssert(m_cache.size() <= m_capacity);
        if (m_cache.size() == m_capacity) {
            if (m_cleanupFunction)
                m_cleanupFunction(m_cache.right.begin()->info);
            // If necessary, make space
            // by purging the least-recently-used element
            m_cache.right.erase(m_cache.right.begin());
        }

        // Create a new record from the key, a dummy and the value
        m_cache.insert(typename cache_type::value_type(k,0,v));
    }

private:
    size_t m_capacity;
    boost::function<V(const K&)> m_generatorFunction;
    boost::function<void(const V&)> m_cleanupFunction;
    cache_type m_cache;
};

MTS_NAMESPACE_END

#endif /* __LRUCACHE_H */
