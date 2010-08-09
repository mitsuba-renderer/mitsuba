/*
 * Licensed to the Apache Software Foundation (ASF) under one or more
 * contributor license agreements.  See the NOTICE file distributed with
 * this work for additional information regarding copyright ownership.
 * The ASF licenses this file to You under the Apache License, Version 2.0
 * (the "License"); you may not use this file except in compliance with
 * the License.  You may obtain a copy of the License at
 *
 *      http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */

/*
 * $Id: CMStateSet.hpp 677430 2008-07-16 21:05:31Z borisk $
 */

#if !defined(XERCESC_INCLUDE_GUARD_CMSTATESET_HPP)
#define XERCESC_INCLUDE_GUARD_CMSTATESET_HPP

//  DESCRIPTION:
//
//  This class is a specialized bitset class for the content model code of
//  the validator. It assumes that its never called with two objects of
//  different bit counts, and that bit sets smaller than 64 bits are far
//  and away the most common. So it can be a lot more optimized than a general
//  purpose utility bitset class
//

#include <xercesc/util/ArrayIndexOutOfBoundsException.hpp>
#include <xercesc/util/RuntimeException.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/MemoryManager.hpp>
#include <string.h>

XERCES_CPP_NAMESPACE_BEGIN

class CMStateSetEnumerator;

class CMStateSet : public XMemory
{
public :
    // -----------------------------------------------------------------------
    //  Constructors and Destructor
    // -----------------------------------------------------------------------
    CMStateSet( const unsigned int bitCount
              , MemoryManager* const manager = XMLPlatformUtils::fgMemoryManager) :

        fBitCount(bitCount)
        , fBitArray(0)
        , fMemoryManager(manager)
    {
        //
        //  See if we need to allocate the byte array or whether we can live
        //  within the 64 bit high performance scheme.
        //
        if (fBitCount > 64)
        {
            fArraySize = fBitCount / 32;
            if (fBitCount % 32)
                fArraySize++;
            fBitArray = (XMLInt32*) fMemoryManager->allocate(fArraySize*sizeof(XMLInt32));
        }
        else
        {
            fArraySize = 2;
            fBitArray = fBits;
        }

        // Init all the bits to zero
        zeroBits();
    }


    /*
     * This method with the 'for' statement (commented out) cannot be made inline
     * because the antiquated CC (CFront) compiler under HPUX 10.20 does not allow
     * the 'for' statement inside any inline method. Unfortunately,
     * we have to support it. So instead, we use memcpy().
     */

    CMStateSet(const CMStateSet& toCopy) :
        XMemory(toCopy)
      , fBitCount(toCopy.fBitCount)
      , fBitArray(0)
      , fMemoryManager(toCopy.fMemoryManager)
    {
        //
        //  See if we need to allocate the byte array or whether we can live
        //  within the 64 bit high performance scheme.
        //
        if (fBitCount > 64)
        {
            fArraySize = fBitCount / 32;
            if (fBitCount % 32)
                fArraySize++;
            fBitArray = (XMLInt32*) fMemoryManager->allocate(fArraySize*sizeof(XMLInt32));
        }
        else
        {
            fArraySize = 2;
            fBitArray = fBits;
        }


        memcpy((void *) fBitArray,
               (const void *) toCopy.fBitArray,
               fArraySize * sizeof(XMLInt32));

        // for (unsigned int index = 0; index < fArraySize; index++)
        //     fBitArray[index] = toCopy.fBitArray[index];
    }

    ~CMStateSet()
    {
        if (fBitArray!=fBits)
            fMemoryManager->deallocate(fBitArray);
    }


    // -----------------------------------------------------------------------
    //  Set manipulation methods
    // -----------------------------------------------------------------------
    void operator|=(const CMStateSet& setToOr)
    {
        for (unsigned int index = 0; index < fArraySize; index++)
            fBitArray[index] |= setToOr.fBitArray[index];
    }

    bool operator==(const CMStateSet& setToCompare) const
    {
        if (fBitCount != setToCompare.fBitCount)
            return false;

        for (unsigned int index = 0; index < fArraySize; index++)
        {
            if (fBitArray[index] != setToCompare.fBitArray[index])
                return false;
        }
        return true;
    }

    CMStateSet& operator=(const CMStateSet& srcSet)
    {
        if (this == &srcSet)
            return *this;

        // They have to be the same size
        if (fBitCount != srcSet.fBitCount)
            ThrowXMLwithMemMgr(RuntimeException, XMLExcepts::Bitset_NotEqualSize, fMemoryManager);

        for (unsigned int index = 0; index < fArraySize; index++)
            fBitArray[index] = srcSet.fBitArray[index];

        return *this;
    }


    bool getBit(const unsigned int bitToGet) const
    {
        if (bitToGet >= fBitCount)
            ThrowXMLwithMemMgr(ArrayIndexOutOfBoundsException, XMLExcepts::Bitset_BadIndex, fMemoryManager);

        const XMLInt32 mask = (0x1UL << (bitToGet % 32));
        const unsigned int byteOfs = bitToGet / 32;
        // And access the right bit and byte
        return (fBitArray[byteOfs]!=0 && (fBitArray[byteOfs] & mask) != 0);
    }

    bool isEmpty() const
    {
        for (unsigned int index = 0; index < fArraySize; index++)
        {
            if (fBitArray[index] != 0)
                return false;
        }
        return true;
    }

    void setBit(const unsigned int bitToSet)
    {
        if (bitToSet >= fBitCount)
            ThrowXMLwithMemMgr(ArrayIndexOutOfBoundsException, XMLExcepts::Bitset_BadIndex, fMemoryManager);

        const XMLInt32 mask = (0x1UL << (bitToSet % 32));
        const unsigned int byteOfs = bitToSet / 32;

        // And access the right bit and byte
        fBitArray[byteOfs] &= ~mask;
        fBitArray[byteOfs] |= mask;
    }

    void zeroBits()
    {
        for (unsigned int index = 0; index < fArraySize; index++)
            fBitArray[index] = 0;
    }

    XMLSize_t hashCode() const
    {
        XMLSize_t hash = 0;
        for (XMLSize_t index = 0; index<fArraySize; index++)
            hash = fBitArray[index] + hash * 31;
        return hash;
    }

private :
    // -----------------------------------------------------------------------
    //  Unimplemented constructors and operators
    // -----------------------------------------------------------------------
    CMStateSet();


    // -----------------------------------------------------------------------
    //  Private data members
    //
    //  fBitCount
    //      The count of bits that the outside world wants to support,
    //      so its the max bit index plus one.
    //
    //  fArraySize
    //      If the bit count is > 64, then we use the fBitArray member to
    //      store the bits, and this indicates its size in bytes. Otherwise
    //      its value is meaningless and unset.
    //
    //  fBits
    //      When the bit count is <= 64 (very common), these hold the bits.
    //      Otherwise, the fBitArray member holds htem.
    //
    //  fBitArray
    //      The array of bytes used when the bit count is > 64. It is
    //      allocated as required.
    // -----------------------------------------------------------------------
    unsigned int    fBitCount;
    unsigned int    fArraySize;
    XMLInt32        fBits[2];
    XMLInt32*       fBitArray;
    MemoryManager*  fMemoryManager;

    friend class CMStateSetEnumerator ;
};

class CMStateSetEnumerator : public XMemory
{
public:
    CMStateSetEnumerator(const CMStateSet* toEnum) :
      fToEnum(toEnum),
      fIndexCount((unsigned int)-1),
      fLastValue(0),
      fByteArrayCursor(0)
    {
        findNext();
    }

    bool hasMoreElements()
    {
        return fLastValue!=0;
    }

    unsigned int nextElement()
    {
        for(int i=0;i<32;i++)
        {
            XMLInt32 mask=(1UL << i);
            if(fLastValue & mask)
            {
                fLastValue &= ~mask;
                unsigned int retVal=fIndexCount+i;
                if(fLastValue==0)
                    findNext();
                return retVal;
            }
        }
        return 0;
    }

private:
    void findNext()
    {
        unsigned int nOffset=((fIndexCount==(unsigned int)-1)?0:(fIndexCount/32)+1), i;
        for(i=nOffset;i<fToEnum->fArraySize;i++)
        {
            if(fToEnum->fBitArray[i]!=0)
            {
                fIndexCount=i*32;
                fLastValue=fToEnum->fBitArray[i];
                break;
            }
        }
    }

    const CMStateSet*                   fToEnum;
    unsigned int                        fIndexCount;
    XMLInt32                            fLastValue;
    unsigned int                        fByteArrayCursor;
};

XERCES_CPP_NAMESPACE_END

#endif
