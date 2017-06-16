/*
 * yasper - A non-intrusive reference counted pointer. 
 *	    Version: 1.02
 *			  
 *  Many ideas borrowed from Yonat Sharon and 
 *  Andrei Alexandrescu.
 *
 * (zlib license)
 * ----------------------------------------------------------------------------------	
 * Copyright (C) 2005 Alex Rubinsteyn
 *
 * This software is provided 'as-is', without any express or implied
 * warranty.  In no event will the authors be held liable for any damages
 * arising from the use of this software.
 *
 * Permission is granted to anyone to use this software for any purpose,
 * including commercial applications, and to alter it and redistribute it
 * freely, subject to the following restrictions:
 *
 * 1. The origin of this software must not be misrepresented; you must not
 *    claim that you wrote the original software. If you use this software
 *    in a product, an acknowledgment in the product documentation would be
 *    appreciated but is not required.
 * 2. Altered source versions must be plainly marked as such, and must not be
 *    misrepresented as being the original software.
 * 3. This notice may not be removed or altered from any source distribution.
 *
 * -----------------------------------------------------------------------------------
 * 
 * Send all questions, comments and bug reports to:
 * Alex Rubinsteyn (rubnstyn@uiuc.edu)
 * 
 * $Revision: 1956 $
 * $Date: 2007-01-12 11:33:49 -0800 (Fri, 12 Jan 2007) $
 * $Author: Alex Rubinsteyn (rubnstyn@uiuc.edu) $
 * $Id: yasper.hpp 1956 2007-01-12 19:33:49Z dbrown $
 * 
 */


#ifndef yasper_ptr_h
#define yasper_ptr_h

namespace yasper
{

#include <exception>
#include <string>

struct NullPointerException : public std::exception
{
    NullPointerException() throw() {}
    
    ~NullPointerException() throw() {}
    
    const char* what() const throw()
    {
        return "[Yasper Exception] Attempted to dereference null pointer";
    }
};

template <typename X> 
class ptr
{

public:
    typedef X element_type;

	/* 
		ptr needs to be its own friend so ptr< X > and ptr< Y > can access
		each other's private data members 
	*/ 
	template <class Y> friend class ptr; 
	/* 
		default constructor
			- don't create Counter
	*/
	ptr() : rawPtr(0), counter(0)
	{
	}
	
	/*
		Construct from a raw pointer 
	*/
	ptr(X* raw) : rawPtr(0), counter(0)
	{
		if (raw)
		{
			rawPtr = raw; 
			counter = new Counter;
		}
	}

	template <typename Y>
 	ptr(Y* raw) : rawPtr(0), counter(0) 
	{
		if (raw)
		{
			rawPtr = static_cast<X*>( raw );
			counter = new Counter; 
		}
	}
	
	/*
		Copy constructor 
	*/
	ptr(const ptr< X >& otherPtr)
	{
		acquire( otherPtr.counter );
		rawPtr = otherPtr.rawPtr;
	}
	
	template <typename Y>
	ptr(const ptr< Y >& otherPtr) : rawPtr(0), counter(0)
    	{
		acquire( (Counter*) (otherPtr.counter) );
		rawPtr = static_cast<X*>( otherPtr.GetRawPointer());
	}
	

	/* 
		Destructor 
	*/ 
	~ptr()
    	{
		release();
	}

/*
	Assignment to another ptr 
*/

ptr& operator=(const ptr< X >& otherPtr)
{
	if (this != &otherPtr)
	{
		release();
		acquire(otherPtr.counter); 
		rawPtr = otherPtr.rawPtr; 
	}
	return *this; 
}

template <typename Y>
ptr& operator=(const ptr< Y >& otherPtr)
{
	if ( this != (ptr< X >*) &otherPtr )
	{
		release();
		acquire( (Counter*)(otherPtr.counter) );
		rawPtr = static_cast<X*> (otherPtr.GetRawPointer()); 
	}
	return *this;
}

/*
	Assignment to raw pointers is really dangerous business.
	If the raw pointer is also being used elsewhere,
	we might prematurely delete it, causing much pain.
	Use sparingly/with caution.
*/

ptr& operator=(X* raw)
{
	if (raw)
	{
		release(); 
		counter = new Counter; 
		rawPtr = raw; 
	}
	return *this;
}

template <typename Y>
ptr& operator=(Y* raw)
{
	if (raw)
	{
		release();
		counter = new Counter; 
		rawPtr = static_cast<X*>(raw);
	}
	return *this;
}

/* 
	assignment to long to allow ptr< X > = NULL, 
	also allows raw pointer assignment by conversion. 
	Raw pointer assignment is really dangerous!
	If the raw pointer is being used elsewhere, 
	it will get deleted prematurely. 
*/ 
ptr& operator=(long num)
{
	if (num == 0)  //pointer set to null
	{
		release(); 
	}

	else //assign raw pointer by conversion
	{
		release();
		counter = new Counter; 
		rawPtr = reinterpret_cast<X*>(num);
	}	

	return *this; 
} 

/*
	Member Access
*/
	X* operator->() const 
	{
		return GetRawPointer(); 
	}


/*
	Dereference the pointer
*/
	X& operator* () const 
	{
		return *GetRawPointer(); 
	}


/*
	Conversion/casting operators
*/
	operator bool() const
	{
		return IsValid();
	}

	template <typename Y>
	operator Y*() const
	{
		return static_cast<Y*>(rawPtr);  
	}

	template <typename Y>
	operator const Y*() const
	{
		return static_cast<const Y*>(rawPtr);
	}


/*
	Provide access to the raw pointer 
*/

	X* GetRawPointer() const         
	{
		if (rawPtr == 0) throw new NullPointerException;
		return rawPtr;
	}

	
/* 
	Is there only one reference on the counter?
*/
	bool IsUnique() const
	{
		if (counter && counter->count == 1) return true; 
		else return false; 
	}
	
	bool IsValid() const
	{
		if (counter && rawPtr) return true;
		else return false; 
	}

	unsigned GetCount() const
	{
		if (counter) return counter->count;
		else return 0;
	}

private:
	X* rawPtr;
	
	struct Counter 
	{
		Counter(unsigned c = 1) : count(c) {}
		unsigned count; 
	};

	Counter* counter;

	// increment the count
	void acquire(Counter* c) 
	{ 
 		counter = c;
		if (c)
		{
			(c->count)++;
		}
	}

	// decrement the count, delete if it is 0
	void release()
	{ 
        if (counter) 
		{			
			(counter->count)--; 	

			if (counter->count == 0) 
			{
   				delete rawPtr;
				delete counter;			
			}
		}
		counter = 0;
		rawPtr = 0; 

	}
};


template <typename X, typename Y>
bool operator==(const ptr< X >& lptr, const ptr< Y >& rptr) 
{
	return lptr.GetRawPointer() == rptr.GetRawPointer(); 
}

template <typename X, typename Y>
bool operator==(const ptr< X >& lptr, Y* raw) 
{
	return lptr.GetRawPointer() == raw ; 
}

template <typename X>
bool operator==(const ptr< X >& lptr, long num)
{
	if (num == 0 && !lptr.IsValid())  //both pointer and address are null
	{
		return true; 
	}

	else //convert num to a pointer, compare addresses
	{
		return lptr == reinterpret_cast<X*>(num);
	}
	 
} 

template <typename X, typename Y>
bool operator!=(const ptr< X >& lptr, const ptr< Y >& rptr) 
{
	return ( !operator==(lptr, rptr) );
}

template <typename X, typename Y>
bool operator!=(const ptr< X >& lptr, Y* raw) 
{
		return ( !operator==(lptr, raw) );
}

template <typename X>
bool operator!=(const ptr< X >& lptr, long num)
{
	return (!operator==(lptr, num) ); 
}

template <typename X, typename Y>
bool operator&&(const ptr< X >& lptr, const ptr< Y >& rptr)
{
	return lptr.IsValid() &&  rptr.IsValid();
}

template <typename X>
bool operator&&(const ptr< X >& lptr, bool rval)
{
	return lptr.IsValid() && rval;
}

template <typename X>
bool operator&&(bool lval, const ptr< X >& rptr)
{
	return lval &&  rptr.IsValid();
}

template <typename X, typename Y>
bool operator||(const ptr< X >& lptr, const ptr< Y >& rptr)
{
	return lptr.IsValid() || rptr.IsValid();
}

template <typename X>
bool operator||(const ptr< X >& lptr, bool rval)
{
	return lptr.IsValid() || rval;
}

template <typename X>
bool operator||(bool lval, const ptr< X >& rptr)
{
	return lval || rptr.IsValid(); 
}

template <typename X>
bool operator!(const ptr< X >& p)
{
	return (!p.IsValid());
}


/* less than comparisons for storage in containers */
template <typename X, typename Y>
bool operator< (const ptr< X >& lptr, const ptr < Y >& rptr)
{
	return lptr.GetRawPointer() < rptr.GetRawPointer();
}

template <typename X, typename Y>
bool operator< (const ptr< X >& lptr, Y* raw)
{
	return lptr.GetRawPointer() < raw;
}

template <typename X, typename Y>
bool operator< (X* raw, const ptr< Y >& rptr)
{
	return raw < rptr.GetRawPointer();
}

} //close yasper namespace
#endif

