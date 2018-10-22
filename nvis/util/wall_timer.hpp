/*************************************************************************
* nvis: "no name visualization": suite of handy data structures and
*       agorithmic building blocks for scientific visualization.
*
* Author: Christoph Garth, University of Kaiserslautern
*
* Copyright (C) 2006 Christoph Garth
*
*************************************************************************/


#ifndef __timer_hpp
#define __timer_hpp

#ifdef _WIN32
	#include <windows.h>
#else
	#include <sys/time.h>
#endif
#include <ctime>

namespace nvis
{

#ifdef _WIN32
    // Windows: --- a wall clock timer modeled after boost::timer ---

    class timer
    {
	private:
		LARGE_INTEGER frequency;
		LARGE_INTEGER t1;
		
    public:

	timer() 
	{
		//get ticks per second
		QueryPerformanceFrequency (&frequency);
		QueryPerformanceCounter(&t1);
	}
	
	void restart() 
	{ 
	    QueryPerformanceCounter(&t1);
	}
	
	double elapsed() const
	{ 
		LARGE_INTEGER t2;
		QueryPerformanceCounter(&t2);
		double elapsedTime = (t2.QuadPart - t1.QuadPart) * 1000.0 /frequency.QuadPart;
	    return elapsedTime;
	}
 	
	};

#else
/* a wall clock timer with microsecond precision,
 * modeled after boost::timer
 */

class timer
{
public:
    
    timer() 
    {
	gettimeofday( &_start_time, 0 );
    }
    
    void restart() 
    { 
	gettimeofday( &_start_time, 0 );
    }
    
    double elapsed() const
    { 
	timeval cur;
	gettimeofday( &cur, 0 );
	
	return timeval_diff( _start_time, cur );
    }
    
private:
    
    static double timeval_diff( const timeval& t1,
				const timeval& t2 )
    {
	timeval r;
	
	r.tv_sec  = t2.tv_sec  - t1.tv_sec;
	r.tv_usec = t2.tv_usec - t1.tv_usec;
	
	if( r.tv_usec < 0 )
	{
	    --r.tv_sec;
	    r.tv_usec += 1000000;
	}
	
	return (double)r.tv_sec + (double)r.tv_usec / 1000000.0;
    }
    
    timeval _start_time;
};

#endif //Linux v. Windows

} // namespace nvis

#endif //__timer_hpp
