/*
Copyright (c) 2026, Michael Kazhdan
All rights reserved.

Redistribution and use in source and binary forms, with or without modification,
are permitted provided that the following conditions are met:

Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer. Redistributions in binary form must reproduce
the above copyright notice, this list of conditions and the following disclaimer
in the documentation and/or other materials provided with the distribution. 

Neither the name of the Johns Hopkins University nor the names of its contributors
may be used to endorse or promote products derived from this software without specific
prior written permission. 

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO THE IMPLIED WARRANTIES 
OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
TO, PROCUREMENT OF SUBSTITUTE  GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
DAMAGE.
*/

#pragma once


#include <array>

#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Misha/MultiThreading.h"
#include "Misha/Atomic.h"
#include "Misha/Atomic.Geometry.h"
#include "Misha/OrthogonalFrame.h"

namespace MishaK::ECD
{
	template< unsigned int HalfRes , typename ... Data >
	struct Descriptor : public std::pair< double , std::array< std::tuple< double , Data ... > , ( 2 * HalfRes + 1 ) * ( 2 * HalfRes + 1 ) > >
	{
		using Base = std::pair< double , std::array< std::tuple< double , Data ... > , ( 2 * HalfRes + 1 ) * ( 2 * HalfRes + 1 ) > >;

		static const unsigned int Width = 2 * HalfRes + 1;
		static const unsigned int Dim = Width * Width;

		using IndexedDescriptor = std::pair< Descriptor , MishaK::template RegularGrid< 2 >::Index >;

		// Default (zero) initialization
		Descriptor( void );

		// Clears the descriptor values
		Descriptor & clear( void );

		// Returns the magnitude of the descriptor
		double magnitude( void ) const;

		// Returns the normalized descriptor
		Descriptor normalized( void ) const;

		// Returns the distance between descriptors
		static double Difference( const Descriptor & desc1 , const Descriptor & desc2 );

		// Returns the descriptor for every pixel
		static MishaK::RegularGrid< 2 , Descriptor > Get( const MishaK::RegularGrid< 2 , MishaK::Point< double , 2 > > & grads , const MishaK::RegularGrid< 2 , Data > & ... values , double radius , bool gather=true );

	protected:
		template< bool ThreadSafe >
		void _splat( double x , double y , double weight , Data ... data );

		template< unsigned int Idx , bool ThreadSafe >
		void _add( size_t idx , std::tuple< double , Data... > entry );

		template< unsigned int Idx >
		void _scale( double scl );

		template< unsigned int Idx >
		static double _Difference( const Descriptor & desc1 , const Descriptor & desc2 );

		template< unsigned int Idx >
		static void _CheckDataGridResolutions( std::tuple< const MishaK::RegularGrid< 2 , Data > & ... > dataGrids , unsigned int width , unsigned int height );
	};

	// Returns the top <num> descriptors (i.e. with largest magnitude)
	template< unsigned int HalfRes , typename ... Data >
	std::vector< Descriptor< HalfRes , Data ... > > GetTop( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int num );

	// Returns the top <num> descriptors with their associated pixel indices
	template< unsigned int HalfRes , typename ... Data >
	std::vector< typename Descriptor< HalfRes , Data ... >::IndexedDescriptor > GetTopIndexed( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int num );

	// Returns the top <num> descriptors (i.e. with largest magnitude)
	template< unsigned int HalfRes , typename ... Data >
	std::vector< Descriptor< HalfRes , Data ... > > GetTopStratified( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int stratificationRes , unsigned int num );

	// Returns the top <num> descriptors with their associated pixel indices
	template< unsigned int HalfRes , typename ... Data >
	std::vector< typename Descriptor< HalfRes , Data ... >::IndexedDescriptor > GetTopStratifiedIndexed( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int stratificationRes , unsigned int num );

#include "ECD.inl"
}
