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

#include "Misha/Geometry.h"
#include "Misha/RegularGrid.h"
#include "Misha/MultiThreading.h"


namespace MishaK::LocalFit
{

	template< typename WeightFunction >
	concept HasWeightFunction = requires( WeightFunction f , unsigned int x , unsigned int y ) { { f(x,y) } -> std::same_as< double >; };

	template< HasWeightFunction WeightFunction >
	MishaK::RegularGrid< 2 , MishaK::Point< double , 3 > > Linear( const MishaK::RegularGrid< 2 , double > & scalars , const WeightFunction & weights , double dev , double radius );

	template< HasWeightFunction WeightFunction >
	MishaK::RegularGrid< 2 , MishaK::Point< double , 2 > > Gradients( const MishaK::RegularGrid< 2 , double > & scalars , const WeightFunction & weights , double dev , double radius );

#include "LocalFit.inl"
};
