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

#include <stdio.h>
#include <stdlib.h>
#include <vector>

#include "Include/PreProcessor.h"

#include "Misha/CmdLineParser.h"
#include "Misha/Geometry.h"
#include "Misha/Image.h"
#include "Misha/Miscellany.h"
#include <Misha/SquaredEDT.h>

#include "Include/LocalFit.h"
#include "Include/ECD.h"


using namespace MishaK;

CmdLineParameterArray< std::string , 2 >
	In( "in" ) ,
	Mask( "mask" );

CmdLineParameterArray< unsigned int , 2 >
	StratifiedSamplingResolution( "stratified" , DefaultStratifiedSamplingResolution ) ,
	ECDNum( "ecdNum" , DefaultECDFeatureNum );

CmdLineParameter< unsigned int >
	ECDHalfResolution( "ecdHalfRes" , DefaultECDHalfResolution );

CmdLineParameter< double > 
	FitDeviation( "fitDev" , DefaultDeviation ),
	ECDRadius( "ecdRadius" );

CmdLineReadable
	NoColor( "noColor" ) ,
	NoNormalize( "noNormalize" ) ,
	KeepUnsupported( "keepUnsupported" ) ,
	Verbose( "verbose" );

std::vector< CmdLineReadable * > params =
{
	&In ,
	&Mask ,
	&FitDeviation ,
	&ECDHalfResolution ,
	&ECDRadius ,
	&ECDNum ,
	&NoColor ,
	&StratifiedSamplingResolution ,
	&NoNormalize ,
	&KeepUnsupported ,
	&Verbose ,
};

void ShowUsage( const char* ex )
{
	printf( "Usage %s:\n" , ex );
	printf( "\t --%s <input source> <input target>\n" , In.name.c_str() );
	printf( "\t[--%s <input source mask> <input target mask>]\n" , Mask.name.c_str() );
	printf( "\t[--%s <ECD radius resolution>=%d]\n" , ECDHalfResolution.name.c_str() , ECDHalfResolution.value );
	printf( "\t[--%s <ECD radius>=<ECD half resolution>*%f]\n" , ECDRadius.name.c_str() , ECDResolutionToDiameter );
	printf( "\t[--%s <fit deviation>=%f]\n" , FitDeviation.name.c_str() , FitDeviation.value );
	printf( "\t[--%s <number of source extended convolution descriptors, number of target extended convolution descriptors>=%d %d]\n" , ECDNum.name.c_str() , ECDNum.values[0] , ECDNum.values[1] );
	printf( "\t[--%s <stratified source sampling resolution, stratified target sampling resolution>=%d %d]\n" , StratifiedSamplingResolution.name.c_str() , StratifiedSamplingResolution.values[0] , StratifiedSamplingResolution.values[1] );
	printf( "\t[--%s]\n" , NoColor.name.c_str() );
	printf( "\t[--%s]\n" , NoNormalize.name.c_str() );
	printf( "\t[--%s]\n" , KeepUnsupported.name.c_str() );
	printf( "\t[--%s]\n" , Verbose.name.c_str() );
}

template< bool UseColor , unsigned int ECDHalfResolution >
using Descriptor = std::conditional_t< UseColor , ECD::Descriptor< ECDHalfResolution , Point< double , 3 > > , ECD::Descriptor< ECDHalfResolution > >;

RegularGrid< 2 , bool > ReadMask( std::string fileName , bool flip=false )
{
	unsigned int res[2];
	unsigned char * clr = ImageReader<>::ReadColor( fileName , res[0] , res[1] );
	RegularGrid< 2 , bool > mask( res );
	for( unsigned int j=0 , idx=0 ; j<res[1] ; j++ ) for( unsigned int i=0 ; i<res[0] ; i++ , idx++ )
	{
		bool value = ( clr[idx*3+0] || clr[idx*3+1] || clr[idx*3+2] ) ? true : false;
		if( flip ) mask(i,mask.res(1)-1-j) = value;
		else       mask(i,j) = value;
	}
	return mask;
}

RegularGrid< 2 , double > ReadGrid( std::string fileName )
{
	unsigned int res[2];
	unsigned char * clr = ImageReader<>::ReadColor( fileName , res[0] , res[1] );
	RegularGrid< 2 , double > grid( res );
	for( unsigned int j=0 , idx=0 ; j<res[1] ; j++ ) for( unsigned int i=0 ; i<res[0] ; i++ , idx++ )
	{
		double value = static_cast< double >( clr[idx*3+0] ) + static_cast< double >( clr[idx*3+1] ) + static_cast< double >( clr[idx*3+2] );
		grid(i,j) = value / ( 255. * 3 );
	}
	return grid;
}

RegularGrid< 2 , Point< double , 3 > > ReadColorGrid( std::string fileName , bool flip=false )
{
	unsigned int res[2];
	unsigned char * clr = ImageReader<>::ReadColor( fileName , res[0] , res[1] );
	RegularGrid< 2 , Point< double , 3 > > grid( res );
	for( unsigned int j=0 , idx=0 ; j<res[1] ; j++ ) for( unsigned int i=0 ; i<res[0] ; i++ , idx++ ) for( unsigned int c=0 ; c<3 ; c++ )
		if( flip ) grid(i,grid.res(1)-1-j)[c] = clr[idx*3+c] / 255.;
		else       grid(i,              j)[c] = clr[idx*3+c] / 255.;
	return grid;
}

RegularGrid< 2 , double > ColorToGray( const RegularGrid< 2 , Point< double , 3 > > & colors )
{
	RegularGrid< 2 , double > gray( colors.res() );
	for( size_t i=0 ; i<colors.resolution() ; i++ ) gray[i] = Point< double , 3 >::Dot( Point< double , 3 >(1.,1.,1.) , colors[i] ) / 3.;
	return gray;
}

template< bool UseColor , unsigned int ECDHalfResolution >
std::vector< Descriptor< UseColor , ECDHalfResolution > > GetFeatures( std::string inFile , std::string maskFile , double fitDev , double fitRadius , double ecdRadius , unsigned int ecdNum , unsigned int stratificationRes , bool removeUnsupported )
{
	using Desc = Descriptor< UseColor , ECDHalfResolution >;

	RegularGrid< 2 , double > values = ReadGrid( inFile );
	RegularGrid< 2 , Point< double , 2 > > grads;

	struct InverseMaskWrapper
	{
		const RegularGrid< 2 , bool > & mask;
		InverseMaskWrapper( const RegularGrid< 2 , bool > & mask ) : mask(mask){}
		unsigned int res( unsigned int d ) const { return mask.res(d); }
		bool operator[]( size_t idx ) const { return mask[idx]==0; };
	};

	RegularGrid< 2 , bool > mask;

	double radius2 = ecdRadius*ecdRadius;

	if( maskFile.size() )
	{
		mask = ReadMask( maskFile );
		if( values.res(0)!=mask.res(0) || values.res(1)!=mask.res(1) )
			MK_ERROR_OUT( "Image and mask resolutions don't match: " , values.res(0) , " x " , values.res(1) , " != " , mask.res(0) , " x " , mask.res(1) );
		grads = LocalFit::Gradients( values , [&]( unsigned int x , unsigned int y ){ return mask(x,y) ? 1. : 0. ; } , fitDev , fitRadius );
		for( size_t i=0 ; i<grads.resolution() ; i++ ) grads[i] *= mask[i] ? 1. : 0.;
	}
	else grads = LocalFit::Gradients( values , [&]( unsigned int , unsigned int ){ return 1.; } , fitDev , fitRadius );

	if constexpr( UseColor )
	{
		RegularGrid< 2 , Desc > descriptors = Desc::Get( grads , ReadColorGrid( inFile ) , ecdRadius );
		if( removeUnsupported && mask.resolution() )
		{
			RegularGrid< 2 , unsigned int > edt = SquaredEDT< bool , 2 >::Saito( InverseMaskWrapper(mask) );
			for( size_t i=0 ; i<edt.resolution() ; i++ ) if( edt[i]<=radius2 ) descriptors[i].clear();
		}
		if( stratificationRes ) return ECD::GetTopStratified< ECDHalfResolution >( descriptors , stratificationRes , ecdNum );
		else                    return ECD::GetTop< ECDHalfResolution >( descriptors , ecdNum );
	}
	else
	{
		RegularGrid< 2 , Desc > descriptors = Desc::Get( grads , ecdRadius );
		if( removeUnsupported && mask.resolution() )
		{
			RegularGrid< 2 , unsigned int > edt = SquaredEDT< bool , 2 >::Saito( InverseMaskWrapper(mask) );
			for( size_t i=0 ; i<edt.resolution() ; i++ ) if( edt[i]<=radius2 ) descriptors[i].clear();
		}
		if( stratificationRes ) return ECD::GetTopStratified< ECDHalfResolution >( descriptors , stratificationRes , ecdNum );
		else                    return ECD::GetTop< ECDHalfResolution >( descriptors , ecdNum );
	}
}

template< bool UseColor , unsigned int ECDHalfResolution >
int Execute( void )
{
	using ECDDescriptor = Descriptor< UseColor , ECDHalfResolution >;

	static const unsigned int DescriptorEntries = sizeof( ECDDescriptor ) / sizeof( double );

	Miscellany::Timer timer;

	std::vector< ECDDescriptor > sourceFeatures;
	{
		timer.reset();
		sourceFeatures = GetFeatures< UseColor , ECDHalfResolution >( In.values[0] , Mask.set ? Mask.values[0] : std::string() , FitDeviation.value , FitDeviation.value * DeviationToRadius , ECDRadius.value , ECDNum.values[0] , StratifiedSamplingResolution.values[0] , !KeepUnsupported.set );
		if( !NoNormalize.set ) for( unsigned int i=0 ; i<sourceFeatures.size() ; i++ ) sourceFeatures[i] = sourceFeatures[i].normalized();
		if( Verbose.set ) std::cout << "Source features ( " << sourceFeatures.size() << " ): " << timer() << std::endl;
	}

	timer.reset();

	std::vector< std::string > targets = Tokenize( In.values[1] , ':' );
	std::vector< std::string > targetMasks;
	if( Mask.set )
	{
		targetMasks = Tokenize( Mask.values[1] , ':' );
		if( targetMasks.size()!=targets.size() ) MK_ERROR_OUT( "Number of targets and masks don't match: " , targets.size() , " != " , targetMasks.size() );
	}
	for( unsigned int i=0 ; i<targets.size() ; i++ )
	{
		timer.reset();

		std::vector< ECDDescriptor > targetFeatures = GetFeatures< UseColor , ECDHalfResolution >( targets[i] , Mask.set ? targetMasks[i] : std::string() , FitDeviation.value , FitDeviation.value * DeviationToRadius , ECDRadius.value , ECDNum.values[1] , StratifiedSamplingResolution.values[1] , !KeepUnsupported.set );
		if( !NoNormalize.set ) for( unsigned int i=0 ; i<targetFeatures.size() ; i++ ) targetFeatures[i] = targetFeatures[i].normalized();
		if( Verbose.set )
			if( targets.size()>1 ) std::cout << "Target features[" << i << "] ( " << targetFeatures.size() << " ): " << timer() << std::endl;
			else                   std::cout << "Target features ( " << targetFeatures.size() << " ): " << timer() << std::endl;

		std::atomic< double > _err2 = 0;
		ThreadPool::ParallelFor
			(
				0 , targetFeatures.size() , 
				[&]( size_t i )
				{
					unsigned int sIdx = -1;
					double diff = std::numeric_limits< double >::infinity();

					for( unsigned int j=0 ; j<sourceFeatures.size() ; j++ ) 
					{
						double d = ECDDescriptor::Difference( sourceFeatures[j] , targetFeatures[i] );
						if( d<diff ) sIdx = j , diff = d;
					}
					if( sIdx==-1 ) MK_ERROR_OUT( "Could not find a source feature" );
					_err2 += diff;
				}
			);

		if( Verbose.set )
			if( targets.size()>1 ) std::cout << "Closest features[" << i << " ] ( " << _err2 << " ) " << timer() << std::endl;
			else                   std::cout << "Closest features ( " << _err2 << " ) " << timer() << std::endl;
		else std::cout << _err2 << std::endl;
	}

	return EXIT_SUCCESS;
}

template< bool UseColor , unsigned int Radius , unsigned int MaxRadius >
void Execute( unsigned int radius )
{
	if constexpr( Radius>MaxRadius ) MK_ERROR_OUT( "Max radius exceeded: " , radius , " : " , MaxRadius );
	else
	{
		if( radius==Radius ) Execute< UseColor , Radius >();
		else                 Execute< UseColor , Radius+1 , MaxRadius >( radius );
	}
}

int main( int argc , char* argv[] )
{
	CmdLineParse( argc-1 , argv+1 , params );

	if( !In.set )
	{
		ShowUsage( argv[0] );
		return EXIT_FAILURE;
	}

	if( !ECDRadius.set ) ECDRadius.value = ECDHalfResolution.value * ECDResolutionToDiameter;

	if( NoColor.set ) Execute< false , 0 , MaxECDHalfResolution >( ECDHalfResolution.value );
	else              Execute< true  , 0 , MaxECDHalfResolution >( ECDHalfResolution.value );

	return EXIT_SUCCESS;
}
