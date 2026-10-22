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

template< unsigned int HalfRes , typename ... Data >
Descriptor< HalfRes , Data ... >::Descriptor( void )
{
	Base::first = 0;
	for( unsigned int i=0 ; i<Dim ; i++ ) Base::second[i] = {};
}

template< unsigned int HalfRes , typename ... Data >
template< unsigned int Idx , bool ThreadSafe >
void Descriptor< HalfRes , Data ... >::_add( size_t idx , std::tuple< double , Data... > entry )
{
	if constexpr( Idx<sizeof...(Data)+1 )
	{
		if constexpr( ThreadSafe ) MishaK::Atomic< std::tuple_element_t< Idx , std::tuple< double , Data... > > >::Add( std::get< Idx >( Base::second.operator[]( idx ) ) , std::get< Idx >( entry ) );
		else std::get< Idx >( Base::second.operator[]( idx ) ) += std::get< Idx >( entry );
		_add< Idx+1 , ThreadSafe >( idx , entry );
	}
}

template< unsigned int HalfRes , typename ... Data >
template< bool ThreadSafe >
void Descriptor< HalfRes , Data ... >::_splat( double x , double y , double weight , Data ... data )
{
	x *= HalfRes , y *= HalfRes;

	int ix1 = static_cast< int >( std::floor( x ) );
	int iy1 = static_cast< int >( std::floor( y ) );
	double dx2 = x - ix1;
	double dx1 = 1. - dx2;
	double dy2 = y - iy1;
	double dy1 = 1. - dy2;

	ix1 += HalfRes;
	iy1 += HalfRes;
	int ix2 = ix1 + 1;
	int iy2 = iy1 + 1;

	auto Add = [&]( int ix , int iy , double dxdy )
		{
			if( ix>=0 && ix<=static_cast< int >(Width) ) if( iy>=0 && iy<=static_cast< int >(Width) )
			{
				double w = weight * dxdy;
				if constexpr( ThreadSafe ) MishaK::AddAtomic( Base::first , w );
				else                       Base::first += w;
				_add< 0 , ThreadSafe >( iy*Width + ix , std::tuple< double , Data... >( w ,  ( data * w )... ) );
			}
		};

	Add( ix1 , iy1 , dx1*dy1 );
	Add( ix1 , iy2 , dx1*dy2 );
	Add( ix2 , iy1 , dx2*dy1 );
	Add( ix2 , iy2 , dx2*dy2 );
}

template< unsigned int HalfRes , typename ... Data >
Descriptor< HalfRes , Data ... > & Descriptor< HalfRes , Data ... >::clear( void )
{
	Base::first = 0.;
	for( unsigned int d=0 ; d<Dim ; d++ ) Base::second[d] = {};
	return *this;
}

template< unsigned int HalfRes , typename ... Data >
double Descriptor< HalfRes , Data ... >::magnitude( void ) const
{
	return Base::first;
}

template< unsigned int HalfRes , typename ... Data >
template< unsigned int Idx >
void Descriptor< HalfRes , Data ... >::_scale( double scl )
{
	if constexpr( Idx<sizeof...(Data)+1 )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) std::get< Idx >( Base::second[d] ) *= scl;
		_scale< Idx+1 >( scl );
	}
}

template< unsigned int HalfRes , typename ... Data >
Descriptor< HalfRes , Data ... > Descriptor< HalfRes , Data ... >::normalized( void ) const
{
	if( Base::first )
	{
		Descriptor out = *this;
		out.template _scale< 0 >( 1./Base::first );
		return out;
	}
	else return Descriptor();
}

template< unsigned int HalfRes , typename ... Data >
template< unsigned int Idx >
double Descriptor< HalfRes , Data ... >::_Difference( const Descriptor & desc1 , const Descriptor & desc2 )
{
	double l2 = 0;
	if constexpr( Idx==0 )
	{
		for( unsigned int d=0 ; d<Dim ; d++ )
		{
			double diff = std::get< Idx >( desc1.second[d] ) - std::get< Idx >( desc2.second[d] );
			l2 += diff * diff;
		}
		l2 += _Difference< Idx+1 >( desc1 , desc2 );
	}
	else if constexpr( Idx<sizeof...(Data)+1 )
	{
		for( unsigned int d=0 ; d<Dim ; d++ ) l2 += ( std::get< Idx >( desc1.second[d] ) - std::get< Idx >( desc2.second[d] ) ).squareNorm();
		l2 += _Difference< Idx+1 >( desc1 , desc2 );
	}
	return l2;
}

template< unsigned int HalfRes , typename ... Data >
double Descriptor< HalfRes , Data ... >::Difference( const Descriptor & desc1 , const Descriptor & desc2 )
{
	return _Difference< 0 >( desc1 , desc2 );
}

template< unsigned int HalfRes , typename ... Data >
template< unsigned int Idx >
void Descriptor< HalfRes , Data ... >::_CheckDataGridResolutions( std::tuple< const MishaK::RegularGrid< 2 , Data > & ... > dataGrids , unsigned int width , unsigned int height )
{
	if constexpr( Idx<sizeof...(Data) )
	{
		if( std::get< Idx >( dataGrids ).res(0)!=width || std::get< Idx >( dataGrids ).res(1)!=height )
			MK_THROW( "Data grid has bad resolution: " , std::get< Idx >( dataGrids ).res(0) , " x " , std::get< Idx >( dataGrids ).res(1) , " != " , width , " x " , height );
		_CheckDataGridResolutions< Idx+1 >( dataGrids , width , height );
	}
}

template< unsigned int HalfRes , typename ... Data >
MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > Descriptor< HalfRes , Data ... >::Get( const MishaK::RegularGrid< 2 , MishaK::Point< double , 2 > > & grads , const MishaK::RegularGrid< 2 , Data > & ... values , double radius , bool gather )
{
	_CheckDataGridResolutions< 0 >( std::tie( values... ) , grads.res(0) , grads.res(1) );

	MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > ecd( grads.res() );

	int iRadius = static_cast< int >( std::ceil( radius ) );

	if( gather )
	{
		MishaK::RegularGrid< 2 , MishaK::ProjectiveData< double , MishaK::SquareMatrix< double , 2 > > > wFrames( grads.res() );
		MishaK::RegularGrid< 2 , unsigned char > mask( grads.res() );
		{
			double radius2 = radius * radius;

			MishaK::ThreadPool::ParallelFor
			(
				0 , grads.resolution() ,
				[&]( size_t idx )
				{
					int x = static_cast< int >( idx % grads.res(0) );
					int y = static_cast< int >( idx / grads.res(0) );
	
					double l2 = MishaK::Point< double , 2 >::SquareNorm( grads[idx] );
					if( l2 )
					{
						double l = sqrt( l2 );
						MishaK::OrthogonalFrame< 2 > frame( &grads[idx] , 1 );
						for( unsigned int i=0 ; i<2 ; i++ ) for( unsigned int j=0 ; j<2 ; j++ ) wFrames[idx].data(i,j) = frame[j][i];
						wFrames[idx].weight = l;
					}
					else wFrames[idx].weight = 0;

					mask[idx] = 0;
					for( int dy=-iRadius ; dy<=iRadius ; dy++ ) if( y+dy>=0 && y+dy<static_cast< int >( grads.res(1) ) )
					{
						size_t _idx = (y+dy) * grads.res(0) + x;
						for( int dx=-iRadius ; dx<=iRadius ; dx++ ) if( x+dx>=0 && x+dx<static_cast< int >( grads.res(0) ) )
							if( dx*dx+dy*dy<=radius2 && grads[_idx+dx].squareNorm() ) mask[idx] = 1;
					}
				}
			);
		}

		MishaK::ThreadPool::ParallelFor
		(
			0 , grads.resolution() ,
			[&]( size_t idx )
			{
				int x = static_cast< int >( idx % grads.res(0) );
				int y = static_cast< int >( idx / grads.res(0) );
				if( mask[idx] )
				{
					for( int dy=-iRadius ; dy<=iRadius ; dy++ ) if( y+dy>=0 && y+dy<static_cast< int >( grads.res(1) ) )
					{
						size_t _idx = (y+dy) * grads.res(0) + x;
						for( int dx=-iRadius ; dx<=iRadius ; dx++ ) if( x+dx>=0 && x+dx<static_cast< int >( grads.res(0) ) )
						{
							MishaK::Point< double , 2 > p( dx , dy );
							p /= radius;
							if( MishaK::Point< double , 2 >::SquareNorm( p )<=1 )
							{
								p = wFrames[_idx+dx].data * p;
								ecd[idx].template _splat< false >( p[0] , p[1] , wFrames[_idx+dx].weight , ( values[_idx+dx] )... );
							}
						}
					}
				}
			}
		);
	}
	else
	{
		MishaK::ThreadPool::ParallelFor
		(
			0 , grads.resolution() ,
			[&]( size_t idx )
			{
				int x = static_cast< int >( idx % grads.res(0) );
				int y = static_cast< int >( idx / grads.res(0) );

				double l2 = MishaK::Point< double , 2 >::SquareNorm( grads[idx] );
				if( l2 )
				{
					double l = sqrt( l2 );
					MishaK::OrthogonalFrame< 2 > frame( &grads[idx] , 1 );

					for( int dx=-iRadius ; dx<=iRadius ; dx++ ) if( x+dx>=0 && x+dx<static_cast< int >( grads.res(0) ) )
						for( int dy=-iRadius ; dy<=iRadius ; dy++ ) if( y+dy>=0 && y+dy<static_cast< int >( grads.res(1) ) )
						{
							MishaK::Point< double , 2 > p( dx , dy );
							p /= radius;
							if( MishaK::Point< double , 2 >::SquareNorm( p )<=1 )
							{
								p = MishaK::Point< double , 2 >( MishaK::Point< double , 2 >::Dot( p , frame[0] ) , MishaK::Point< double , 2 >::Dot( p , frame[1] ) );

								size_t _idx = (x+dx) + (y+dy) * grads.res(0);
								ecd[_idx].template _splat< true >( -p[0] , -p[1] , l , ( values(x,y) )... );
							}
						}
				}
			}
		);
	}

	return ecd;
}

template< unsigned int HalfRes , typename ... Data >
std::vector< Descriptor< HalfRes , Data ... > > GetTop( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int num )
{
	using Desc = Descriptor< HalfRes , Data ... >;
	auto Cmp = []( const Desc & f1 , const Desc & f2 ){ return f1.magnitude() > f2.magnitude(); };

	std::vector< Desc > top;

	std::atomic< size_t > count = 0;

	MishaK::ThreadPool::ParallelFor( 0 , ecd.resolution() , [&]( size_t idx ){ if( ecd[idx].magnitude() ) count++; } );
	top.resize( count );

	count = 0;
	MishaK::ThreadPool::ParallelFor( 0 , ecd.resolution() , [&]( size_t idx ){ if( ecd[idx].magnitude() ) top[ count++ ] = ecd[idx]; } );

	num = std::min< unsigned int >( num , static_cast< unsigned int >( top.size() ) );
	std::nth_element( top.begin() , top.begin() + num , top.end() , Cmp );
	std::sort( top.begin() , top.begin() + num , Cmp );
	top.resize( num );

	return top;
}

template< unsigned int HalfRes , typename ... Data >
std::vector< typename Descriptor< HalfRes , Data ... >::IndexedDescriptor > GetTopIndexed( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int num )
{
	using IDesc = typename Descriptor< HalfRes , Data ... >::IndexedDescriptor;

	auto Index = [&]( size_t idx ){ return MishaK::template RegularGrid< 2 >::Index( static_cast< int >( idx % ecd.res(0) ) , static_cast< int >( idx / ecd.res(0) ) ); };
	auto Cmp = []( const IDesc & f1 , const IDesc & f2 ){ return f1.first.magnitude() > f2.first.magnitude(); };

	std::vector< IDesc > top;

	std::atomic< size_t > count = 0;

	MishaK::ThreadPool::ParallelFor( 0 , ecd.resolution() , [&]( size_t idx ){ if( ecd[idx].magnitude() ) count++; } );
	top.resize( count );

	count = 0;
	MishaK::ThreadPool::ParallelFor( 0 , ecd.resolution() , [&]( size_t idx ){ if( ecd[idx].magnitude() ) top[ count++ ] = IDesc( ecd[idx] , Index(idx) ); } );

	num = std::min< unsigned int >( num , static_cast< unsigned int >( top.size() ) );
	std::nth_element( top.begin() , top.begin() + num , top.end() , Cmp );
	std::sort( top.begin() , top.begin() + num , Cmp );
	top.resize( num );

	return top;
}

template< unsigned int HalfRes , typename ... Data >
std::vector< Descriptor< HalfRes , Data ... > > GetTopStratified( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int stratificationRes , unsigned int num )
{
	using Desc = Descriptor< HalfRes , Data ... >;
	auto Cmp = []( const Desc & f1 , const Desc & f2 ){ return f1.magnitude() > f2.magnitude(); };

	std::vector< Desc > top;

	unsigned xStratification = ( ecd.res(0)+stratificationRes-1 ) / stratificationRes;
	unsigned yStratification = ( ecd.res(1)+stratificationRes-1 ) / stratificationRes;

	for( unsigned int y=0 ; y<yStratification ; y++ )
	{
		unsigned int yStart = ( y*ecd.res(1) ) / yStratification , yEnd = ( (y+1)*ecd.res(1) ) / yStratification;
		for( unsigned int x=0 ; x<xStratification ; x++ )
		{
			unsigned int xStart = ( x*ecd.res(0) ) / xStratification , xEnd = ( (x+1)*ecd.res(0) ) / xStratification;

			unsigned int idx[] = { static_cast< unsigned int >(-1) , static_cast< unsigned int >(-1) };
			for( unsigned int j=yStart ; j<yEnd ; j++ ) for( unsigned int i=xStart ; i<xEnd ; i++ )
				if( ecd(i,j).magnitude()>0 )
					if( idx[0]==-1 || ecd(i,j).magnitude() > ecd(idx).magnitude() )
						idx[0] = static_cast< int >( i ) , idx[1] = static_cast< int >( j );
			if( idx[0]!=-1 ) top.push_back( ecd(idx) );
		}
	}

	num = std::min< unsigned int >( num , static_cast< unsigned int >( top.size() ) );
	std::nth_element( top.begin() , top.begin() + num , top.end() , Cmp );
	std::sort( top.begin() , top.begin() + num , Cmp );
	top.resize( num );

	return top;
}

template< unsigned int HalfRes , typename ... Data >
std::vector< typename Descriptor< HalfRes , Data ... >::IndexedDescriptor > GetTopStratifiedIndexed( const MishaK::RegularGrid< 2 , Descriptor< HalfRes , Data ... > > & ecd , unsigned int stratificationRes , unsigned int num )
{
	using IDesc = typename Descriptor< HalfRes , Data ... >::IndexedDescriptor;
	auto Cmp = []( const IDesc & f1 , const IDesc & f2 ){ return f1.first.magnitude() > f2.first.magnitude(); };

	std::vector< IDesc > top;

	unsigned xStratification = ( ecd.res(0)+stratificationRes-1 ) / stratificationRes;
	unsigned yStratification = ( ecd.res(1)+stratificationRes-1 ) / stratificationRes;

	for( unsigned int y=0 ; y<yStratification ; y++ )
	{
		unsigned int yStart = ( y*ecd.res(1) ) / yStratification , yEnd = ( (y+1)*ecd.res(1) ) / yStratification;
		for( unsigned int x=0 ; x<xStratification ; x++ )
		{
			unsigned int xStart = ( x*ecd.res(0) ) / xStratification , xEnd = ( (x+1)*ecd.res(0) ) / xStratification;

			MishaK::template RegularGrid< 2 >::Index idx(-1,-1);
			for( unsigned int j=yStart ; j<yEnd ; j++ ) for( unsigned int i=xStart ; i<xEnd ; i++ )
				if( ecd(i,j).magnitude()>0 )
					if( idx[0]==-1 || ecd(i,j).magnitude() > ecd(idx).magnitude() )
						idx[0] = static_cast< int >( i ) , idx[1] = static_cast< int >( j );
			if( idx[0]!=-1 ) top.push_back( IDesc( ecd(idx) , idx ) );
		}
	}

	num = std::min< unsigned int >( num , static_cast< unsigned int >( top.size() ) );
	std::nth_element( top.begin() , top.begin() + num , top.end() , Cmp );
	std::sort( top.begin() , top.begin() + num , Cmp );
	top.resize( num );

	return top;
}
