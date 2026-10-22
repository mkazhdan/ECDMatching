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

template< HasWeightFunction WeightFunction >
MishaK::RegularGrid< 2 , MishaK::Point< double , 3 > > Linear( const MishaK::RegularGrid< 2 , double > & scalars , const WeightFunction & weights , double dev , double radius )
{
	MishaK::RegularGrid< 2 , MishaK::Point< double , 3 > > fit( scalars.res() );

	int iRadius = static_cast< unsigned int >( std::ceil( radius ) );
	MishaK::ThreadPool::ParallelFor
	(
		0 , scalars.resolution() ,
		[&]( size_t idx )
		{
			int x = static_cast< int >( idx % scalars.res(0) );
			int y = static_cast< int >( idx / scalars.res(0) );

			MishaK::SquareMatrix< double , 3 > A;
			MishaK::Point< double , 3 > b;
			for( int dx=-iRadius ; dx<=iRadius ; dx++ ) if( x+dx>=0 && x+dx<static_cast< int >( scalars.res(0) ) )
				for( int dy=-iRadius ; dy<=iRadius ; dy++ ) if( y+dy>=0 && y+dy<static_cast< int >( scalars.res(1) ) )
				{
					// E( v , x , y ) = \sum_i w_i * [ v_i - ( v + x * x_i + y * y_i ) ]^2
					// dE / dv = 2 * \sum_i w_i * [ v_i - ( v + x * x_i + y * y_i ) ] * ( - 1 )
					//         = 2 * [ - \sum_i w_i * v_i       + \sum_i w_i * v       + \sum_i w_i * x_i       * x + \sum_i w_i *       y_i * y ]
					// dE / dx = 2 * \sum_i w_i * [ v_i - ( v + x * x_i + y * y_i ) ] * ( - x_i )
					//         = 2 * [ - \sum_i w_i * v_i * x_i + \sum_i w_i * v * x_i + \sum_i w_i * x_i * x_i * x + \sum_i w_i * x_i * y_i * y ]
					// dE / dy = 2 * \sum_i w_i * [ v_i - ( v + x * x_i + y * y_i ) ] * ( - y_i )
					//         = 2 * [ - \sum_i w_i * v_i * y_i + \sum_i w_i * v * y_i + \sum_i w_i * x_i * y_i * x + \sum_i w_i * y_i * y_i * y ]
					// Setting the gradient to zero:
					//  \sum_i w_i * v_i       = \sum_i w_i       * v + \sum_i w_i * x_i       * x + \sum_i w_i * y_i       * y
					//  \sum_i w_i * v_i * x_i = \sum_i w_i * x_i * v + \sum_i w_i * x_i^2     * x + \sum_i w_i * x_i * y_i * y
					//  \sum_i w_i * v_i * y_i = \sum_i w_i * y_i * v + \sum_i w_i * x_i * y_i * x + \sum_i w_i * y_i^2     * y

					double dist2 = dx*dx + dy*dy;
					if( dist2<=radius*radius )
					{
						double w = weights( x+dx , y+dy ) * exp( -dist2 / ( 2. * dev * dev ) );
						double v = scalars( x+dx , y+dy );

						MishaK::Point< double , 3 > _b( 1. , dx , dy );
						for( unsigned int i=0 ; i<3 ; i++ )
						{
							b[i] += w * v * _b[i];
							for( unsigned int j=0 ; j<3 ; j++ ) A(i,j) += w * _b[i] * _b[j];
						}
					}
				}		
			bool invertible;
			MishaK::Point< double , 3 > g = A.inverse(invertible) * b;
			if( invertible ) fit[idx] = g;
		}
	);

	return fit;
}


template< HasWeightFunction WeightFunction >
MishaK::RegularGrid< 2 , MishaK::Point< double , 2 > > Gradients( const MishaK::RegularGrid< 2 , double > & scalars , const WeightFunction & weights , double dev , double radius )
{
	MishaK::RegularGrid< 2 , MishaK::Point< double , 3 > > linearFit = Linear( scalars , weights , dev , radius );
	MishaK::RegularGrid< 2 , MishaK::Point< double , 2 > > grads;

	grads.resize( linearFit.res() );
	for( size_t i=0 ; i<linearFit.resolution() ; i++ ) grads[i] = MishaK::Point< double , 2 >( linearFit[i][1] , linearFit[i][2] );

	return grads;
}
