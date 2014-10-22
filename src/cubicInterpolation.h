
#include <vector>
#include <iostream>

class Spline3Interp 
{

    private:


        double M0,    Mn;
        double coefs[1024][4];
		std::vector<double> xi,yi;

	public:
    
		Spline3Interp( std::vector<double> &xn,  std::vector<double> &yn )
		{
			xi = xn; yi = yn ;  M0 = Mn = 0;
		}

		~Spline3Interp()
		{
		}


		/**
		 * Compute the second derivative at interpolated points.
		 */
		void derivative2( std::vector<double> &dx, std::vector<double> &d1,  std::vector<double> &d2 )
		{
			int N = xi.size(),	M = N-1;
			std::vector<double> b(M),
						 v(M),
						 y(M),
						 alpha(M),
						 beta(M-1);

			for( int i=1; i<M; ++i )
				b[i] = 2 * (dx[i]+dx[i-1]);

			v[1] = 6*(d1[1]-d1[0]) - dx[0]*d2[0];
			for( int i=1; i<M-1; ++i )
				v[i] = 6 * (d1[i]-d1[i-1]);
			v[M-1] = 6*(d1[M-1]-d1[M-2]) - dx[M-1]*d2[M];

			alpha[1] = b[1];
			for( int i=2; i<M; ++i )
				alpha[i] = b[i] - dx[i]*dx[i-1]/alpha[i-1];

			for( int i=1; i<M-1; ++i )
				beta[i] = dx[i]/alpha[i];

			y[1] = v[1]/alpha[1];
			for( int i=2; i<M; ++i )
				y[i] = (v[i]-dx[i]*y[i-1]) / alpha[i];

			d2[M-1] = y[M-1];
			for( int i=M-2; i>0; --i )
				d2[i] = y[i] - beta[i]*d2[i+1];
		}


		/**
		 * Compute the polynomial' coefsficient in each interval.
		 */
		void calcCoefs()
		{
			int N = xi.size(),
				M = N-1;

			std::vector<double> m(N),
						 h(M),
						 d(M);

			m[0] = M0;
			m[M] = Mn;
			for( int i=0; i<M; ++i )
			{
				h[i] = xi[i+1]-xi[i];
				d[i] = (yi[i+1]-yi[i]) / h[i];
			}

			derivative2( h, d, m );

			for( int i=0; i<M; ++i )
			{
				coefs[i][0] = yi[i];
				coefs[i][1] = d[i] - h[i]*(2*m[i]+m[i+1])/6;
				coefs[i][2] = m[i] / 2;
				coefs[i][3] = (m[i+1]-m[i]) / (6*h[i]);
			}
		}


		/**
		 * Compute the value of polynomial at given "x".
		 */
		double evaluate( double x )
		{
			int k = -1,
				N = xi.size(),
				M = N-1;

			double dx,
				 y;

			for( int i=0; i<M; ++i )
			{
				if( (xi[i]<=x) && (xi[i+1]>=x) )
				{
					k = i;
					dx = x-xi[i];
					break;
				}
			}
			if(k!=-1)
			{
				y = ( ( coefs[k][3]*dx + coefs[k][2] ) * dx + coefs[k][1] ) * dx
				  + coefs[k][0];
				return y;
			}
			else
			{
				std::cerr << "The value is out of range!" << std::endl;
				return double(0);
			}
		}


};