# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include "../MT19937.h"

/******************************************************************************/

double r8_exp ( /* unsigned int64_t int *jsr, */ int64_t ke[256], double fe[256], 
  double we[256] )

/******************************************************************************/
/*
  Purpose:

    R4_EXP returns an exponentially distributed single precision real value.

  Discussion:

    The underlying algorithm is the ziggurat method.

    Before the first call to this function, the user must call R4_EXP_SETUP
    to determine the values of KE, FE and WE.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 20080

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Input/output, unsigned int64_t int *JSR, the seed.

    Input, int KE[256], data computed by R4_EXP_SETUP.

    Input, float FE[256], WE[256], data computed by R4_EXP_SETUP.

    Output, float R4_EXP, an exponentially distributed random value.
*/
{
  int iz;
  int64_t jz;
  double x;

  /* Using SFMT in lieu of native uniform random number generator */
  MT_FLUSH();
  jz = /* shr3 ( jsr ) */ Rand++->l;
  iz = ( jz & 255 );

  if ( abs ( jz  ) < ke[iz] )
  {
    return ( double ) ( abs ( jz ) ) * we[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        return 7.69711 - log ( /* r4_uni ( jsr ) */ uniform_double_PRN() );
      }

      x = ( double ) ( abs ( jz ) ) * we[iz];

      if ( fe[iz] + /* r4_uni ( jsr ) */ uniform_double_PRN() * ( fe[iz-1] - fe[iz] ) < exp ( - x ) )
      {
        return x;
      }

      MT_FLUSH();
      jz = /* shr3 ( jsr ) */ Rand++->i[0];
      iz = ( jz & 255 );

      if ( abs ( jz ) < ke[iz] )
      {
        return ( double ) ( abs ( jz ) ) * we[iz];
      }
    }
  }
}
/******************************************************************************/

void r8_exp_setup ( int64_t ke[256], double fe[256], double we[256] )

/******************************************************************************/
/*
  Purpose:

    R8_EXP_SETUP sets data needed by R8_EXP.

    The constants in this algorithm are not necessarily calculated to full double
    precision, although rounding errors should remain small. This algorithm is
    intended for benchmarking only. 

*/
{
  double de = 7.697117470131487;
  int i;
  const double m2 = 2147483648.0; 
  double q;
  double te = 7.697117470131487;
  const double ve = 3.949659822581572E-03;

  q = ve / exp ( - de );

  ke[0] = ( int64_t ) ( ( de / q ) * m2 );
  ke[1] = 0;

  we[0] = ( double ) ( q / m2 );
  we[255] = ( double ) ( de / m2 );

  fe[0] = 1.0;
  fe[255] = ( double ) ( exp ( - de ) );

  for ( i = 254; 1 <= i; i-- )
  {
    de = - log ( ve / de + exp ( - de ) );
    ke[i+1] = ( double ) ( ( de / te ) * m2 );
    te = de;
    fe[i] = ( double ) ( exp ( - de ) );
    we[i] = ( double ) ( de / m2 );
  }
  return;
}

