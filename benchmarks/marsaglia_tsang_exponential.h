# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>

# include "MT19937.h"

/******************************************************************************/

float r4_exp ( unsigned long int *jsr, int ke[256], float fe[256], 
  float we[256] )

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

    Input/output, unsigned long int *JSR, the seed.

    Input, int KE[256], data computed by R4_EXP_SETUP.

    Input, float FE[256], WE[256], data computed by R4_EXP_SETUP.

    Output, float R4_EXP, an exponentially distributed random value.
*/
{
  int iz;
  int jz;
  float value;
  float x;

  /* Using SFMT in lieu of native uniform random number generator */
  MT_FLUSH();
  jz = /* shr3 ( jsr ) */ Rand++->i[0];
  iz = ( jz & 255 );

  if ( abs ( jz  ) < ke[iz] )
  {
    value = ( float ) ( abs ( jz ) ) * we[iz];
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        value = 7.69711 - log ( /* r4_uni ( jsr ) */ uniform_double() );
        break;
      }

      x = ( float ) ( abs ( jz ) ) * we[iz];

      if ( fe[iz] + /* r4_uni ( jsr ) */ uniform_double() * ( fe[iz-1] - fe[iz] ) < exp ( - x ) )
      {
        value = x;
        break;
      }

      MT_FLUSH();
      jz = /* shr3 ( jsr ) */ Rand++->i[0];
      iz = ( jz & 255 );

      if ( abs ( jz ) < ke[iz] )
      {
        value = ( float ) ( abs ( jz ) ) * we[iz];
        break;
      }
    }
  }
  return value;
}
/******************************************************************************/

void r4_exp_setup ( int ke[256], float fe[256], float we[256] )

/******************************************************************************/
/*
  Purpose:

    R4_EXP_SETUP sets data needed by R4_EXP.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    09 December 2008

  Author:

    John Burkardt

  Reference:

    George Marsaglia, Wai Wan Tsang,
    The Ziggurat Method for Generating Random Variables,
    Journal of Statistical Software,
    Volume 5, Number 8, October 2000, seven pages.

  Parameters:

    Output, int KE[256], data needed by R4_EXP.

    Output, float FE[256], WE[256], data needed by R4_EXP.
*/
{
  double de = 7.697117470131487;
  int i;
  const double m2 = 2147483648.0;
  double q;
  double te = 7.697117470131487;
  const double ve = 3.949659822581572E-03;

  q = ve / exp ( - de );

  ke[0] = ( int ) ( ( de / q ) * m2 );
  ke[1] = 0;

  we[0] = ( float ) ( q / m2 );
  we[255] = ( float ) ( de / m2 );

  fe[0] = 1.0;
  fe[255] = ( float ) ( exp ( - de ) );

  for ( i = 254; 1 <= i; i-- )
  {
    de = - log ( ve / de + exp ( - de ) );
    ke[i+1] = ( int ) ( ( de / te ) * m2 );
    te = de;
    fe[i] = ( float ) ( exp ( - de ) );
    we[i] = ( float ) ( de / m2 );
  }
  return;
}

