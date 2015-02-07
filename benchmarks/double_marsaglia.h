# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include "../MT19937.h"


typedef struct ke_we_t {
        int64_t ke;
            double we;
} ke_we_t;

#ifdef BIT_OPERATIONS
#define RECURSIVE_TAIL
#endif
#ifdef RECURSIVE_TAIL
#define PAIRED
#endif

static double fe[256];
#ifdef PAIRED
static ke_we_t ke_we[256];
#define KE(i) (ke_we[(i)].ke)
#define WE(i) (ke_we[(i)].we)
#else
static int64_t ke[256];
static double we[256];
#define KE(i) (ke[(i)])
#define WE(i) (we[(i)])
#endif

#ifdef RECURSIVE_TAIL
#define EXP_TAIL() ( r8_exp() )
#else
#define EXP_TAIL() ( -log(uniform_double_PRN()) )
#endif

#ifdef BIT_OPERATIONS
#define labs( jz ) ( (jz) & 0x7fffffffffffffff ) 
#endif



/******************************************************************************/

double r8_exp ( )

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

  MT_FLUSH();
  jz = Rand++->sl;
  iz = ( jz & 255 );

  if ( labs ( jz  ) < KE(iz) )
  {
    return ( double ) ( labs ( jz ) ) * WE(iz);
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 ) return 7.69711 + EXP_TAIL();

      x = ( double ) ( labs ( jz ) ) * WE(iz);

      if ( fe[iz] + uniform_double_PRN() * ( fe[iz-1] - fe[iz] ) < exp ( - x ) )
      {
        return x;
      }

      MT_FLUSH();
      jz =  Rand++->sl;
      iz = ( jz & 255 );

      if ( labs ( jz ) < KE(iz) )
      {
        return ( double ) ( labs ( jz ) ) * WE(iz);
      }
    }
  }
}
/******************************************************************************/

void r8_exp_setup ( )

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
  const double m4 = 9223372036854775808.0; 
  double q;
  double te = 7.697117470131487;
  const double ve = 3.949659822581572E-03;

  q = ve / exp ( - de );

  KE(0) = ( int64_t ) ( ( de / q ) * m4 );
  KE(1) = 0;

  WE(0) = ( double ) ( q / m4 );
  WE(255) = ( double ) ( de / m4 );

  fe[0] = 1.0;
  fe[255] = ( double ) ( exp ( - de ) );

  for ( i = 254; 1 <= i; i-- )
  {
    de = - log ( ve / de + exp ( - de ) );
    KE(i+1) = ( double ) ( ( de / te ) * m4 );
    te = de;
    fe[i] = ( double ) ( exp ( - de ) );
    WE(i) = ( double ) ( de / m4 );
  }
  return;
}
