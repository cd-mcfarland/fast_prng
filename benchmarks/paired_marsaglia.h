# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
# include "../MT19937.h"

typedef struct ke_we_t {        // The name for a haplotype in c is haplo-types, get it?
    int64_t ke;       // Every haplotype has a parent haplotype, thus defining the structure of the tree. 
    double we;                    // Every haplotype has a fitness that differs from its parent by ds. 
} ke_we_t;

/******************************************************************************/

double r8_exp ( ke_we_t ke_we[256], double fe[256] )

/******************************************************************************/
/*
  Purpose:

    R4_EXP returns an exponentially distributed single precision real value.
    
*/
{
  int iz;
  int64_t jz;
  double x;

  MT_FLUSH();
  jz = Rand++->l;
  iz = ( jz & 255 );

  if ( abs ( jz  ) < ke_we[iz].ke )
  {
    return ( double ) ( abs ( jz ) ) * ke_we[iz].we;
  }
  else
  {
    for ( ; ; )
    {
      if ( iz == 0 )
      {
        return 7.69711 - log ( uniform_double_PRN() );
      }

      x = ( double ) ( abs ( jz ) ) * ke_we[iz].we;

      if ( fe[iz] + uniform_double_PRN() * ( fe[iz-1] - fe[iz] ) < exp ( - x ) )
      {
        return x;
      }

      MT_FLUSH();
      jz = Rand++->i[0];
      iz = ( jz & 255 );

      if ( abs ( jz ) < ke_we[iz].ke )
      {
        return ( double ) ( abs ( jz ) ) * ke_we[iz].we;
      }
    }
  }
}
/******************************************************************************/

void r8_exp_setup ( ke_we_t ke_we[256], double fe[256] )

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

  ke_we[0].ke = ( int64_t ) ( ( de / q ) * m2 );
  ke_we[1].ke = 0;

  ke_we[0].we = ( double ) ( q / m2 );
  ke_we[255].we = ( double ) ( de / m2 );

  fe[0] = 1.0;
  fe[255] = ( double ) ( exp ( - de ) );

  for ( i = 254; 1 <= i; i-- )
  {
    de = - log ( ve / de + exp ( - de ) );
    ke_we[i+1].ke = ( double ) ( ( de / te ) * m2 );
    te = de;
    fe[i] = ( double ) ( exp ( - de ) );
    ke_we[i].we = ( double ) ( de / m2 );
  }
  return;
}

