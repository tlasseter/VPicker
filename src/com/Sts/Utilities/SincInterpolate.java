package com.Sts.Utilities;
import java.lang.Math;

public class SincInterpolate
{
/* Copyright (c) Colorado School of Mines, 2003.*/
/* All rights reserved.                       */

/*********************** self documentation **********************/
/*****************************************************************************
INTSINC8 - Functions to interpolate uniformly-sampled data via 8-coeff. sinc
		approximations:

ints8r	Interpolation of a uniformly-sampled real function y(x) via a
		table of 8-coefficient sinc approximations

******************************************************************************
Function Prototypes:
void ints8r (int nxin, float dxin, float fxin, float yin[],
	float yinl, float yinr, int nxout, float xout[], float yout[]);

******************************************************************************
Input:
nxin		number of x values at which y(x) is input
dxin		x sampling interval for input y(x)
fxin		x value of first sample input
yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
yinl		value used to extrapolate yin values to left of yin[0]
yinr		value used to extrapolate yin values to right of yin[nxin-1]
nxout		number of x values a which y(x) is output
xout		array[nxout] of x values at which y(x) is output

Output:
yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.

******************************************************************************
Notes:
Because extrapolation of the input function y(x) is defined by the
left and right values yinl and yinr, the xout values are not restricted
to lie within the range of sample locations defined by nxin, dxin, and
fxin.

The maximum error for frequiencies less than 0.6 nyquist is less than
one percent.

******************************************************************************
Author:  Dave Hale, Colorado School of Mines, 06/02/89
*****************************************************************************/
/**************** end self doc ********************************/

/* these are used by both ints8c and ints8r */
    static final int K8 = 8;
    static final int NTABLE = 513;
    static float table[][];
    static boolean tabled=false;

    static public void ints8r (int nxin, float dxin, float fxin, double[] yin,
    	float yinl, float yinr, int nxout, float[] xout, float[] yout) {
    /*****************************************************************************
    Interpolation of a uniformly-sampled real function y(x) via a
    table of 8-coefficient sinc approximations; maximum error for frequiencies
    less than 0.6 nyquist is less than one percent.
    ******************************************************************************
    Input:
    nxin		number of x values at which y(x) is input
    dxin		x sampling interval for input y(x)
    fxin		x value of first sample input
    yin		array[nxin] of input y(x) values:  yin[0] = y(fxin), etc.
    yinl		value used to extrapolate yin values to left of yin[0]
    yinr		value used to extrapolate yin values to right of yin[nxin-1]
    nxout		number of x values a which y(x) is output
    xout		array[nxout] of x values at which y(x) is output

    Output:
    yout		array[nxout] of output y(x):  yout[0] = y(xout[0]), etc.
    ******************************************************************************
    Notes:
    Because extrapolation of the input function y(x) is defined by the
    left and right values yinl and yinr, the xout values are not restricted
    to lie within the range of sample locations defined by nxin, dxin, and
    fxin.
    ******************************************************************************
    Author:  Dave Hale, Colorado School of Mines, 06/02/89
    *****************************************************************************/

    	int jtable;
    	float frac;

    	/* tabulate sinc interpolation coefficients if not already tabulated */
    	if (!tabled) {
                table = new float[513][];
                for (int ia=0; ia < 513; ia++)
                      table[ia]=new float[8];
    		for (jtable=1; jtable<513-1; jtable++) {
    			frac = (float)jtable/(float)(513-1);
    			mksinc(frac,8,jtable);
    		}
    		for (jtable=0; jtable<8; jtable++) {
    			table[0][jtable] = 0.0f;
    			table[513-1][jtable] = 0.0f;
    		}
    		table[0][8/2-1] = 1.0f;
    		table[513-1][8/2] = 1.0f;
    		tabled = true;
    	}

    	/* interpolate using tabulated coefficients */
    	intt8r(513,table,nxin,dxin,fxin,yin,yinl,yinr,nxout,xout,yout);
    }

    /* Copyright (c) Colorado School of Mines, 2003.*/
    /* All rights reserved.                       */

    /*********************** self documentation **********************/
    /*****************************************************************************
    INTTABLE8 -  Interpolation of a uniformly-sampled complex function y(x)
    		via a table of 8-coefficient interpolators

    intt8r	interpolation of a uniformly-sampled real function y(x) via a
    		table of 8-coefficient interpolators

    ******************************************************************************
    Function Prototype:
    void intt8r (int ntable, float table[][8],
    	int nxin, float dxin, float fxin, float yin[],
    	float yinl, float yinr, int nxout, float xout[], float yout[]);

    ******************************************************************************
    Input:
    ntable		number of tabulated interpolation operators; ntable>=2
    table		array of tabulated 8-point interpolation operators
    nxin		number of x values at which y(x) is input
    dxin		x sampling interval for input y(x)
    fxin		x value of first sample input
    yin		array of input y(x) values:  yin[0] = y(fxin), etc.
    yinl		value used to extrapolate yin values to left of yin[0]
    yinr		value used to extrapolate yin values to right of yin[nxin-1]
    nxout		number of x values a which y(x) is output
    xout		array of x values at which y(x) is output

    Output:
    yout		array of output y(x) values:  yout[0] = y(xout[0]), etc.

    ******************************************************************************
    NOTES:
    ntable must not be less than 2.

    The table of interpolation operators must be as follows:

    Let d be the distance, expressed as a fraction of dxin, from a particular
    xout value to the sampled location xin just to the left of xout.  Then,
    for d = 0.0,

    table[0][0:7] = 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0

    are the weights applied to the 8 input samples nearest xout.
    Likewise, for d = 1.0,

    table[ntable-1][0:7] = 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0

    are the weights applied to the 8 input samples nearest xout.  In general,
    for d = (float)itable/(float)(ntable-1), table[itable][0:7] are the
    weights applied to the 8 input samples nearest xout.  If the actual sample
    distance d does not exactly equal one of the values for which interpolators
    are tabulated, then the interpolator corresponding to the nearest value of
    d is used.

    Because extrapolation of the input function y(x) is defined by the left
    and right values yinl and yinr, the xout values are not restricted to lie
    within the range of sample locations defined by nxin, dxin, and fxin.

    ******************************************************************************
    AUTHOR:  Dave Hale, Colorado School of Mines, 06/02/89
    *****************************************************************************/
    /**************** end self doc ********************************/

    static public void intt8r (int ntable, float[][] table,
    	int nxin, float dxin, float fxin, double[] yin, float yinl, float yinr,
    	int nxout, float[] xout, float[] yout)
    {
    	int ioutb,nxinm8,ixout,ixoutn,kyin,ktable,itable;
    	float xoutb,xoutf,xouts,xoutn,frac,fntablem1,yini,sum;


    	/* compute constants */
    	ioutb = -3-8;
    	xoutf = fxin;
    	xouts = 1.0f/dxin;
    	xoutb = 8.0f-xoutf*xouts;
    	fntablem1 = (float)(ntable-1);
    	nxinm8 = nxin-8;
    	//yin0 = &yin[0];
    	//table00 = &table[0][0];

    	/* loop over output samples */
    	for (ixout=0; ixout<nxout; ixout++) {

    		/* determine pointers into table and yin */
    		xoutn = xoutb+xout[ixout]*xouts;
    		ixoutn = (int)xoutn;
    		kyin = ioutb+ixoutn;
    		//pyin = yin0+kyin;
    		frac = xoutn-(float)ixoutn;
                if (frac >= 0.0f)
                    ktable = (int) (frac*fntablem1+0.5f);
                else
                    ktable = (int) ((frac+1.0f)*fntablem1-0.5f);
    		//ptable = table00+ktable*8;

    		/* if totally within input array, use fast method */
    		if (kyin>=0 && kyin<=nxinm8) {
    			yout[ixout] =
    				((float)yin[0+kyin])*table[0][ktable]+
    				((float)yin[1+kyin])*table[1][ktable]+
    				((float)yin[2+kyin])*table[2][ktable]+
    				((float)yin[3+kyin])*table[3][ktable]+
    				((float)yin[4+kyin])*table[4][ktable]+
    				((float)yin[5+kyin])*table[5][ktable]+
    				((float)yin[6+kyin])*table[6][ktable]+
    				((float)yin[7+kyin])*table[7][ktable];

    		/* else handle end effects with care */
    		} else {

    			/* sum over 8 tabulated coefficients */
    			for (itable=0,sum=0.0f; itable<8; itable++,kyin++) {
    				if (kyin<0)
    					yini = yinl;
    				else if (kyin>=nxin)
    					yini = yinr;
    				else
    					yini = (float)yin[kyin];
    				sum += yini*(table[itable][ktable]);
    			}
    			yout[ixout] = sum;
    		}
    	}
    }

    /* Copyright (c) Colorado School of Mines, 2003.*/
    /* All rights reserved.                       */

    /*********************** self documentation **********************/
    /*****************************************************************************
    MKSINC - Compute least-squares optimal sinc interpolation coefficients.

    mksinc		Compute least-squares optimal sinc interpolation coefficients.

    ******************************************************************************
    Function Prototype:
    static void mksinc (float d, int lsinc, float sinc[]);

    ******************************************************************************
    Input:
    d		fractional distance to interpolation point; 0.0<=d<=1.0
    lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20

    Output:
    sinc		array[lsinc] containing interpolation coefficients

    ******************************************************************************
    Notes:
    The coefficients are a least-squares-best approximation to the ideal
    sinc function for frequencies from zero up to a computed maximum
    frequency.  For a given interpolator length, lsinc, mksinc computes
    the maximum frequency, fmax (expressed as a fraction of the nyquist
    frequency), using the following empirically derived relation (from
    a Western Geophysical Technical Memorandum by Ken Larner):

    	fmax = min(0.066+0.265*log(lsinc),1.0)

    Note that fmax increases as lsinc increases, up to a maximum of 1.0.
    Use the coefficients to interpolate a uniformly-sampled function y(i)
    as follows:

                lsinc-1
        y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
                  j=0

    Interpolation error is greatest for d=0.5, but for frequencies less
    than fmax, the error should be less than 1.0 percent.

    ******************************************************************************
    Author:  Dave Hale, Colorado School of Mines, 06/02/89
    *****************************************************************************/
    /**************** end self doc ********************************/


    static void mksinc (float d, int lsinc, int jtable)
    /*****************************************************************************
    Compute least-squares optimal sinc interpolation coefficients.
    ******************************************************************************
    Input:
    d		fractional distance to interpolation point; 0.0<=d<=1.0
    lsinc		length of sinc approximation; lsinc%2==0 and lsinc<=20

    Output:
    sinc		array[lsinc] containing interpolation coefficients
    ******************************************************************************
    Notes:
    The coefficients are a least-squares-best approximation to the ideal
    sinc function for frequencies from zero up to a computed maximum
    frequency.  For a given interpolator length, lsinc, mksinc computes
    the maximum frequency, fmax (expressed as a fraction of the nyquist
    frequency), using the following empirically derived relation (from
    a Western Geophysical Technical Memorandum by Ken Larner):

    	fmax = min(0.066+0.265*log(lsinc),1.0)

    Note that fmax increases as lsinc increases, up to a maximum of 1.0.
    Use the coefficients to interpolate a uniformly-sampled function y(i)
    as follows:

                lsinc-1
        y(i+d) =  sum  sinc[j]*y(i+j+1-lsinc/2)
                  j=0

    Interpolation error is greatest for d=0.5, but for frequencies less
    than fmax, the error should be less than 1.0 percent.
    ******************************************************************************
    Author:  Dave Hale, Colorado School of Mines, 06/02/89
    *****************************************************************************/
    {
    	int j;
    	double[]s;
        double[]a;
        double[]c;
        double[]work;
        double fmax;

        s = new double[20];
        a = new double[20];
        c = new double[20];
        work = new double[20];
    	/* compute auto-correlation and cross-correlation arrays */
    	fmax = 0.066+0.265*Math.log((double)lsinc);
    	fmax = (fmax<1.0)?fmax:1.0;
    	for (j=0; j<lsinc; j++) {
    		a[j] = dsinc(fmax*j);
    		c[j] = dsinc(fmax*(lsinc/2-j-1+d));
    	}

    	/* solve symmetric Toeplitz system for the sinc approximation */
    	stoepd(lsinc,a,c,s,work);
    	for (j=0; j<lsinc; j++)
    		table[jtable][j] = (float) s[j];
    }
    static double dsinc (double x)
    {
        double pix;

        if (x==0.0f) {
                return 1.0f;
        } else {
                pix = Math.PI*x;
                return Math.sin(pix)/pix;
        }
    }
    static void stoepd (int n, double[] r, double[] g, double[] f, double[] a)
    {
        int i,j;
        double v,e,c,w,bot;

        if (r[0] == 0.0) return;

        a[0] = 1.0;
        v = r[0];
        f[0] = g[0]/r[0];

        for (j=1; j<n; j++) {

                /* solve Ra=v as in Claerbout, FGDP, p. 57 */
                a[j] = 0.0;
                f[j] = 0.0;
                for (i=0,e=0.0; i<j; i++)
                        e += a[i]*r[j-i];
                c = e/v;
                v -= c*e;
                for (i=0; i<=j/2; i++) {
                        bot = a[j-i]-c*a[i];
                        a[i] -= c*a[j-i];
                        a[j-i] = bot;
                }

                /* use a and v above to get f[i], i = 0,1,2,...,j */
                for (i=0,w=0.0; i<j; i++)
                        w += f[i]*r[j-i];
                c = (w-g[j])/v;
                for (i=0; i<=j; i++)
                        f[i] -= c*a[j-i];
        }
    }
	static public void resamp(double[] in, float []out, float dtin, float dtout, int ntin, int ntout, float[] rt)
	{
	  int i;
	  int mult;
	  int dti, dto;

	  float tmin;
	  float tv;

	  dti =(int) (dtin  * 1000.f);
	  dto =(int) (dtout * 1000.f);

	  for(i=0,tv=0.0f; i<ntout; i++,tv+=dtout)
		rt[i] = tv;

	  /* Distribute input samples? */
	  if((dto%dti) == 0) {
		mult = dto / dti;
		for(i=0; i<ntout && i*mult<ntin ; i++)
		  out[i] = (float)in[i*mult];
	  }
	  /* sinc interpolate new data */
	  else
		ints8r(ntin,dtin,0.f,in,0.0f,0.0f,ntout,rt,out);

	  if((dti%dto) == 0) {
		mult = dti / dto;
		for(i=0; i<ntin && i*mult<ntout; i++)
		  out[i*mult] = (float)in[i];
	  }

	}

}
