package com.Sts.Utilities;

/**
 * <p>Title: S2S development</p>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.J.Lasseter
 * @version c51c
 */

/** Given an interval from 0 to 1 with points at -1, 0, 1, and 2, fit a cubic curve to these four points and
 *  generate equally-spaced values on a specified number of intervals between 0 and 1.
 *  This class is initialized as either a hermite cubic with value and slope at points 0 and 1, or
 *  as a spline cubic with values at -1, 0, 1, and 2.  From these the coefficients a, b, c, and d of the cubic equation
 *  y = ax**3 + bx**2 + cx + d are computed.
 *  The evaluation of the interval 0 to 1 for N subIntervals, returns N+1 points including the first and last points of the interval.
 */
public class StsCubicForwardDifference
{
	/** number of intervals */
	int nIntervals;
	/** parameter step size raised to powers 1, 2, 3 */
	double du, du2, du3;
	/** initial value */
	double p0;
	/** forward difference parameters */
	double d1, d2, d3;
	/** cubic parameters */
	double a, b, c, d;

	/** if true, compares this solution with cubic calculation for accuracy */
	static final boolean debug = false;

	public StsCubicForwardDifference(int nIntervals)
	{
		this.nIntervals = nIntervals;
		du = 1.0/nIntervals;
		du2 = du*du;
		du3 = du2*du;
	}

    /** initialize the coefficients of the equation y = ax**3 + bx**2 + cx + d */
    public void initialize(double a, double b, double c, double d)
	{
		this.a = a;
		this.b = b;
		this.c = c;
		this.d = d;
		p0 = d;
		d1 = c*du;
		d2 = 2*b*du2;
		d3 = 6*a*du3;
	}

    /** initialize a hermite cubic spline with point and slope at points 0 and 1
     *
     * @param p0 value at point 0
     * @param p1 value at point 1
     * @param s0 slope at point 0
     * @param s1 slope at point 1
     */
    public void hermiteInitialize(double p0, double p1, double s0, double s1)
	{
		double a = 2 * (p0 - p1) + s0 + s1;
		double b = 3 * (p1 - p0) - s1 - 2 * s0;
		double c = s0;
		double d = p0;
		initialize(a, b, c, d);
	}

    /** hermite initialization from an array of doubles */
    public void hermiteInitialize(double[] controlPoints)
	{
		hermiteInitialize(controlPoints[0], controlPoints[1], controlPoints[2], controlPoints[3]);
	}

    /** Initialize a cubic spline with points at -1, 0, 1, and 2
     *
     * @param pm1 value at -1
     * @param p0  value at 0
     * @param p1  value at 1
     * @param p2  value at 2
     */
    public void splineInitialize(double pm1, double p0, double p1, double p2)
	{
		double a = 1.5*(p0 - p1) + 0.5*(p2 - pm1);
		double b = -2.5*p0 + 2*p1 + pm1 -0.5*p2;
		double c = 0.5*(p1 - pm1);
		double d = p0;
		initialize(a, b, c, d);
	}

    /** Evaluate the cubic over the interval 0 to 1 for a given number of subIntervals, nIntervals.
     *
     * @return N+1 values for N subIntervals from 0 to 1 including first and last points
     */
    public double[] evaluate()
	{
		double[] p =  new double[nIntervals + 1];
		p[0] = p0;
		for(int n = 0; n < nIntervals; n++)
		{
			p[n+1] = p[n] + d1 + d2/2 + d3/6;
			d1 = d1 + d2 + d3/2;
			d2 = d2 + d3;
		}
		if(debug)
		{
			double u = du;
			for(int n = 0; n < nIntervals; n++, u += du)
			{
				double value = evalCubicSpline(u);
				if (!StsMath.sameAsTol( (float) value, (float) p[n + 1], 0.001f))
					System.out.println("Error: should be " + value + " is " + p[n + 1]);
			}
		}
		return p;
	}

    /** for debugging and testing */
    double evalCubicSpline(double u)
	{
		return d + u*(c + u*(b + u*a));
	}

    /** for debugging and testing */
    static double evalHermiteSpline(double[] controlPoints, double parameter)
	{
		double y0 = controlPoints[0];
		double y1 = controlPoints[1];
		double s0 = controlPoints[2];
		double s1 = controlPoints[3];
		double a = 2*(y0 - y1) + s0 + s1;
		double b = 3*(y1 - y0) - s1 - 2*s0;
		double c = s0;
		double d = y0;
		double dx = parameter; // *(controlPoints[1][0] - controlPoints[0][0]);
		return d + dx*(c + dx*(b + dx*a));
	}

	static public void main(String[] args)
	{
		int nIntervals = 10;
		double[] controlPoints = new double[] { 10000, 10001, 1, -1 };
		StsCubicForwardDifference fd = new StsCubicForwardDifference(nIntervals);
		fd.hermiteInitialize(controlPoints);
		double[] values = fd.evaluate();
		double du = 1.0f/(nIntervals);
		double u = 0.0f;
		for(int n = 0; n <= nIntervals; n++, u += du)
		{
			double value = fd.evalHermiteSpline(controlPoints, u);
			System.out.println(" n " + n + " x " + values[n] + " check " + value);
		}
	}
}
