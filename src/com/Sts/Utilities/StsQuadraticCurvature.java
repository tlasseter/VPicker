package com.Sts.Utilities;

import Jama.*;


/** Computes Coefficients of Quadratic surface using Singular value decomposition */
public class StsQuadraticCurvature
{
    static public final float badCurvature =  99999999;
    static public final float curvatureTest = 9999999;
    
    static final private double chiSqrMultiplyer = 2;
    static final private double stdDevFactor = 2.5;


    static float points[][];
    static int nPoints;
    static Matrix A;
    static double[] B;
    static SingularValueDecomposition svd;
    static Matrix S;
    static Matrix V;
    static Matrix U;

    static final public int nCoefs = 6;
    static final public int minNPoints = 12;

    static double[] X = new double[nCoefs];
    static double[] s = new double[nCoefs];
    static double[] w = new double[nCoefs];

    //Curvature Attribute Types
    static public final byte CURVDip 	= 0;
	static public final byte CURVStrike = 1;
	static public final byte CURVMean 	= 2;
	static public final byte CURVGauss  = 3;
	static public final byte CURVMax	= 4;
	static public final byte CURVMin  	= 5;
	static public final byte CURVPos  	= 6;
	static public final byte CURVNeg  	= 7;

    public static final String CURVDipString = "Localized Dip";
    public static final String CURVStrikeString = "Localized Strike";
    public static final String CURVMeanString = "Mean Curvature";
    public static final String CURVGaussString = "Gaussian Curvature";
    public static final String CURVPosString = "Max Positive";
    public static final String CURVNegString = "Max Negative";
    public static final String CURVMinString = "Min Curvature";
    public static final String CURVMaxString = "Max Curvature";

    public static final String CURVDipName = "DIP_CURV";
    public static final String CURVStrikeName = "STRIKE_CURV";
    public static final String CURVMeanName = "MEAN_CURV";
    public static final String CURVGaussName = "GAUSS_CURV";
    public static final String CURVPosName = "MAX_POS_CURV";
    public static final String CURVNegName = "MIN_NEG_CURV";
    public static final String CURVMinName = "MIN_CURV";
    public static final String CURVMaxName = "MAX_CURV";

    /*  load A matrix with m x n values to fit
    *  ax2 + by2 + cxy + dx + ey + f = z
    *  A: [x2 y2 xy x y 1]
    *  b: [z]
    */

    static public boolean computeSVD(float fitPoints[][], int nFitPoints)
    {
        points = fitPoints;
        nPoints = nFitPoints;
        if(nPoints < nCoefs) return false;
        
        A = new Matrix(nPoints, nCoefs);
        B = new double[nPoints];
        for (int i = 0; i < nPoints; i++)
        {
            float x = points[i][0];
            float y = points[i][1];
            float z = points[i][2];
            A.set(i, 0, x * x);
            A.set(i, 1, y * y);
            A.set(i, 2, x * y);
            A.set(i, 3, x);
            A.set(i, 4, y);
            A.set(i, 5, 1);
            B[i] = z;
        }
        try
        {
            svd = A.svd();
            S = svd.getS();
            V = svd.getV();
            U = svd.getU();
            if(!svdOk()) return false;

            w = svd.getSingularValues();
            //double wMax = w[0];  //values are returned in largest -> smallest order
            //double cutOff = wMax * 1e-8;
            double cutOff = 1e-8;

            for (int j = 0; j < nCoefs; j++)
            {
                s[j] = 0;
                // Filter out low singular vals
                if (w[j] < cutOff) continue;

                for (int i = 0; i < nPoints; i++)
                {
                    s[j] += U.get(i, j) * B[i];
                }
                s[j] /= w[j];
            }
            // Do back substitution to get solution vector
            for (int j = 0; j < nCoefs; j++)
            {
                X[j] = 0;
                for (int jj = 0; jj < nCoefs; jj++)
                {
                    X[j] += V.get(j, jj) * s[jj];
                }
            }
            // coefficients are in X
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsQuadraticCurvature.class, "getQuadraticCoefficients", e);
            return false;
        }
    }

    static public double computeChiSquared()
    {
        double chiSqr = 0;
        for (int i = 0; i < nPoints; i++)
        {
            float x = points[i][0];
            float y = points[i][1];
            float z = points[i][2];
            double zCalc = X[0] * x * x + X[1] * y * y + X[2] * x * y + X[3] * x + X[4] * y + X[5];
            chiSqr += (z - zCalc) * (z - zCalc);
        }
        return chiSqr;
    }

    static private boolean svdOk()
    {
        Matrix Y = U.times(S.times(V.transpose()));
        return check(A, Y);
    }

    static private boolean check(Matrix X, Matrix Y)
    {
        double eps = Math.pow(2.0, -52.0);
        if (X.norm1() == 0. & Y.norm1() < 10 * eps) return true;
        if (Y.norm1() == 0. & X.norm1() < 10 * eps) return true;
        if (X.minus(Y).norm1() > 1000 * eps * Math.max(X.norm1(), Y.norm1()))
        {
            StsException.systemError("Singular value decomposition error: The norm of (X-Y) is too large: " + Double.toString(X.minus(Y).norm1()));
            return false;
        }
        return true;
    }

    static public double getNormXY()
    {
        Matrix Y = U.times(S.times(V.transpose()));
        return A.minus(Y).norm1();
    }

    static public float getCurvatureComponent(byte curveType)
    {
        double kDip = 0;
        double kStrike = 0;
        double kMean = 0;
        double kGauss = 0;

        double a = X[0];
        double b = X[1];
        double c = X[2];
        double d = X[3];
        double e = X[4];

        if (curveType == CURVDip)
        {
            double dipDenom = (e * e + d * d) * Math.pow((1 + d * d + e * e), 1.5);
            if (dipDenom == 0)                 //avoid divide by zero
            {
                if ((a > 0 && b > 0) || (a < 0 && b < 0))
                {
                    kDip = (a + b);            //summit or depression -(a+b) for Positive up
                }
                else
                {
                    kDip = 0;                    //saddle
                }
            }
            else
            {
                //assumes positive down
                kDip = 2 * (a * d * d + b * e * e + c * d * e) / dipDenom;
            }
        }
        else if (curveType == CURVStrike)
        {
            double strikeDenom = Math.pow((e * e + d * d), 1.5);
            if (strikeDenom == 0)                 //avoid divide by zero
            {
                if ((a > 0 && b > 0) || (a < 0 && b < 0))
                {
                    kStrike = StsParameters.nullValue;
                }
                else
                {
                    kStrike = 0;
                }
            }
            else
            {
                //assumes positive down
                kStrike = 2 * (b * d * d + a * e * e - c * d * e) / strikeDenom;
            }
        }
        if (curveType == CURVMean || curveType == CURVMin || curveType == CURVMax)
        {
            kMean = ((a * (1 + e * e) + b * (1 + d * d) - c * d * e)) / Math.pow((1 + d * d + e * e), 1.5);
        }
        if (curveType == CURVGauss || curveType == CURVMin || curveType == CURVMax)
        {
            kGauss = (4 * a * b - c * c) / Math.pow((1 + d * d + e * e), 2.0);
        }

        float val = 0.0f;
        switch (curveType)
        {
            case CURVDip:
                val = (float) kDip;
                break;
            case CURVStrike:
                val = (float) kStrike;
                break;
            case CURVMean:
                val = (float) kMean;
                break;
            case CURVGauss:
                val = (float) kGauss;
                break;
            case CURVPos:
                val = (float) ((a + b) + Math.sqrt((a - b) * (a - b) + c * c));
                break;
            case CURVNeg:
                val = (float) ((a + b) - Math.sqrt((a - b) * (a - b) + c * c));
                break;
            case CURVMin:
                val = (float) (kMean - Math.sqrt(kMean * kMean - kGauss));
                break;
            case CURVMax:
                val = (float) (kMean + Math.sqrt(kMean * kMean - kGauss));
                break;
            default:
                val = 0.0f;
                StsException.outputException(new StsException(StsException.
                    WARNING,
                    "SurfaceCurvTexture.", "Undefined type: " + curveType));

        }
        return val;
    }

    static public String coefficientsToString()
    {
        return StsToolkit.toString(X, "    coefficients: ");
    }


    static public String eigenValuesToString()
    {
        return StsToolkit.toString(w, "    eigenValues: ");
    }
}
