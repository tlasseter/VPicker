//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import com.Sts.Types.*;

import java.lang.reflect.*;
import java.nio.*;
import java.text.*;
import java.util.*;

public class StsMath
{
    //-------------------------
    // constants
    //-------------------------

    public final static double PI2 = Math.PI * 2.0;
    public final static double PI_HALF = Math.PI / 2.0;
    public final static double RADperDEG = Math.PI / 180.0;
    public final static double DEGperRAD = 1.0 / RADperDEG;
    public final static float roundOff = StsParameters.roundOff;
    public final static double doubleRoundOff = StsParameters.doubleRoundOff;
    public final static double doubleNullValue = StsParameters.doubleNullValue;

    //-------------------------
    // trig methods
    //-------------------------
    public static final double asinh(double x)
    {
        return Math.log(x + Math.sqrt(x * x + 1.0));
    }

    public static final double acosh(double x)
    {
        return Math.log(x + Math.sqrt(x * x - 1.0));
    }

    public static final double atanh(double x)
    {
        return Math.log((1.0 + x) / (1.0 - x)) / 2.0;
    }

    // given x and y return the angle between 0 and 360 degrees
    public static final double atan2(double x, double y)
    {
        double deg = DEGperRAD * Math.atan2(y, x);
        if (deg >= 0.0)
        {
            return deg;
        }
        else
        {
            return 360.0 + deg;
        }
    }

	static public double getAzimuthFromCCWEast(double ccwEast)
	{
		if(ccwEast < 0 || ccwEast > 360)
		{
			ccwEast = ccwEast%360;
			if(ccwEast < 0) ccwEast += 360;
		}
		double azimuth = (360 - ccwEast) + 90;
		if(azimuth > 360) azimuth -= 360;
		return azimuth;
	}

    public static final float atan2(float x, float y)
    {
        float deg = (float) (DEGperRAD * Math.atan2((double) y, (double) x));
        if (deg >= 0.0f)
        {
            return deg;
        }
        else
        {
            return 360.0f + deg;
        }
    }

    public static final float atan2(float[] xy)
    {
        return atan2(xy[0], xy[1]);
    }

    public static final double sind(double x)
    {
        return Math.sin(x * RADperDEG);
    }

    public static final double cosd(double x)
    {
        return Math.cos(x * RADperDEG);
    }

    public static final double tand(double x)
    {
        return Math.tan(x * RADperDEG);
    }

    public static final double asind(double x)
    {
        return Math.asin(x) * DEGperRAD;
    }

    public static final double acosd(double x)
    {
        return Math.acos(x) * DEGperRAD;
    }

    public static final double atand(double x)
    {
        return Math.atan(x) * DEGperRAD;
    }

    public static final double atan2d(double a, double b)
    {
        return Math.atan2(a, b) * DEGperRAD;
    }

    //-----------------------------------
    // truncation and fraction extraction
    //-----------------------------------

    public static final double trunc(double x)
    {
        if (x < 0.0)
        {
            return Math.ceil(x);
        }
        else
        {
            return Math.floor(x);
        }
    }

    public static final double frac(double x)
    {
        return x - trunc(x);
    }

    //-----------------------------------
    // Min or max routines
    //-----------------------------------

    public static final int minMax(int value, int min, int max)
    {
        if (value <= min) return min;
        else if (value >= max) return max;
        else return value;
    }

    public static float minMax(float value, float min, float max)
    {
        if (value <= min) return min;
        else if (value >= max) return max;
        else return value;
    }

    public static final double minMax(double value, double min, double max)
    {
        if (value <= min) return min;
        else if (value >= max) return max;
        else return value;
    }

    public static final long minMax(long value, long min, long max)
    {
        if (value <= min) return min;
        else if (value >= max) return max;
        else return value;
    }

    public static final int min2(int v1, int v2)
    {
        if (v1 < v2)
        {
            return v1;
        }
        else
        {
            return v2;
        }
    }

    public static final short min2(short v1, short v2)
    {
        if (v1 < v2)
        {
            return v1;
        }
        else
        {
            return v2;
        }
    }

    public static final int max3(int v1, int v2, int v3)
    {
        return (v1 > v2 ? (v1 > v3 ? v1 : v3) : (v2 > v3 ? v2 : v3));
    }

    public static final int min3(int v1, int v2, int v3)
    {
        return (v1 < v2 ? (v1 < v3 ? v1 : v3) : (v2 < v3 ? v2 : v3));
    }

    public static final int max4(int v1, int v2, int v3, int v4)
    {
        return (v1 > v2 ? max3(v1, v3, v4) : max3(v2, v3, v4));
    }

    public static final int min4(int v1, int v2, int v3, int v4)
    {
        return (v1 < v2 ? min3(v1, v3, v4) : min3(v2, v3, v4));
    }

    public static final float max3(float v1, float v2, float v3)
    {
        return (v1 > v2 ? (v1 > v3 ? v1 : v3) : (v2 > v3 ? v2 : v3));
    }

    public static final float min3(float v1, float v2, float v3)
    {
        return (v1 < v2 ? (v1 < v3 ? v1 : v3) : (v2 < v3 ? v2 : v3));
    }

    public static final float max4(float v1, float v2, float v3, float v4)
    {
        return (v1 > v2 ? max3(v1, v3, v4) : max3(v2, v3, v4));
    }

    public static final float min4(float v1, float v2, float v3, float v4)
    {
        return (v1 < v2 ? min3(v1, v3, v4) : min3(v2, v3, v4));
    }

    public static final float min(float[] values)
    {
        if (values == null)
        {
            return StsParameters.nullValue;
        }
        int nValues = values.length;
        float min = values[0];
        for (int n = 1; n < nValues; n++)
        {
            if (values[n] < min)
            {
                min = values[n];
            }
        }
        return min;
    }

    public static final float minExcludeNull(float[] values)
    {
        if (values == null)
        {
            return StsParameters.nullValue;
        }
        int nValues = values.length;
        float min = StsParameters.largeFloat;
        for (int n = 0; n < nValues; n++)
        {
            if (values[n] == StsParameters.nullValue) continue;
            if (values[n] < min)
            {
                min = values[n];
            }
        }
        return min;
    }

    public static final float max(float[] values)
    {
        if (values == null || values.length == 0) return 0.0f;
        int nValues = values.length;
        float max = values[0];
        for (int n = 1; n < nValues; n++)
        {
            if (values[n] > max)
            {
                max = values[n];
            }
        }
        return max;
    }

    public static final double max(double[] values)
    {
        if (values == null || values.length == 0) return 0.0;
        int nValues = values.length;
        double max = values[0];
        for (int n = 1; n < nValues; n++)
        {
            if (values[n] > max)
            {
                max = values[n];
            }
        }
        return max;
    }

    public static final float maxExcludeNull(float[] values)
    {
        if (values == null)
        {
            return StsParameters.nullValue;
        }
        int nValues = values.length;
        float max = -StsParameters.largeFloat;
        for (int n = 0; n < nValues; n++)
        {
            if (values[n] == StsParameters.nullValue) continue;
            if (values[n] > max)
            {
                max = values[n];
            }
        }
        return max;
    }

    public static double mean(double[] vals)
    {
        if (vals == null) return 0.0;
        double mean = 0.0;
        for (int i = 0; i < vals.length; i++)
            mean += vals[i];
        return mean / vals.length;
    }

    /**
     * calculates mean of array "vals", but ignores zeroed values
     *
     * @param vals
     * @return
     */
    public static double meanIgnoreZeroes(double[] vals)
    {
        if (vals == null) return 0.0;
        int totalVals = 0;
        double sum = 0.0;
        for (int i = 0; i < vals.length; i++)
        {
            if (vals[i] == 0.0) continue;
            totalVals++;
            sum += vals[i];
        }
        if (totalVals > 0)
            return sum / totalVals;
        else
            return 0.0;
    }

    public static final float[] range(float[] values)
    {
        if (values == null) return null;
        int nValues = values.length;
        float min = values[0];
        float max = values[0];
        for (int n = 1; n < nValues; n++)
        {
            if (values[n] > max)
                max = values[n];
            else if (values[n] < min)
                min = values[n];
        }
        return new float[]{min, max};
    }

    public static final double[] range(double[] values)
    {
        if (values == null) return null;
        int nValues = values.length;
        double min = values[0];
        double max = values[0];
        for (int n = 1; n < nValues; n++)
        {
            if (values[n] > max)
                max = values[n];
            else if (values[n] < min)
                min = values[n];
        }
        return new double[]{min, max};
    }

    public static final boolean betweenInclusive(float value, float limit1, float limit2)
    {
        return limit1 <= value && value <= limit2 || limit2 <= value && value <= limit1;
    }

    public static final boolean betweenExclusive(float value, float limit1, float limit2)
    {
        return limit1 < value && value < limit2 || limit2 < value && value < limit1;
    }

    static public int BELOW = -2;
    static public int ON_BELOW = -1;
    static public int BETWEEN = 0;
    static public int ON_ABOVE = 1;
    static public int ABOVE = 2;

    public static final int belowBetweenAbove(float value, float limit1, float limit2)
    {
        if (value < limit1) return BELOW;
        else if (value == limit1) return ON_BELOW;
        else if (value == limit2) return ON_ABOVE;
        else if (value > limit2) return ABOVE;
        else return BETWEEN;
    }

    public static final boolean betweenInclusive(double value, double limit1, double limit2)
    {
        return limit1 <= value && value <= limit2 || limit2 <= value && value <= limit1;
    }

    public static final boolean betweenExclusive(double value, double limit1, double limit2)
    {
        return limit1 < value && value < limit2 || limit2 < value && value < limit1;
    }

    public static final float remainder(float value, float divisor)
    {
        if (divisor == 0.0f)
        {
            return 0.0f;
        }

        int n = (int) (value / divisor);
        float r = value - n * divisor;
        if (r < 0.0f)
        {
            r += divisor;
        }
        return r;
    }

    public static final int sign(int number)
    {
        if (number > 0)
        {
            return 1;
        }
        else if (number < 0)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }

    public static final int sign(float number)
    {
        if (number > 0.0f)
        {
            return 1;
        }
        else if (number < 0.0f)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }

    public static final int sign(double number)
    {
        if (number > 0)
        {
            return 1;
        }
        else if (number < 0)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }

    public static final int signRoundOff(float number)
    {
        if (number > roundOff)
        {
            return 1;
        }
        else if (number < -roundOff)
        {
            return -1;
        }
        else
        {
            return 0;
        }
    }

    public static final boolean sameAsOld(double value1, double value2)
    {
        double diff = Math.abs(value1 - value2);
        if (diff == 0.0)
        {
            return true;
        }
        double largest = Math.max(Math.abs(value1), Math.abs(value2));
        double scaledRoundOff = largest * doubleRoundOff;
        return diff <= scaledRoundOff;
    }

    /**
     * The smallest positive number that can be added to a float 1
     * before the number is considered different from 1.
     */
    public static final float FLT_EPSILON = 1.192092896e-07f;

    /** Default epsilon assumes math has destroyed one digit of accuracy. */
    private static final double fEpsilon = 10. * FLT_EPSILON;

    /**
     * The absolutely smallest positive number that can be added to a double 1
     * before the number is considered different from 1.
     */
    public static final double DBL_EPSILON = 2.2204460492503131e-016;

    /** Default epsilon assumes math has destroyed one digit of accuracy. */
    private static final double dEpsilon = 10. * DBL_EPSILON;

    /**
     * This default minValue may be too small.
     * If subtraction is involved, to same value as epsilon
     */
    private static final double fMinValue = 100. * Float.MIN_VALUE;

    /**
     * This default minValue may be too small.
     * If subtraction is involved, to same value as epsilon
     */
    private static final double dMinValue = 100. * Double.MIN_VALUE;

    public static final boolean sameAs(float f1, float f2)
    {
        return sameAsRoundoff(f1, f2, fEpsilon);
    }

    public static final boolean sameAsRoundoff(float f1, float f2, float roundOff)
    {
        return sameAsRoundoff(f1, f2, (double)roundOff);
    }

    public static final boolean sameAsRoundoff(float f1, float f2, double roundOff)
    {
        double v1 = f1;
        double v2 = f2;
        if (v1 == v2) return true; // no need for more careful test

        double av1, av2; // absolute values
        if (v1 < 0.) av1 = -v1;
        else av1 = v1;
        if (v2 < 0.) av2 = -v2;
        else av2 = v2;
        // easy heuristic check for values close to zero
        if (av1 < fMinValue && av2 < fMinValue) return true;

        double ev1 = roundOff * av1;
        double ev2 = roundOff * av2;
        if (v1 - ev1 > v2 + ev2) return false;
        else if (v1 + ev1 < v2 - ev2) return false;
        else return true;
    }

    public static final boolean sameAs(float[] f1, float[] f2)
    {
        if (f1 == null || f2 == null) return false;
        if (f1.length != f2.length) return false;
        for (int n = 0; n < f1.length; n++)
            if (f1[n] != f2[n]) return false;
        return true;
    }

    public static final boolean sameAs(double v1, double v2)
    {
        if (v1 == v2) return true; // no need for more careful test

        double av1, av2; // absolute values
        if (v1 < 0.) av1 = -v1;
        else av1 = v1;
        if (v2 < 0.) av2 = -v2;
        else av2 = v2;
        // easy heuristic check for values close to zero
        if (av1 < dMinValue && av2 < dMinValue) return true;

        double ev1 = dEpsilon * av1;
        double ev2 = dEpsilon * av2;
        if (v1 - ev1 > v2 + ev2) return false;
        else if (v1 + ev1 < v2 - ev2) return false;
        else return true;
    }

    public static final boolean sameAsOld(float value1, float value2)
    {
        float diff = Math.abs(value1 - value2);
        if (diff == 0.0f)
        {
            return true;
        }
        float largest = Math.max(Math.abs(value1), Math.abs(value2));
        float scaledRoundOff = largest * roundOff;
        return diff <= scaledRoundOff;
    }

    public static final boolean sameAsTol(float value1, float value2, float tol)
    {
        return Math.abs(value1 - value2) <= tol;
    }

    public static final boolean sameAs(double value1, double value2, double tol)
    {
        return Math.abs(value1 - value2) <= tol;
    }

    public static final boolean near(int value1, int value2, int tol)
    {
        if (value1 == value2)
        {
            return true;
        }
        int dif = value1 - value2;
        if (dif > 0)
        {
            return dif <= tol;
        }
        else
        {
            return -dif <= tol;
        }
    }

    public static final boolean isIntegral(float value)
    {
        if (value == 0.0f)
        {
            return true;
        }
        int ivalue = (int) value;
        float remainder = Math.abs(value - (float) ivalue);
        float scaledRoundOff = Math.abs(value) * roundOff;
        return remainder <= scaledRoundOff || 1.0f - remainder <= scaledRoundOff;
    }

    public static final boolean isIntegral(float value, float tol)
    {
        if (value == 0.0f)
        {
            return true;
        }
        int ivalue = (int) value;
        float remainder = Math.abs(value - (float) ivalue);
        return remainder <= tol || 1.0f - remainder <= tol;
    }

    public static final boolean isExactIntegral(float value)
    {

        if (value == 0.0f)
        {
            return true;
        }
        int ivalue = (int) value;
        return value == (float) ivalue;
    }

    public static final boolean isIntegralRatio(float value1, float value2, float scaledError)
    {
        float ratio;

        if (value1 == 0.0f || value2 == 0.0f)
        {
            return false;
        }
        if (sameAs(value1, value2))
        {
            return true;
        }
        if (value1 > value2)
        {
            ratio = value1 / value2;
        }
        else
        {
            ratio = value2 / value1;
        }
        return Math.abs(ratio - Math.round(ratio)) <= scaledError;
    }

    public static final float checkRoundOffInteger(float value)
    {
        return checkRoundOffIntegerAdjust(value, roundOff);
    }

    public static final float checkRoundOffIntegerAdjust(float value, float roundOff)
    {
        if (value == 0.0f)
        {
            return value;
        }
        float scaledRoundOff = roundOff * Math.abs(value);
        int ivalue = (int) value;
        float remainder = Math.abs(value - (float) ivalue);
        if (remainder <= scaledRoundOff)
        {
            return (float) ivalue;
        }
        else if (1.0f - remainder <= scaledRoundOff)
        {
            return (float) (ivalue + 1);
        }
        else
        {
            return value;
        }
    }

    public static final int roundOffInteger(float value)
    {
        int ivalue = (int) value;
        float remainder = value - (float) ivalue;
        if (remainder <= 0.5f)
        {
            return ivalue;
        }
        else
        {
            return (ivalue + 1);
        }
    }

    public static final int roundOffInteger(double value)
    {
        int ivalue = (int) value;
        double remainder = value - (double) ivalue;
        if (remainder <= 0.5)
        {
            return ivalue;
        }
        else
        {
            return (ivalue + 1);
        }
    }

    /** Return integral value >= value */
    public static final int ceiling(float value)
    {
        int ivalue = (int) value;
        if ((float) ivalue == value)
        {
            return ivalue;
        }
        else
        {
            if (value > 0.0f)
            {
                return ivalue + 1;
            }
            else
            {
                return ivalue;
            }
        }
    }

    // if value is within roundOff of an integer value, return it;
    // otherwise return one below
    public static final int ceilingRoundOff(float value)
    {
        if (value == 0.0f)
        {
            return (int) value;
        }
        if (value < 0.0f)
        {
            return -floorRoundOff(-value);
        }

        int ivalue = (int) value;
        float remainder = value - (float) ivalue;
        if (remainder <= roundOff)
        {
            return ivalue;
        }
        else
        {
            return ivalue + 1;
        }
    }

    public static final int ceiling(double value)
    {
        int ivalue = (int) value;
        if ((double) ivalue == value)
        {
            return ivalue;
        }
        else
        {
            if (value > 0.0)
            {
                return ivalue + 1;
            }
            else
            {
                return ivalue;
            }
        }
    }

    /** Return integral value > value */
    public static final int above(float value)
    {
        if (value >= 0.0f)
        {
            int ivalue = (int) value;
            return ivalue + 1;
        }
        else
        {
            int ivalue = (int) value;
            return ivalue;
        }
    }

    public static final int above(double value)
    {
        return above((float) value);
    }

    public static final int aboveRoundoff(float value)
    {
        if (value >= 0.0f)
        {
            int ivalue = (int) (value + roundOff);
            return ivalue + 1;
        }
        else
        {
            int ivalue = (int) (value + roundOff);
            return ivalue;
        }
    }

    /** Return float of integral value >= value within roundOff */
    public static final float ceilingF(float value)
    {
        int i = (int) value;
        float iF = (float) i;
        if (iF == value)
        {
            return iF;
        }
        else
        {
            return iF + 1.0f;
        }
    }

    public static final float ceilingFRoundoff(float value)
    {
        int i = (int) (value - roundOff);
        float iF = (float) i;
        if (iF == value)
        {
            return iF;
        }
        else
        {
            return iF + 1.0f;
        }
    }

    /** Return integral value <= value */
    public static final int floor(float value)
    {
        if (value >= 0.0f)
        {
            return (int) value;
        }
        else
        {
            int ivalue = (int) value;
            if ((float) ivalue == value)
            {
                return ivalue;
            }
            else
            {
                return ivalue - 1;
            }
        }
    }

    public static final int floor(double value)
    {
        if (value >= 0.0)
        {
            return (int) value;
        }
        else
        {
            int ivalue = (int) value;
            if ((double) ivalue == value)
            {
                return ivalue;
            }
            else
            {
                return ivalue - 1;
            }
        }
    }

    // if value is withing roundOff of an integer value, return it;
    // otherwise return one below
    public static final int floorRoundOff(float value)
    {
        if (value == 0.0f)
        {
            return (int) value;
        }
        if (value < 0.0f)
        {
            return -ceilingRoundOff(-value);
        }

        int ivalue = (int) value;
        float remainder = value - (float) ivalue;
        if (1.0f - remainder <= roundOff)
        {
            return ivalue + 1;
        }
        else
        {
            return ivalue;
        }
    }

    public static final int below(double value)
    {
        return below((float) value);
    }

    /** Return integral value < value */
    public static final int below(float value)
    {
        int ivalue = (int) value;
        if ((float) ivalue == value)
        {
            return ivalue - 1;
        }
        else
        {
            return ivalue;
        }
    }

    /** Return float of integral value <= value */
    public static final float floorF(float value)
    {
        int i = (int) value;
        return (float) i;
    }

    public static final float intervalRoundUp(float value, float interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundDown(value, -interval);
        }
        return interval * (StsMath.ceiling(value / interval));
    }

    public static final double intervalRoundUp(double value, double interval)
    {
        if (interval == 0.0)
        {
            return value;
        }
        else if (interval < 0.0)
        {
            return intervalRoundDown(value, -interval);
        }
        return interval * (StsMath.ceiling(value / interval));
    }

    public static final float intervalRoundDown(float value, float interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundUp(value, -interval);
        }
        return interval * (StsMath.floor(value / interval));
    }

    public static final double intervalRoundDown(double value, double interval)
    {
        if (interval == 0.0)
        {
            return value;
        }
        else if (interval < 0.0)
        {
            return intervalRoundUp(value, -interval);
        }
        return interval * (StsMath.floor(value / interval));
    }

    public static final int intervalRound(int value, int interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRound(value, -interval);
        }
        int remainder = value % interval;
        if (remainder == 0)
            return value;
        else if (remainder <= interval / 2)
            return value - remainder;
        else
            return value - remainder + interval;
    }

    public static final int intervalRound(int value, int min, int interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRound(value, min, -interval);
        }
        return min + interval * Math.round((float) (value - min) / interval);
    }

    public static final float intervalRound(float value, float min, float interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRound(value, min, -interval);
        }
        return min + interval * Math.round((float) (value - min) / interval);
    }

    public static final int intervalRoundUp(int value, int min, int interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRoundDown(value, min, -interval);
        }
        return min + interval * (StsMath.ceiling((float) (value - min) / interval));
    }

    public static final int intervalRoundDown(int value, int min, int interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRoundUp(value, min, -interval);
        }
        return min + interval * (StsMath.floor((float) (value - min) / interval));
    }

    public static final int intervalRoundUp(int value, int interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundDown(value, -interval);
        }
        int remainder = value % interval;
        if (remainder == 0)
        {
            return value;
        }
        else
        {
            return value - remainder + interval;
        }
    }

    public static final int intervalRoundDown(int value, int interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRoundUp(value, -interval);
        }
        return interval * (value / interval);
    }

    public static final float intervalRoundUp(float value, float min, float interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundDown(value, min, -interval);
        }
        return min + interval * (StsMath.ceiling((value - min) / interval));
    }

    public static final double intervalRoundUp(double value, double min, double interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundDown(value, min, -interval);
        }
        return min + interval * (StsMath.ceiling((value - min) / interval));
    }

    public static final float intervalRoundDown(float value, float min, float interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundUp(value, min, -interval);
        }
        return min + interval * (StsMath.floor((value - min) / interval));
    }

    public static final double intervalRoundDown(double value, double min, double interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundUp(value, min, -interval);
        }
        return min + interval * (StsMath.floor((value - min) / interval));
    }

    /** given a value and interval, round to the next integral value (evenly divisible by interval) which is above ( > value ).
     *  Thus if value is integral, it will return value + interval (the next one above).
     */
    public static final int intervalRoundUpAbove(int value, int interval)
    {
        if (interval == 0.0f)
        {
            return value;
        }
        else if (interval < 0.0f)
        {
            return intervalRoundDown(value, -interval);
        }
        int remainder = value % interval;
        if (remainder == 0)
        {
            return value + interval;
        }
        else
        {
            return value - remainder + interval;
        }
    }

    /** given a value and interval, round to the next integral value (evenly divisible by interval) which is below ( < value ).
     *  Thus if value is integral, it will return value - interval (the next one below).
     */
    public static final int intervalRoundDownBelow(int value, int interval)
    {
        if (interval == 0)
        {
            return value;
        }
        else if (interval < 0)
        {
            return intervalRoundUpAbove(value, -interval);
        }
        int remainder = value % interval;
        if (remainder == 0)
        {
            return value - interval;
        }
        else
        {
            return value - remainder - interval;
        }
    }

    public static final boolean equals(int i, float f)
    {
        int fi = (int) f;
        return (float) i == f;
    }

    public static final boolean equals(float f, int i)
    {
        int fi = (int) f;
        return (float) i == f;
    }

    public static final boolean equals(float[] f1, float[] f2)
    {
        int length = Math.min(f1.length, f2.length);
        for (int n = 0; n < length; n++)
        {
            if (f1[n] != f2[n])
            {
                return false;
            }
        }
        return true;
    }

    public static final boolean isEven(int number)
    {
        return number % 2 == 0;
    }

    /**
     * Java Arrays.binarySearch is busted!  Use this one.
     * returns index of element in vector "a" just below value "x".
     */
    static public int binarySearch(long[] a, long x)
    {
        int low = 0;
        int high = a.length - 1;
        while (low <= high)
        {
            int mid = (low + high) / 2;
            if (a[mid] < x) low = mid + 1;
            else if (a[mid] > x) high = mid - 1;
            else return mid;
        }
        return low;
    }

    /**
     * Java Arrays.binarySearch is busted!  Use this one.
     * returns index of element in vector "a" just below value "x".
     */
    static public int binarySearch(double[] a, double x)
    {
        int low = 0;
        int high = a.length - 1;
        while (low <= high)
        {
            int mid = (low + high) / 2;
            if (a[mid] < x) low = mid + 1;
            else if (a[mid] > x) high = mid - 1;
            else return mid;
        }
        return low;
    }

    /**
     * Java Arrays.binarySearch is busted!  Use this one.
     * returns index of element in vector "a" just below value "x".
     */
    static public int binarySearch(int[] a, int x)
    {
        int low = 0;
        int high = a.length - 1;
        while (low <= high)
        {
            int mid = (low + high) / 2;
            if (a[mid] < x) low = mid + 1;
            else if (a[mid] > x) high = mid - 1;
            else return mid;
        }
        return low;
    }

    /** returns interpolating index of x in vector a */
    static public float binarySearchF(int[] a, int x)
    {
        int low = 0;
        int high = a.length - 1;
        while (low <= high)
        {
            int mid = (low + high) / 2;
            if (a[mid] < x) low = mid + 1;
            else if (a[mid] > x) high = mid - 1;
            else
            {
                return mid + (x - a[mid]) / (a[mid + 1] - a[mid]);
            }
        }
        return low;
    }

    /** returns interpolating index of x in vector a */
    static public float binarySearchF(double[] a, double x)
    {
        int low = 0;
        int high = a.length - 1;
        while (low <= high)
        {
            int mid = (low + high) / 2;
            if (a[mid] < x) low = mid + 1;
            else if (a[mid] > x) high = mid - 1;
            else
            {
                return (float) (mid + (x - a[mid]) / (a[mid + 1] - a[mid]));
            }
        }
        return low;
    }

    /** returns interpolating index of object in monotonically increasing/decreasing array of Comparable objects just before Comparable key object */

    static public int binarySearchObject(Object[] a, Object key)
    {
        int low = 0;
        int high = a.length - 1;

        while (low <= high)
        {
            int mid = (low + high) / 2;
            Comparable midVal = (Comparable) a[mid];
            int cmp = midVal.compareTo(key);

            if (cmp < 0)
                low = mid + 1;
            else if (cmp > 0)
                high = mid - 1;
            else
                return mid; // key found
        }
        return low;
    }

    //-----------------------------------
    // Geometry operations
    //-----------------------------------

    public static final float distance(StsPoint point1, StsPoint point2)
    {
        int length = Math.min(point1.v.length, point2.v.length);
        float dv, sum = 0.0f;

        for (int n = 0; n < length; n++)
        {
            dv = point1.v[n] - point2.v[n];
            sum += dv * dv;
        }

        return (float) Math.sqrt(sum);
    }

    public static final float distanceSq(float[] p1, float[] p2, int length)
    {
        double dv, sum = 0.0;

        for (int n = 0; n < length; n++)
        {
            dv = p1[n] - p2[n];
            sum += dv * dv;
        }
        return (float) sum;
    }

    public static final double distanceSq(double[] p1, double[] p2, int length)
    {
        double dv, sum = 0.0;

        for (int n = 0; n < length; n++)
        {
            dv = p1[n] - p2[n];
            sum += dv * dv;
        }
        return sum;
    }

    public static final float distance(float[] p1, float[] p2, int length)
    {
        return (float) Math.sqrt(distanceSq(p1, p2, length));
    }

    public static final double distance(double[] p1, double[] p2, int length)
    {
        return Math.sqrt(distanceSq(p1, p2, length));
    }

    public static final float length(float[] v)
    {
        int dim = v.length;

        float dv, sum = 0.0f;

        for (int n = 0; n < dim; n++)
        {
            dv = v[n];
            sum += dv * dv;
        }
        return (float) Math.sqrt(sum);
    }

    public static final double length(double[] v)
    {
        int dim = v.length;

        double dv, sum = 0.0f;

        for (int n = 0; n < dim; n++)
        {
            dv = v[n];
            sum += dv * dv;
        }
        return Math.sqrt(sum);
    }

    public static final float lengthSq(float[] v)
    {
        int dim = v.length;

        float dv, sum = 0.0f;

        for (int n = 0; n < dim; n++)
        {
            dv = v[n];
            sum += dv * dv;
        }
        return sum;
    }

    /**
     * Given an interval 0 to 1, slopes and values at the end points, compute
     * the interval location for a maximum of this curve.
     *
     * @param v0 value at 0
     * @param v1 value at 1
     * @param s0 slope at 0
     * @param s1 slope at 1
     * @return null if no max, otherwise return max in array of length 1
     */
    public static final double findCubicMax(double v0, double v1, double s0, double s1, boolean debug)
    {
        double a = 0, b = 0, c = 0, d = 0;
        double[] roots = null;

        if (debug)
        {
            System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);

            // handle simple cases when one or both slopes are 0.0
            // assuming curve is parabolic

        }
        if (s0 == 0.0f)
        {
            if (s1 == 0.0)
            {
                if (v1 < v0)
                {
                    return v0;
                }
                else if (v1 > v0)
                {
                    return v1;
                }
                else
                {
                    return doubleNullValue;
                }
            }
            else if (s1 < 0.0)
            {
                if (v1 < v0)
                {
                    return v0;
                }
            }
            else // s1 > 0.0f
            {
                if (v1 <= v0)
                {
                    return v0;
                }
                else
                {
                    return doubleNullValue;
                }
            }
        }
        else if (s1 == 0.0)
        {
            if (s0 < 0.0)
            {
                if (v1 >= v0)
                {
                    return v1;
                }
            }
            else // s0 > 0.0f; s0 == 0.0f handled above
            {
                if (v0 < v1)
                {
                    return v1;
                }
            }
        }

        // for other cases handle as a cubic

        a = 2 * (v0 - v1) + s0 + s1;
        b = 3 * (v1 - v0) - 2 * s0 - s1;
        c = s0;
        roots = new double[2];
        int nRoots = findQuadraticRoots(3 * a, 2 * b, c, roots);
        for (int n = 0; n < nRoots; n++)
        {
            if (StsMath.betweenInclusive(roots[n], 0.0, 1.0) && 6 * a * roots[n] + 2 * b < 0)
            {
                double x = roots[n];
                double y = (x * (x * (a * x + b) + c) + d);
                if (!debug) return y;

                // mainDebug output
                d = v0;
                System.out.print(" a: " + a + " b: " + b + " c: " + c + " d: " + d);
                double yp = 3 * a * x * x + 2 * b * x + c;
                double ypp = 6 * a * x + 2 * b;
                System.out.println(" Max at: " + x + " y: " + (float) y + " y': " + (float) yp + " y'': " + (float) ypp);
                return y;
            }
        }
        return doubleNullValue;
    }

    public static final int findQuadraticRoots(double a, double b, double c, double[] roots)
    {
        int nRoots = 0;
        double d, q;

        d = b * b - 4 * a * c;
        if (d < 0)
        {
            return 0;
        }

        if (b >= 0)
        {
            q = -(b + Math.sqrt(d)) / 2;
        }
        else
        {
            q = -(b - Math.sqrt(d)) / 2;

        }
        if (a != 0)
        {
            roots[nRoots++] = q / a;
        }
        if (q != 0)
        {
            roots[nRoots++] = c / q;
        }
        return nRoots;
    }

    public static final double findQuadraticRootBetween(double a, double b, double c, double min, double max)
    {
        int nRoots = 0;
        double d, q;

        d = b * b - 4 * a * c;
        if (d < 0)
        {
            return 0;
        }

        if (b >= 0)
        {
            q = -(b + Math.sqrt(d)) / 2;
        }
        else
        {
            q = -(b - Math.sqrt(d)) / 2;

        }
        if (a != 0)
        {
            double root = q / a;
            if (StsMath.betweenInclusive(root, min, max)) return root;
        }
        if (q != 0)
        {
            double root = c / q;
            if (StsMath.betweenInclusive(root, min, max)) return root;
        }
        return doubleNullValue;
    }

    /**
     * Given quadratic y = ax**2 + bx + c, y(-1) = v0, y(0) = v1, y(1) = v2 and
     * knowing that a max or min exists at y'(?) = 0, return ?
     */
    public static final float findQuadraticMinMax(float v0, float v1, float v2, boolean debug)
    {
        if (debug)
        {
            System.out.print("v0: " + v0 + " v1: " + v1 + " v2: " + v2);
        }
        float max = (v0 - v2) / (v0 + v2 - 2 * v1) / 2;
        if (debug)
        {
            System.out.println(" max: " + max);
        }
        return max;
    }

    public static final float hermiteCubic(float v0, float v1, float s0, float s1, float f)
    {
        double a = 0, b = 0, c = 0, d = 0;

        a = 2 * (v0 - v1) + s0 + s1;
        b = 3 * (v1 - v0) - 2 * s0 - s1;
        c = s0;
        d = v0;
        return (float) (f * (f * (a * f + b) + c) + d);
    }

    public static final double hermiteCubic(double v0, double v1, double s0, double s1, double f)
    {
        double a = 0, b = 0, c = 0, d = 0;

        a = 2 * (v0 - v1) + s0 + s1;
        b = 3 * (v1 - v0) - 2 * s0 - s1;
        c = s0;
        d = v0;
        return (f * (f * (a * f + b) + c) + d);
    }

    public static final double[] hermiteValueAndSlope(double z, double v0, double v1, double s0, double s1, double f)
    {
        if (f == 0.0)
            return new double[]{z, v0, s0};
        else if (f == 1.0)
            return new double[]{z, v1, s1};
        else
        {
            double a = 2 * (v0 - v1) + s0 + s1;
            double b = 3 * (v1 - v0) - 2 * s0 - s1;
            double c = s0;
            double d = v0;
            double value = (f * (f * (a * f + b) + c) + d);
            double slope = (f * (3 * a * f + 2 * b) + c);
            return new double[]{z, value, slope};
        }
    }

    public static final float cubic(float v0, float v1, float v0m, float v1p, float f)
    {
        double a, b, c, d;
        double s0 = (v1 - v0m) / 2;
        double s1 = (v1p - v0) / 2;
        a = 2 * (v0 - v1) + s0 + s1;
        b = 3 * (v1 - v0) - 2 * s0 - s1;
        c = s0;
        d = v0;
        return (float) (f * (f * (a * f + b) + c) + d);
    }

    public static final double cubic(double v0, double v1, double v0m, double v1p, double f)
    {
        double a, b, c, d;
        double s0 = (v1 - v0m) / 2;
        double s1 = (v1p - v0) / 2;
        a = 2 * (v0 - v1) + s0 + s1;
        b = 3 * (v1 - v0) - 2 * s0 - s1;
        c = s0;
        d = v0;
        return (f * (f * (a * f + b) + c) + d);
    }

    /**
     * find the solution to sum{a[n]x^n} for n = 0,d where d is degree of
     * polynomial using Newton-Raphson plus bisection with the solution between x1
     * and x2. Return the number of real roots found and put roots in array roots.
     *
     * @param a     coefficients of polynomial with a[0] being the constant
     * @param x1    minimum bracket value
     * @param x2    maximum bracket value
     * @param error max error allowed in convergence
     * @return root value if found or doubleNullValue if not.
     */
    public static final double findRoot(double[] a, double x1, double x2, double error, boolean debug)
    {
        int maxIterations = 100;
        double[] e = evalPoly(a, 1, x1);
        double f1 = e[0];
        e = evalPoly(a, 1, x2);
        double f2 = e[0];
        if (f1 > 0.0 && f2 > 0.0 || f1 < 0.0 && f2 < 0.0)
        {
            if (debug) System.out.println("Root doesn't lie between bracketed values.");
            return doubleNullValue;
        }
        if (f1 == 0.0)
        {
            if (debug) System.out.println("Value is lower end of bracket: zero iterations.");
            return x1;
        }
        else if (f2 == 0.0)
        {
            if (debug) System.out.println("Value is upper end of bracket: zero iterations.");
            return x2;
        }
        double xl, xh;
        if (f1 <= 0.0)
        {
            xl = x1;
            xh = x2;
        }
        else
        {
            xl = x2;
            xh = x1;
        }
        double rts = 0.5 * (x1 + x2); // classInitialize guess of root
        double dxold = Math.abs(x2 - x1);
        double dx = dxold;
        e = evalPoly(a, 1, rts);
        double f = e[0];
        double df = e[1];
        for (int j = 1; j <= maxIterations; j++)
        {
            if ((((rts - xh) * df - f) * ((rts - xl) * df - f) >= 0.0) || (Math.abs(2 * f) > Math.abs(dxold * df)))
            {
                dxold = dx;
                dx = 0.5 * (xh - xl);
                rts = xl + dx;
                if (xl == rts)
                {
                    if (debug) System.out.println("Converged in " + j + " iterations.");
                    return rts;
                }
            }
            else
            {
                dxold = dx;
                dx = f / df;
                double temp = rts;
                rts -= dx;
                if (temp == rts)
                {
                    if (debug) System.out.println("Converged in " + j + " iterations.");
                    return rts;
                }
            }

            if (Math.abs(dx) < error)
            {
                if (debug) System.out.println("Converged in " + j + " iterations.");
                return rts;
            }
            e = evalPoly(a, 1, rts);
            f = e[0];
            df = e[1];
            if (f < 0.0)
                xl = rts;
            else
                xh = rts;
        }
        if (debug) System.out.println("Failed to converge in " + maxIterations + " iterations.");
        return doubleNullValue;
    }

    /**
     * evalPoly computes the value and a specified number of derivates of a
     * polynomial. From "Numerical Recipes", page 174.
     *
     * @param a  coefficients of polynomial with a[0] being the constant
     * @param nd number of derivatives desired
     * @param x  value at which polynomial is computed
     * @return array with first entry being the value of polynomial followed by the number of derivatives desired.
     */
    static public final double[] evalPoly(double[] a, int nd, double x)
    {
        int nc = a.length - 1; // degree of polynomial
        double[] pd = new double[1 + nd];
        pd[0] = a[nc];
        for (int i = 1; i <= nd; i++)
            pd[i] = 0.0;
        for (int i = nc - 1; i >= 0; i--)
        {
            int nnd = (nd < (nc - i) ? nd : nc - i);
            for (int j = nnd; j >= 1; j--)
                pd[j] = pd[j] * x + pd[j - 1];
            pd[0] = pd[0] * x + a[i];
        }
        double c = 1.0;
        // after first derivate, factorial constants come in
        for (int i = 2; i <= nd; i++)
        {
            c *= i;
            pd[i] *= c;
        }
        return pd;
    }

    //-----------------------------------
    // Array operations
    //-----------------------------------

    public static final int arrayLength(Object array)
    {
        if (array == null)
        {
            return 0;
        }
        Class c = array.getClass();
        if (!c.isArray())
        {
            return 0;
        }
        return Array.getLength(array);
    }

    public static Object arrayGrow(Object array, float growFraction)
    {
        Class componentType = componentType(array);
        int length = Array.getLength(array);
        int newLength = (int) (length * (1.0f + growFraction));
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, 0, newArray, 0, length);
        return newArray;
    }

    static private Class componentType(Object array)
    {
        if (array == null)
        {
            return null;
        }
        return array.getClass().getComponentType();
    }

    public static Object arrayGrow(Object array, int increment)
    {
        Class componentType = componentType(array);
        Object newArray = null;
        if (array != null)
        {
            int length = Array.getLength(array);
            int newLength = length + increment;
            newArray = Array.newInstance(componentType, newLength);
            System.arraycopy(array, 0, newArray, 0, length);
        }
        else
        {
            newArray = Array.newInstance(componentType, increment);
        }
        return newArray;
    }

    public static Object arrayAddUniqueElement(Object array, Object element)
    {
        for (int i = 0; i < ((Object[]) array).length; i++)
            if (((Object[]) array)[i] == element)
                return array;
        return arrayAddElement(array, element);
    }

    public static Object arrayAddElement(Object array, Object element)
    {
        if (array == null)
            return oneElementArray(element);
        if (element == null)
            return array;

        Class arrayComponentType = componentType(array);
        Class elementComponentType = element.getClass();
        Class componentType;
        if (elementComponentType.isAssignableFrom(arrayComponentType))
            componentType = elementComponentType;
        else if (arrayComponentType.isAssignableFrom(elementComponentType))
            componentType = arrayComponentType;
        else
        {
            StsException.systemError(StsMath.class, "arrayAddElement(array, element)",
                "arrayComponentType " + StsToolkit.getSimpleClassname(arrayComponentType) +
                    " and elementComponentType " + StsToolkit.getSimpleClassname(elementComponentType) + " are not compatible.");
            return array;
        }
        int length = Array.getLength(array);
        int newLength = length + 1;
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, 0, newArray, 0, length);
        Array.set(newArray, length, element);
        return newArray;
    }

    public static Object arrayAddElement(Object array, Object element, Class componentType)
    {
        if (array == null) return oneElementArray(element);
        int length = Array.getLength(array);
        int newLength = length + 1;
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, 0, newArray, 0, length);
        Array.set(newArray, length, element);
        return newArray;
    }

    public static Object arrayAddArray(Object oldArray, Object newArray)
    {
        if (oldArray == null)
        {
            return newArray;
        }
        if (newArray == null)
        {
            return oldArray;
        }
        int oldLength = Array.getLength(oldArray);
        int newLength = Array.getLength(newArray);
        Class componentType = componentType(oldArray);
        int length = oldLength + newLength;
        Object array = Array.newInstance(componentType, length);
        System.arraycopy(oldArray, 0, array, 0, oldLength);
        System.arraycopy(newArray, 0, array, oldLength, newLength);
        return array;
    }

    public static Object arrayAddArray(Object oldArray, Object newArray, int newLength)
    {
        if (oldArray == null)
        {
            return newArray;
        }
        if (newArray == null)
        {
            return oldArray;
        }
        int oldLength = Array.getLength(oldArray);
        Class componentType = componentType(oldArray);
        int length = oldLength + newLength;
        Object array = Array.newInstance(componentType, length);
        System.arraycopy(oldArray, 0, array, 0, oldLength);
        System.arraycopy(newArray, 0, array, oldLength, newLength);
        return array;
    }

    public static Object arrayAddArray(Object oldArray, Object newArray, Class componentType)
    {
        if (oldArray == null)
        {
            if (newArray != null && Array.getLength(newArray) > 0)
                return newArray;
            else
                return Array.newInstance(componentType, 0);
        }
        if (newArray == null)
        {
            if (oldArray != null && Array.getLength(oldArray) > 0)
                return oldArray;
            else
                return Array.newInstance(componentType, 0);
        }
        int oldLength = Array.getLength(oldArray);
        int newLength = Array.getLength(newArray);
        int length = oldLength + newLength;
        Object array = Array.newInstance(componentType, length);
        System.arraycopy(oldArray, 0, array, 0, oldLength);
        System.arraycopy(newArray, 0, array, oldLength, newLength);
        return array;
    }

    public static Object arrayAddSortedElement(Object array, Object element)
    {
        if (array == null)
        {
            return oneElementArray(element);
        }

        Comparable newObj = (Comparable) element;
        int length = Array.getLength(array);
        for (int n = 0; n < length; n++)
        {
            Object obj = Array.get(array, n);
            if (newObj.compareTo(obj) <= 0)
            {
                return arrayInsertElementBefore(array, element, n);
            }
        }
        return arrayInsertElementBefore(array, element, length);
    }

    /** insert this element before array[index] */
    public static Object arrayInsertElementBefore(Object array, Object element, int index)
    {
        if (index < 0)
        {
            StsException.systemError("StsMath.arrayInsertElememnt() failed. Index to insert (" + index + ") is invalid.");
            return array;
        }
        if (array == null)
        {
            return oneElementArray(element);
        }

        Class componentType = componentType(array);
        int length = Array.getLength(array);
        int newLength = length + 1;
        if (index > newLength)
        {
            StsException.systemError("StsMath.arrayInsertElememnt() failed. Index to insert (" + index + ") is invalid.");
            return null;
        }
        Object newArray = Array.newInstance(componentType, newLength);

        if (index > 0)
        {
            System.arraycopy(array, 0, newArray, 0, index);
        }
        if (length > index)
        {
            System.arraycopy(array, index, newArray, index + 1, length - index);
        }
        Array.set(newArray, index, element);
        return newArray;
    }

    public static Object arrayInsertElementAfter(Object array, Object element, Object prevElement)
    {
        if(array == null) return null;
        Object[] list = (Object[])array;
        int nObjects = list.length;
        for(int n = 0; n < nObjects; n++)
        {
            if(list[n] == prevElement)
                return arrayInsertElementBefore(array, element, n+1);
        }
        return null;
    }

    public static Object arrayInsertElementBefore(Object array, Object element, Object nextElement)
    {
        if(array == null) return null;
        Object[] list = (Object[])array;
        for(int n = 0; n < list.length; n++)
        {
            if(list[n] == nextElement)
                return arrayInsertElementBefore(array, element, n);
        }
        return null;
    }

    public static int[] intListInsertValue(int[] list, int value, int index)
    {
        if (index < 0)
        {
            StsException.systemError(StsMath.class, "intListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        if (list == null)
        {
            return new int[]{value};
        }

        int length = list.length;
        int newLength = length + 1;
        if (index >= newLength)
        {
            StsException.systemError(StsMath.class, "intListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        int[] newList = new int[newLength];

        if (index > 0)
        {
            System.arraycopy(list, 0, newList, 0, index);
        }
        if (length > index)
        {
            System.arraycopy(list, index, newList, index + 1, length - index);
        }
        newList[index] = value;
        return newList;
    }

    public static int[] intListAddSortedValue(int[] list, int value)
    {
        if (list == null)
        {
            return new int[]{value};
        }

        int length = list.length;
        for (int n = 0; n < length; n++)
        {
            if (value <= list[n])
                return intListInsertValue(list, value, n);
        }
        return intListInsertValue(list, value, length);
    }

    public static float[] floatListInsertValue(float[] list, float value, int index)
    {
        if (index < 0)
        {
            StsException.systemError(StsMath.class, "floatListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        if (list == null)
        {
            return new float[]{value};
        }

        int length = list.length;
        int newLength = length + 1;
        if (index >= newLength)
        {
            StsException.systemError(StsMath.class, "floatListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        float[] newList = new float[newLength];

        if (index > 0)
        {
            System.arraycopy(list, 0, newList, 0, index);
        }
        if (length > index)
        {
            System.arraycopy(list, index, newList, index + 1, length - index);
        }
        newList[index] = value;
        return newList;
    }

    public static long[] longListAddSortedValue(long[] list, long value)
    {
        if (list == null)
        {
            return new long[]{value};
        }

        int length = list.length;
        for (int n = 0; n < length; n++)
        {
            if (value <= list[n])
                return longListInsertValue(list, value, n);
        }
        return longListInsertValue(list, value, length);
    }

    public static String[] stringListInsertValue(String[] list, String value, int index)
    {
        if (index < 0)
        {
            StsException.systemError(StsMath.class, "stringListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        if (list == null)
        {
            return new String[]{value};
        }

        int length = list.length;
        int newLength = length + 1;
        if (index >= newLength)
        {
            StsException.systemError(StsMath.class, "stringListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        String[] newList = new String[newLength];

        if (index > 0)
        {
            for (int i = 0; i < index; i++)
                newList[i] = list[i];
        }
        if (length > index)
        {
            for (int i = index + 1; i < newLength; i++)
                newList[i] = list[i - 1];
        }
        newList[index] = value;
        return newList;
    }

    public static long[] longListInsertValue(long[] list, long value, int index)
    {
        if (index < 0)
        {
            StsException.systemError(StsMath.class, "longListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        if (list == null)
        {
            return new long[]{value};
        }

        int length = list.length;
        int newLength = length + 1;
        if (index >= newLength)
        {
            StsException.systemError(StsMath.class, "longListInsertValue", "Index to insert (" + index + ") is invalid.");
            return list;
        }
        long[] newList = new long[newLength];

        if (index > 0)
        {
            for (int i = 0; i < index; i++)
                newList[i] = list[i];
        }
        if (length > index)
        {
            for (int i = index + 1; i < newLength; i++)
                newList[i] = list[i - 1];
        }
        newList[index] = value;
        return newList;
    }

    public static float[] floatListAddSortedValue(float[] list, float value)
    {
        if (list == null)
        {
            return new float[]{value};
        }

        int length = list.length;
        for (int n = 0; n < length; n++)
        {
            if (value <= list[n])
                return floatListInsertValue(list, value, n);
        }
        return floatListInsertValue(list, value, length);
    }

    // array is initially all nulls; array is assumed to be filled if last entry is not null
    public static Object arrayAddSortedElementNullArray(Object array, int incLength, Object element)
    {
        if (array == null)
        {
            return createArrayAddElement(element, incLength);
        }

        Comparable newObj = (Comparable) element;
        int length = Array.getLength(array);
        for (int n = 0; n < length; n++)
        {
            Object obj = Array.get(array, n);
            if (newObj.compareTo(obj) <= 0)
            {
                return arrayInsertElement(array, element, n, incLength);
            }
        }
        return arrayInsertElement(array, element, length, incLength);
    }

    public static Object createOneElementArray(Object element)
    {
        if (element == null)
        {
            return null;
        }
        Class componentType = element.getClass();
        Object newArray = Array.newInstance(componentType, 1);
        Array.set(newArray, 0, element);
        return newArray;
    }

    public static Object createArrayAddElement(Object element, int incLength)
    {
        if (element == null)
        {
            return null;
        }
        Class componentType = element.getClass();
        Object newArray = Array.newInstance(componentType, incLength);
        Array.set(newArray, 0, element);
        return newArray;
    }

    /** insert this element before array[index]; lengthen array if last element is not null */
    public static Object arrayInsertElement(Object array, Object element, int index, int incLength)
    {
        if (index < 0)
        {
            StsException.systemError("StsMath.arrayInsertElememnt() failed. Index to insert (" + index + ") is invalid.");
            return array;
        }
        if (array == null)
        {
            return createArrayAddElement(element, incLength);
        }

        Class componentType = componentType(array);
        int length = Array.getLength(array);

        if (index >= length + 1)
        {
            StsException.systemError("StsMath.arrayInsertElememnt() failed. Index to insert (" + index + ") is invalid.");
            return null;
        }
        if (Array.get(array, length - 1) != null)
        {
            int newLength = length + incLength;

            Object newArray = Array.newInstance(componentType, newLength);

            if (index > 0)
            {
                System.arraycopy(array, 0, newArray, 0, index);
            }
            if (length > index)
            {
                System.arraycopy(array, index, newArray, index + 1, length - index);
            }
            array = newArray;
        }
        Array.set(array, index, element);
        return array;
    }

    static final boolean debug = false;

    public static Object arrayAddSortedElement(Object array, Object element, Comparator comparator)
    {
        if (array == null)
        {
            return createOneElementArray(element);
        }
        int length = Array.getLength(array);
        if (length == 0)
        {
//			if (mainDebug)System.out.println("Inserted element in null array.");
            return createOneElementArray(element);
        }

        for (int n = 0; n < length; n++)
        {
            Object object = Array.get(array, n);
            int compare = comparator.compare(object, element);
            if (compare < 0)
            {
                if (debug) System.out.println("Inserted element at index " + n + ".");
                array = arrayInsertElementBefore(array, element, n);
                return array;
            }
        }
        if (debug) System.out.println("Inserted element at end of array.");
        return arrayAddElement(array, element);
    }

    public static Object arrayAddElementNoRepeat(Object array, Object element, Class componentType)
    {
        if (array == null)
        {
            return oneElementArray(element);
        }

        int length = Array.getLength(array);

        for (int n = 0; n < length; n++)
        {
            if (Array.get(array, n) == element)
            {
                return array;
            }
        }

        int newLength = length + 1;
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, 0, newArray, 0, length);
        Array.set(newArray, length, element);
        return newArray;

    }

    public static Object arrayAddElementNoRepeat(Object array, Object element)
    {
        if (array == null)
        {
            return oneElementArray(element);
        }

        Class componentType = componentType(array);
        int length = Array.getLength(array);

        for (int n = 0; n < length; n++)
        {
            if (Array.get(array, n) == element)
            {
                return array;
            }
        }

        int newLength = length + 1;
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, 0, newArray, 0, length);
        Array.set(newArray, length, element);
        return newArray;
    }

    public static Object oneElementArray(Object element)
    {
        if (element == null)
        {
            return null;
        }
        Class componentType = element.getClass();
        Object newArray = Array.newInstance(componentType, 1);
        Array.set(newArray, 0, element);
        return newArray;
    }

    /**
     * Add element to next empty position in array of same class type. If array is
     * full, increment size by incSize. Routine requires that you keep track of
     * nElements
     *
     * @param array     Object array coming in
     * @param element   Object object being added
     * @param nElements int Current number of elements in the array.
     * @param incSize   int size to increment array by if it needs to be enlarged
     * @return Object  array returned
     */
    public static Object arrayAddElement(Object array, Object element, int nElements, int incSize)
    {
        Object newArray;
        int len = -1;

        try
        {
            if (element == null) return array;
            Class componentType = checkComponentType(array, element);
            if (componentType == null) return array;

            if (array == null)
                len = 0;
            else
                len = Array.getLength(array);

            if (len <= nElements)
            {
                newArray = Array.newInstance(componentType, len + incSize);
                if (nElements > 0)
                {
                    System.arraycopy(array, 0, newArray, 0, nElements);
                }
            }
            else
            {
                newArray = array;

            }
            Object currentElement = Array.get(newArray, nElements);
            if (currentElement != null)
            {
                StsException.systemError("StsMath.arrayAddElement() failed." +
                    " Adding element at a non-null position: have not incremented nElements properly.");
                return array;
            }
            Array.set(newArray, nElements, element);
            return newArray;
        }
        catch (Exception e)
        {
            StsException.outputException("StsMath.arrayAddElement() failed. array.len " + len + " nElements " +
                nElements, e, StsException.WARNING);
            return null;
        }
    }

    private static Class checkComponentType(Object array, Object element)
    {
        Class elementComponentType = element.getClass();
        if (array == null) return elementComponentType;
        Class arrayComponentType = array.getClass().getComponentType();
        if (!arrayComponentType.isAssignableFrom(elementComponentType))
        {
            StsException.systemError("StsMath.arrayAddElement() failed. New element class " + elementComponentType.getName() +
                " is not class or subClass of array class: " + arrayComponentType.getName());
            return null;
        }
        return arrayComponentType;
    }

    public static Object arrayDeleteLastElement(Object array)
    {
        Class componentType = componentType(array);
        int length = Array.getLength(array);
        if (length < 1)
        {
            return array;
        }
        int newLength = length - 1;
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, 0, newArray, 0, newLength);
        return newArray;
    }

    /** removes elements fromIndex thru toIndex-1 by moving nElements starting at toIndex down to start at fromIndex */
    public static Object arrayRemoveRange(Object array, int fromIndex, int toIndex, int nElements)
    {
        Class componentType = componentType(array);
        int length = Array.getLength(array);
        if (length < 1)
        {
            return array;
        }
        try
        {
            System.arraycopy(array, toIndex, array, fromIndex, nElements);
        }
        catch (Exception e)
        {
            System.out.println("StsMath.arrayRemoveRange failed. arrayLength: " +
                length + " fromIndex " + fromIndex + " toIndex " + toIndex + " nElements " + nElements);
        }
        for (int n = 0, i = toIndex; n < nElements; n++, i++)
        {
            Array.set(array, i, null);
        }
        return array;
    }

    public static Object arrayDeleteElement(Object array, Object deleteElement)
    {
        if (array == null || Array.getLength(array) == 0)
        {
            return array;
        }
        Class componentType = componentType(array);
        int length = Array.getLength(array);
        int index = arrayGetIndex(array, deleteElement);
        if (index == -1) return array;
        Object newArray = Array.newInstance(componentType, length - 1);
        System.arraycopy(array, 0, newArray, 0, index);
        System.arraycopy(array, index + 1, newArray, index, length - index - 1);
        return newArray;
    }

    public static int arrayGetIndex(Object array, Object element)
    {
        if (array == null || Array.getLength(array) == 0)
            return StsParameters.NO_MATCH;
        int length = Array.getLength(array);
        for (int n = 0; n < length; n++)
            if (Array.get(array, n) == element) return n;
        return StsParameters.NO_MATCH;
    }

    public static Object arrayDeleteElement(Object array, int index)
    {
        int length = Array.getLength(array);
        if (index < 0 || index >= length) return array;
        Class componentType = componentType(array);
        Object newArray = Array.newInstance(componentType, length - 1);
        if (index > 0) System.arraycopy(array, 0, newArray, 0, index);
        if (index < length - 1) System.arraycopy(array, index + 1, newArray, index, length - index - 1);
        return newArray;
    }

    public static Object arrayDeleteElementRange(Object array, int min, int max)
    {
        int length = Array.getLength(array);
        min = Math.max(0, min);
        max = Math.min(max, length - 1);
        int nDeletedElements = max - min + 1;
        if (nDeletedElements <= 0) return array;
        int nRemainingElements = length - nDeletedElements;
        Class componentType = componentType(array);
        Object newArray = Array.newInstance(componentType, nRemainingElements);
        if (min > 0) System.arraycopy(array, 0, newArray, 0, min);
        if (max < length - 1) System.arraycopy(array, max, newArray, min, length - max - 1);
        return newArray;
    }


    public static Object arrayDeleteElements(Object array, int[] indices)
    {
        int length = Array.getLength(array);
        Class componentType = componentType(array);
        int nDeleted = 0;
        for (int n = 0; n < indices.length; n++)
        {
            int index = indices[n];
            if (index >= 0 || index < length)
            {
                nDeleted++;
                Array.set(array, index, null);
            }
        }
        int nRemainingElements = length - nDeleted;
        Object newArray = Array.newInstance(componentType, nRemainingElements);
        int nElement = 0;
        for (int n = 0; n < length; n++)
        {
            Object element = Array.get(array, n);
            if (element != null) Array.set(newArray, nElement++, element);
        }
        return newArray;
    }

    public static int stringArrayContainsIndex(String[] array, String object)
    {
        if (array == null || object == null) return -1;
        int length = Array.getLength(array);
        for (int i = 0; i < length; i++)
            if (array[i].equalsIgnoreCase(object)) return i;
        return -1;
    }

    public static boolean arrayContains(Object array, Object object)
    {
        if (array == null || object == null) return false;
        int length = Array.getLength(array);
        for (int i = 0; i < length; i++)
            if (Array.get(array, i) == object) return true;
        return false;
    }

    public static boolean arrayContains(int[] array, int object)
    {
        if (array == null) return false;
        int length = array.length;
        for (int i = 0; i < length; i++)
            if (array[i] == object) return true;
        return false;
    }

    public static int index(Object[] objects, Object element)
    {
        for (int n = 0; n < objects.length; n++)
            if (objects[n] == element) return n;
        return StsParameters.NO_MATCH;
    }

    public static Object trimArray(Object array, int nElements)
    {
        Class componentType = componentType(array);
        if (componentType == null)
        {
            return null;
        }
        if (nElements <= 0)
        {
            nElements = 0;
        }
        Object newArray = Array.newInstance(componentType, nElements);
        System.arraycopy(array, 0, newArray, 0, nElements);
        return newArray;
    }

    public static Object trimArrayEnds(Object array, int iFirst, int iLast)
    {
        int newLength = iLast - iFirst + 1;
        if (newLength <= 0)
        {
            return null;
        }
        Class componentType = componentType(array);
        Object newArray = Array.newInstance(componentType, newLength);
        System.arraycopy(array, iFirst, newArray, 0, newLength);
        return newArray;
    }

    public static Object arraycopy(Object array)
    {
        Class componentType = null;

        try
        {
            componentType = componentType(array);
            if (componentType == null)
            {
                return null;
            }
            int length = Array.getLength(array);
            Object newArray = Array.newInstance(componentType, length);
            System.arraycopy(array, 0, newArray, 0, length);
            return newArray;
        }
        catch (Exception e)
        {
            StsException.outputException("StsMath.arraycopy() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public static Object arraycopy(Object array, Class componentType)
    {
        try
        {
            if (componentType == null) return null;
            if (array == null) return Array.newInstance(componentType, 0);
            int length = Array.getLength(array);
            if (length == 0) return Array.newInstance(componentType, length);
            Object newArray = Array.newInstance(componentType, length);
            System.arraycopy(array, 0, newArray, 0, length);
            return newArray;
        }
        catch (Exception e)
        {
            StsException.outputException("StsMath.arraycopy() failed.", e, StsException.WARNING);
            return Array.newInstance(componentType, 0);
        }
    }

    public static Object arraycopy(Object array, int start, int end)
    {
        Class componentType = null;

        try
        {
            int newLength = end - start + 1;
            if (newLength < 0)
            {
                return null;
            }
            int length = Array.getLength(array);
            if (length == 0)
            {
                return null;
            }
            for (int n = start; n < end; n++)
            {
                Object element = Array.get(array, n);
                if (element != null)
                {
                    componentType = element.getClass();
                    break;
                }
            }
            if (componentType == null)
            {
                return null;
            }
            Object newArray = Array.newInstance(componentType, newLength);
            System.arraycopy(array, start, newArray, 0, newLength);
            return newArray;
        }
        catch (Exception e)
        {
            StsException.outputException("StsMath.arraycopy() failed.", e, StsException.WARNING);
            return Array.newInstance(componentType, 0);
        }
    }

    public static Object arraycopy(Object array, int newLength)
    {
        Class componentType = null;

        try
        {
            if (newLength < 0)
            {
                return null;
            }
            int length = Array.getLength(array);
            if (length == 0)
            {
                return null;
            }

            for (int n = 0; n < length; n++)
            {
                Object element = Array.get(array, n);
                if (element != null)
                {
                    componentType = element.getClass();
                    break;
                }
            }
            if (componentType == null)
            {
                return null;
            }
            Object newArray = Array.newInstance(componentType, newLength);
            System.arraycopy(array, 0, newArray, 0, Math.min(length, newLength));
            return newArray;
            /*
                  componentType = componentType(array);
                  if(componentType == null) return null;

                  if(newLength <= 0)
                return Array.newInstance(componentType, 0);
                  int length = Array.getLength(array);
                  Object newArray = Array.newInstance(componentType, newLength);
                  System.arraycopy(array, 0, newArray, 0, Math.min(length, newLength));
                  return newArray;
                */
        }
        catch (Exception e)
        {
            StsException.outputException("StsMath.arraycopy() failed.", e, StsException.WARNING);
            return Array.newInstance(componentType, 0);
        }
    }

    public static Object arraycopy(Object array, int newLength, Class componentType)
    {
        try
        {
            if (componentType == null) return null;
            if (newLength <= 0) return Array.newInstance(componentType, 0);
            int length = Array.getLength(array);
            Object newArray = Array.newInstance(componentType, newLength);
            System.arraycopy(array, 0, newArray, 0, Math.min(length, newLength));
            return newArray;
        }
        catch (Exception e)
        {
            StsException.outputException("StsMath.arraycopy() failed.", e, StsException.WARNING);
            return Array.newInstance(componentType, 0);
        }
    }

    /**
     * toCastArray outputs an array of objects cast to the type of the first.
     * Requires that all objects be of the same type.
     * If this is array list of instances of class Orange, then the object returned
     * will be an Orange[] array.  Note that return type is Object, so caller must
     * cast to (Orange[]) in this case.
     *
     * @param arrayList ArrayList input of instance objects of the same class.
     * @return Object Output array of cast objects.  Caller must cast to an array
     *         of this type.
     */
    public static Object toCastArray(ArrayList arrayList)
    {
        int nObjects = arrayList.size();
        if (nObjects == 0) return null;
        Object firstObject = arrayList.get(0);
        Class componentType = firstObject.getClass();
        Object array = Array.newInstance(componentType, nObjects);
        arrayListCopy(arrayList, array, nObjects);
        return array;
    }

    private static void arrayListCopy(ArrayList arrayList, Object array, int nObjects)
    {
        Object[] listObjects = arrayList.toArray();
        System.arraycopy(listObjects, 0, array, 0, nObjects);
    }

    public static Iterator arrayIterator2D(Object[][] array)
    {
        return new ArrayIterator2D(array);
    }

    static public final class ArrayIterator2D implements Iterator
    {
        Object[][] array;
        int length1, length2;
        Object[] nextList;
        Object next;
        int index1, index2;

        ArrayIterator2D(Object[][] array)
        {
            this.array = array;
            length1 = array.length;
            index1 = 0;
            setNext();
        }

        private void setNext()
        {
            for (; index1 < length1; index1++)
            {
                nextList = array[index1];
                if (nextList == null)
                {
                    continue;
                }
                length2 = nextList.length;
                index2 = 0;
                next = nextList[0];
                return;
            }
            next = null;
        }

        public boolean hasNext()
        {
            return next != null;
        }

        public Object next()
        {
            Object current = next;
            index2++;
            if (index2 < length2)
            {
                next = array[index1][index2];
            }
            else
            {
                index1++;
                setNext();
            }
            return current;
        }

        public void remove()
        {
        }
    }

    public static Iterator arrayIterator(Object[] array)
    {
        return new ArrayIterator(array);
    }

    static public final class ArrayIterator implements Iterator
    {
        Object[] array;
        int length;
        Object next;
        int index;

        ArrayIterator(Object[] array)
        {
            this.array = array;
            length = array.length;
            index = 0;
            setNext();
        }

        private void setNext()
        {
            for (; index < length; index++)
            {
                next = array[index];
                return;
            }
            next = null;
        }

        public boolean hasNext()
        {
            return next != null;
        }

        public Object next()
        {
            Object current = next;
            index++;
            setNext();
            return current;
        }

        public void remove()
        {
        }
    }

    public static final int[] arrayMinMax(int[] array)
    {
        int min = StsParameters.largeInt;
        int max = -StsParameters.largeInt;
        min = array[0];
        max = array[0];
        for (int n = 0; n < array.length; n++)
        {
            if (array[n] < min)
                min = array[n];
            else if (array[n] > max)
                max = array[n];
        }
        return new int[]{min, max};
    }

    public static final float[] arrayMinMax(float[] array)
    {
        float min = StsParameters.largeFloat;
        float max = -StsParameters.largeFloat;
        min = array[0];
        max = array[0];
        for (int n = 0; n < array.length; n++)
        {
            if (array[n] < min)
                min = array[n];
            else if (array[n] > max)
                max = array[n];
        }
        return new float[]{min, max};
    }

    public static final double[] arrayMinMax(double[] array)
    {
        double min = StsParameters.largeDouble;
        double max = -StsParameters.largeDouble;
        min = array[0];
        max = array[0];
        for (int n = 0; n < array.length; n++)
        {
            if (array[n] < min)
                min = array[n];
            else if (array[n] > max)
                max = array[n];
        }
        return new double[]{min, max};
    }

    public static final long[] arrayMinMax(long[] array)
    {
        long min = 9999999999999l;
        long max = 0l;
        min = array[0];
        max = array[0];
        for (int n = 0; n < array.length; n++)
        {
            if (array[n] < min)
                min = array[n];
            else if (array[n] > max)
                max = array[n];
        }
        return new long[]{min, max};
    }

    //-----------------------------------
    // Vector operations
    //-----------------------------------

    public static final float[] interpolate(float[] a, float[] b, float f)
    {
        if (a == null || b == null) return null;
        if(f == 0.0f) return a;
        else if(f == 1.0f) return b;
        int length = Math.min(a.length, b.length);
        float[] c = new float[length];

        for (int n = 0; n < length; n++)
            c[n] = a[n] + f * (b[n] - a[n]);
        return c;
    }

    public static final float[] interpolate(StsPoint a, StsPoint b, float f)
    {
        if (a == null || b == null) return null;
        return interpolate(a.v, b.v, f);
    }

    public static final boolean interpolate(float[] a, float[] b, int length, float f, float[] result)
    {
        if (a == null || b == null) return false;
        for (int n = 0; n < length; n++)
            result[n] = a[n] + f * (b[n] - a[n]);
        return true;
    }

    public static final boolean interpolateBilinear(float[] v00, float[] v01, float[] v10, float[] v11, int length, float fx, float fy, float[] result)
    {
        if (v00 == null || v01 == null || v10 == null || v11 == null) return false;
        for (int n = 0; n < length; n++)
        {
            float v0 = v00[n] + fx * (v01[n] - v00[n]);
            float v1 = v10[n] + fx * (v11[n] - v10[n]);
            result[n] = v0 + fy * (v1 - v0);
        }
        return true;
    }

    public static final double[] interpolate(double[] a, double[] b, double f)
    {
        if (a == null || b == null) return null;
        int length = Math.min(a.length, b.length);
        double[] c = new double[length];

        for (int n = 0; n < length; n++)
            c[n] = a[n] + f * (b[n] - a[n]);
        return c;
    }

    public static final boolean interpolate(double[] a, double[] b, double f, int length, double[] result)
    {
        if (a == null || b == null) return false;
        for (int n = 0; n < length; n++)
            result[n] = a[n] + f * (b[n] - a[n]);
        return true;
    }

    public static final boolean interpolateUnsignedBytes254(byte[] a, byte[] b, float f, int length, byte[] result)
    {
        if (a == null || b == null) return false;
        for (int n = 0; n < length; n++)
        {
            int au = a[n] + 127;
            int bu = b[n] + 127;
            int r = Math.round(StsMath.minMax((au + f * (bu - au) - 127), 0, 254));
            result[n] = (byte) r;
        }
        return true;
    }

    public static final boolean interpolateVector(double[] origin, double[] vector, double f, int length, double[] result)
    {
        if (origin == null || vector == null) return false;
        for (int n = 0; n < length; n++)
            result[n] = origin[n] + f * vector[n];
        return true;
    }

    public static final boolean interpolateVectorsAddOrigin(double[] origin, double[] vector0, double[] vector1, double f, int length, double[] result)
    {
        if (origin == null || vector0 == null || vector1 == null) return false;
        for (int n = 0; n < length; n++)
            result[n] = vector0[n] + f * (vector1[n] - vector0[n]) + origin[n];
        return true;
    }

    public static final float[] add(float[] a, float[] b)
    {
        int length = Math.min(a.length, b.length);
        float[] c = new float[length];

        for (int n = 0; n < length; n++)
        {
            c[n] = a[n] + b[n];

        }
        return c;
    }

    public static final double[] add(double[] a, double[] b)
    {
        int length = Math.min(a.length, b.length);
        double[] c = new double[length];

        for (int n = 0; n < length; n++)
        {
            c[n] = a[n] + b[n];

        }
        return c;
    }

    public static final float[] subtract(float[] a, float[] b)
    {
        int length = Math.min(a.length, b.length);
        float[] c = new float[length];

        for (int n = 0; n < length; n++)
        {
            c[n] = a[n] - b[n];

        }
        return c;
    }

    public static final double[] subtract(double[] a, double[] b)
    {
        int length = Math.min(a.length, b.length);
        double[] c = new double[length];

        for (int n = 0; n < length; n++)
        {
            c[n] = a[n] - b[n];

        }
        return c;
    }

    public static final double[] subtractDivide(double[] a, double[] b, double f)
    {
        int length = Math.min(a.length, b.length);
        double[] c = new double[length];

        for (int n = 0; n < length; n++)
        {
            c[n] = (a[n] - b[n]) / f;

        }
        return c;
    }

    public static void sum(float[][] a, float[] sum)
    {
        int length = sum.length;
        int nPoints = a.length;
        for(int n = 0; n < nPoints; n++)
            for(int i = 0; i < length; i++)
                sum[i] += a[n][i];
    }

    public static void sum(double[][] a, float[] sum)
    {
        int length = sum.length;
        int nPoints = a.length;
        for(int n = 0; n < nPoints; n++)
            for(int i = 0; i < length; i++)
                sum[i] += (float)a[n][i];
    }

    public static void sum(float[][] a, float[] sum, float nullValue)
    {
        int length = sum.length;
        int nPoints = a.length;
        for(int n = 0; n < nPoints; n++)
            for(int i = 0; i < length; i++)
                if(a[n][i] != nullValue) sum[i] += a[n][i];
    }

    public static final double[] vector2(double x0, double y0, double x1, double y1)
    {
        return new double[]{x1 - x0, y1 - y0};
    }

    /** add this vector to existing vector */
    public static final void vectorAdd(float[] a, float[] da)
    {
        if (da == null || a == null) return;
        int length = Math.min(a.length, da.length);
        for (int n = 0; n < length; n++)
            a[n] += da[n];
    }


    static public double[] multByConstantAddPointStatic(double[] v0, double f, double[] v1)
    {
        int len = Math.min(v0.length, v1.length);
        double[] p = new double[len];
        for (int i = 0; i < len; i++)
            p[i] = v0[i] * f + v1[i];
        return p;
    }

    public static final long[] addPrimativeArrays(long[] a, long[] b)
    {
        if (a == null && b == null) return null;
        if (b == null) return a;
        if (a == null) return b;
        int length = a.length + b.length - 1;
        long[] c = new long[length];
        for (int n = 0; n < a.length; n++)
            c[n] = a[n];
        int cnt = a.length - 1;
        for (int n = 0; n < b.length; n++)
            c[cnt++] = b[n];
        return c;
    }

    public static final float[] addPrimativeArrays(float[] a, float[] b)
    {
        if (a == null && b == null) return null;
        if (b == null) return a;
        if (a == null) return b;
        int length = a.length + b.length - 1;
        float[] c = new float[length];
        for (int n = 0; n < a.length; n++)
            c[n] = a[n];
        int cnt = a.length - 1;
        for (int n = 0; n < b.length; n++)
            c[cnt++] = b[n];
        return c;
    }

    public static final float[] addVectorsNormalize(float[] a, float[] b)
    {
        if (b == null || a == null) return new float[3];
        int length = Math.min(a.length, b.length);
        float[] c = new float[length];
        float sumSq = 0;
        for (int n = 0; n < length; n++)
        {
            c[n] = a[n] + b[n];
            sumSq += c[n] * c[n];
        }
        if (sumSq == 0.0f) return c;
        float vectorLength = (float) Math.sqrt(sumSq);
        for (int n = 0; n < length; n++)
            c[n] /= vectorLength;
        return c;
    }

    public static final float[] getAddVectors(float[] a, float[] b)
    {
        if (b == null || a == null) return new float[3];
        int length = Math.min(a.length, b.length);
        float[] c = new float[length];
        for (int n = 0; n < length; n++)
            c[n] = a[n] + b[n];
        return c;
    }
    /**
     * add the first length components of these four vectors and normalize.  If the number of nonNull vectors is less
     * than min, then return a null vector
     *
     * @param a      first vector
     * @param b      second vector
     * @param c      third vector
     * @param d      fourth vector
     * @param length number of elements to sum
     * @param min    minimum number of nonNull vectors allowed
     * @return normalized vector or null if the number of nonNull vectors is less than min
     */
    public static final float[] addVectorsNormalize(float[] a, float[] b, float[] c, float[] d, int length, int min)
    {
        float[][] vectors = new float[4][];
        vectors[0] = a;
        vectors[1] = b;
        vectors[2] = c;
        vectors[3] = d;
        return addVectorsNormalize(vectors, length, min);
    }

    public static final float[] addVectorsNormalize(float[][] vectors, int length, int min)
    {
        int nNonNull = 0;
        for (int n = 0; n < 4; n++)
            if (vectors[n] != null) nNonNull++;
        if (nNonNull < min) return null;

        float[] sum = new float[length];
        float sumSq = 0;
        for (int i = 0; i < length; i++)
        {
            for (int n = 0; n < 4; n++)
                if (vectors[n] != null) sum[i] += vectors[n][i];
            sumSq += sum[i] * sum[i];
        }
        if (sumSq == 0.0f) return sum;
        float vectorLength = (float) Math.sqrt(sumSq);
        for (int n = 0; n < length; n++)
            sum[n] /= vectorLength;
        return sum;
    }

    public static float average(float[] values, float nullValue, int min)
    {
        int nNonNull = 0;
        for (int n = 0; n < 4; n++)
            if (values[n] != nullValue) nNonNull++;
        if (nNonNull < min) return nullValue;
        return average(values, nullValue);
    }

    public static float average(float[] values, float nullValue)
    {
        int nValues = 0;
        float sum = 0.0f;
        for (int i = 0; i < values.length; i++)
        {
            if (values[i] == nullValue)
                continue;
            sum += values[i];
            nValues++;
        }
        return sum / nValues;
    }

    public static float[] average(float[][] values, int length)
    {
        float[] average = new float[length];
        sum(values, average);
        int nPoints = values.length;
        for(int n = 0; n < length; n++)
            average[n] /= nPoints;
        return average;
    }

    public static float[] average(double[][] values, int length)
    {
        float[] average = new float[length];
        sum(values, average);
        int nPoints = values.length;
        for(int n = 0; n < length; n++)
            average[n] /= nPoints;
        return average;
    }

    public static float[] average(float[][] values, int length, float nullValue)
    {
        float[] average = new float[length];
        sum(values, average, nullValue);
        int nPoints = values.length;
        for(int n = 0; n < length; n++)
            average[n] /= nPoints;
        return average;
    }

    public static float average(float a, float b, float c, float d, float nullValue)
    {
        int nValues = 0;
        float sum = 0.0f;
        if (a != nullValue)
        {
            sum += a;
            nValues++;
        }
        if (b != nullValue)
        {
            sum += b;
            nValues++;
        }
        if (c != nullValue)
        {
            sum += c;
            nValues++;
        }
        if (d != nullValue)
        {
            sum += d;
            nValues++;
        }
        if (nValues == 0)
            return nullValue;
        else
            return sum / nValues;
    }

    public static final float[] copy(float[] a)
    {
        int length = a.length;
        float[] b = new float[length];
        System.arraycopy(a, 0, b, 0, length);
        return b;
    }

    public static final double[] copy(double[] a)
    {
        int length = a.length;
        double[] b = new double[length];
        System.arraycopy(a, 0, b, 0, length);
        return b;
    }

    public static final float[][] copy(float[][] a)
    {
        int length0 = a.length;
        float[][] b = new float[length0][];
        for (int n = 0; n < length0; n++)
        {
            int length1 = a[n].length;
            b[n] = new float[length1];
            System.arraycopy(a[n], 0, b[n], 0, a[n].length);
        }
        return b;
    }

    public static final double[][] copy(double[][] a)
    {
        int length0 = a.length;
        double[][] b = new double[length0][];
        for (int n = 0; n < length0; n++)
        {
            int length1 = a[n].length;
            b[n] = new double[length1];
            System.arraycopy(a[n], 0, b[n], 0, a[n].length);
        }
        return b;
    }

    public static final double[] copyDouble(float[] a)
    {
        int length = a.length;
        double[] b = new double[length];
        for (int n = 0; n < length; n++)
        {
            b[n] = (double) a[n];
        }
        return b;
    }

    public static final double[] copyDouble(float[] a, int length)
    {
        double[] b = new double[length];
        for (int n = 0; n < length; n++)
        {
            b[n] = (double) a[n];
        }
        return b;
    }

    public static final void increment(double[] x, double[] dx)
    {
        int length = x.length;
        for (int n = 0; n < length; n++)
            x[n] += dx[n];
    }

    public static final void scale(float[] a, float f)
    {
        for (int n = 0; n < a.length; n++)
        {
            a[n] *= f;
        }
    }

    public static final void scaleNormalize(float[] a)
    {
        float f = 0.0f;
        for (int n = 0; n < a.length; n++)
            f = Math.max(f, Math.abs(a[n]));
        if (f == 0.0f) return;
        float rf = 1.0f / f;
        for (int n = 0; n < a.length; n++)
            a[n] *= rf;
    }

    public static final void scale(double[] a, double f)
    {
        for (int n = 0; n < a.length; n++)
        {
            a[n] *= f;
        }
    }

    public static final float[] cross(float[] a, float[] b)
    {
        int minLength = Math.min(a.length, b.length);

        if (minLength != 3)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.cross3 failed.",
                "Vector length wrong: " + minLength));
            return null;
        }

        float[] cross = new float[3];

        cross[0] = a[1] * b[2] - a[2] * b[1];
        cross[1] = a[0] * b[2] - a[2] * b[0];
        cross[2] = a[0] * b[1] - a[1] * b[0];

        return cross;
    }

    static public float[] crossProduct(float[] vA, float[] vB)
    {
        if (vA == null || vB == null) return null;
        float[] cross = new float[3];
        cross[0] = vA[1] * vB[2] - vA[2] * vB[1];
        cross[1] = vA[2] * vB[0] - vA[0] * vB[2];
        cross[2] = vA[0] * vB[1] - vA[1] * vB[0];

        return cross;
    }

    static public float[] leftCrossProduct(float[] vA, float[] vB)
    {
        if (vA == null || vB == null) return null;
        float[] cross = new float[3];
        cross[0] = -vA[1] * vB[2] + vA[2] * vB[1];
        cross[1] = -vA[2] * vB[0] + vA[0] * vB[2];
        cross[2] = -vA[0] * vB[1] + vA[1] * vB[0];

        return cross;
    }

	static public float[] leftCrossProduct(double[] vA, double[] vB)
	{
		if (vA == null || vB == null) return null;
		float[] cross = new float[3];
		cross[0] = (float)(-vA[1] * vB[2] + vA[2] * vB[1]);
		cross[1] = (float)(-vA[2] * vB[0] + vA[0] * vB[2]);
		cross[2] = (float)(-vA[0] * vB[1] + vA[1] * vB[0]);
		return cross;
	}

	static public double[] leftCrossProductDouble(double[] vA, double[] vB)
	{
		if (vA == null || vB == null) return null;
		double[] cross = new double[3];
		cross[0] = (-vA[1] * vB[2] + vA[2] * vB[1]);
		cross[1] = (-vA[2] * vB[0] + vA[0] * vB[2]);
		cross[2] = (-vA[0] * vB[1] + vA[1] * vB[0]);
		return cross;
	}

    public static final double[] cross(double[] a, double[] b)
    {
        int minLength = Math.min(a.length, b.length);

        if (minLength != 3)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.cross3 failed.",
                "Vector length wrong: " + minLength));
            return null;
        }

        double[] cross = new double[3];

        cross[0] = a[1] * b[2] - a[2] * b[1];
        cross[1] = a[2] * b[0] - a[0] * b[2];
        cross[2] = a[0] * b[1] - a[1] * b[0];

        return cross;
    }

    public static final float[] cross3(float[] origin, float[] a, float[] b)
    {
        int minLength = StsMath.min3(origin.length, a.length, b.length);

        if (minLength != 3)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.cross3 failed.",
                "Vector length wrong: " + minLength));
            return null;
        }

        float[] cross = new float[3];
        float[] da = new float[3];
        float[] db = new float[3];

        for (int i = 0; i < 3; i++)
        {
            da[i] = a[i] - origin[i];
            db[i] = b[i] - origin[i];
        }

        cross[0] = da[1] * db[2] - da[2] * db[1];
        cross[1] = da[2] * db[0] - da[0] * db[2];
        cross[2] = da[0] * db[1] - da[1] * db[0];

        return cross;
    }

    public static final double[] cross3(double[] origin, double[] a, double[] b)
    {
        int minLength = StsMath.min3(origin.length, a.length, b.length);

        if (minLength != 3)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.cross3 failed.",
                "Vector length wrong: " + minLength));
            return null;
        }

        double[] cross = new double[3];
        double[] da = new double[3];
        double[] db = new double[3];

        for (int i = 0; i < 3; i++)
        {
            da[i] = a[i] - origin[i];
            db[i] = b[i] - origin[i];
        }

        cross[0] = da[1] * db[2] - da[2] * db[1];
        cross[1] = da[2] * db[0] - da[0] * db[2];
        cross[2] = da[0] * db[1] - da[1] * db[0];

        return cross;
    }

    public static final float[] leftCross3(float[] origin, float[] a, float[] b)
    {
        if (origin == null || a == null || b == null)
        {
            return null;
        }

        int minLength = StsMath.min3(origin.length, a.length, b.length);

        if (minLength != 3)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.leftCross3 failed.",
                "Vector length wrong: " + minLength));
            return null;
        }

        float[] cross = new float[3];
        float[] da = new float[3];
        float[] db = new float[3];

        for (int i = 0; i < 3; i++)
        {
            da[i] = a[i] - origin[i];
            db[i] = b[i] - origin[i];
        }

        cross[0] = -da[1] * db[2] + da[2] * db[1];
        cross[1] = -da[2] * db[0] + da[0] * db[2];
        cross[2] = -da[0] * db[1] + da[1] * db[0];

        return cross;
    }

    public static final float cross2(float[] a, float[] b)
    {
        int minLength = Math.min(a.length, b.length);

        if (minLength != 2)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.cross3 failed.",
                "Vector length wrong: " + minLength));
            return StsParameters.nullValue;
        }

        return a[0] * b[1] - a[1] * b[0];
    }

    public static final double cross2(double[] a, double[] b)
    {
        int minLength = Math.min(a.length, b.length);

        if (minLength != 2)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.cross3 failed.",
                "Vector length wrong: " + minLength));
            return StsParameters.nullValue;
        }

        return a[0] * b[1] - a[1] * b[0];
    }

    static public double[] crossNormalize(double[] a, double[] b)
    {
        double[] cross = cross(a, b);
        normalize(cross);
        return cross;
    }

    static public float dot(float[] a, float[] b)
    {
        int minLength = Math.min(a.length, b.length);
        float dot = a[0] * b[0];
        for (int n = 1; n < minLength; n++)
        {
            dot += a[n] * b[n];

        }
        return dot;
    }

    static public double dot(double[] a, double[] b)
    {
        int minLength = Math.min(a.length, b.length);
        double dot = a[0] * b[0];
        for (int n = 1; n < minLength; n++)
        {
            dot += a[n] * b[n];

        }
        return dot;
    }

    static public double dot(double[] a, float[] b)
    {
        int minLength = Math.min(a.length, b.length);
        double dot = a[0] * b[0];
        for (int n = 1; n < minLength; n++)
        {
            dot += a[n] * b[n];

        }
        return dot;
    }


    static public double dotAbs(double[] a, double[] b)
    {
        int minLength = Math.min(a.length, b.length);
        double dot = Math.abs(a[0] * b[0]);
        for (int n = 1; n < minLength; n++)
        {
            dot += Math.abs(a[n] * b[n]);

        }
        return dot;
    }

    static public float[] horizontalNormal(float[] origin, float[] a, int direction)
    {
        if (origin == null || a == null)
        {
            return null;
        }

        int minLength = Math.min(origin.length, a.length);

        if (minLength < 2)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.horizontalNormal failed.", "Vector length wrong: " + minLength));
            return null;
        }

        float dx = a[0] - origin[0];
        float dy = a[1] - origin[1];
        float length = (float) Math.sqrt(dx * dx + dy * dy);
        dx /= length;
        dy /= length;

        if (direction >= 0)
        {
            return new float[]{dy, -dx, 0.0f};
        }
        else
        {
            return new float[]{-dy, dx, 0.0f};
        }
    }

    static public float[] horizontalNormal2D(float[] origin, float[] a)
    {
        if (origin == null || a == null)
        {
            return null;
        }

        int minLength = Math.min(origin.length, a.length);

        if (minLength < 2)
        {
            StsException.outputException(new StsException(StsException.WARNING, "StsMath.horizontalNormal failed.", "Vector length wrong: " + minLength));
            return null;
        }

        float dx = a[0] - origin[0];
        float dy = a[1] - origin[1];
        float length = (float) Math.sqrt(dx * dx + dy * dy);
        dx /= length;
        dy /= length;

        return new float[]{-dy, dx};
    }

    static public float normalizeVectorReturnLength(float[] v)
    {
        return normalizeVector(v);
    }

    static public float normalizeVector(float[] v)
    {
        int length = v.length;

        double sumSq = 0.0;
        for (int n = 0; n < length; n++)
        {
            sumSq += v[n] * v[n];

        }
        float vLength = (float) Math.sqrt(sumSq);

        if (vLength > 0.0f)
        {
            for (int n = 0; n < length; n++)
            {
                v[n] /= vLength;
            }
        }

        return vLength;
    }

    static public double[] normalizeSum(double[] v)
    {
        int length = v.length;
        double[] normalized = new double[length];
        double sum = 0.0;
        for (int n = 0; n < length; n++)
            sum += v[n];
        for (int n = 0; n < length; n++)
            normalized[n] = v[n] / sum;
        return normalized;
    }

    static public double normalize(double[] v)
    {
        int length = v.length;

        double sumSq = 0.0;
        for (int n = 0; n < length; n++)
        {
            sumSq += v[n] * v[n];

        }
        double vLength = Math.sqrt(sumSq);

        if (vLength > 0.0)
        {
            for (int n = 0; n < length; n++)
            {
                v[n] /= vLength;
            }
        }

        return vLength;
    }

    static public double normalize(float[] v)
    {
        int length = v.length;

        double sumSq = 0.0;
        for (int n = 0; n < length; n++)
        {
            sumSq += v[n] * v[n];

        }
        double vLength = Math.sqrt(sumSq);

        if (vLength > 0.0)
        {
            for (int n = 0; n < length; n++)
            {
                v[n] /= vLength;
            }
        }

        return (float)vLength;
    }

    static public boolean normalizeAmplitude(double[] v, double max)
    {
        if (max == 0.0) return false;
        double rmax = 1.0 / max;
        for (int n = 0; n < v.length; n++)
            v[n] *= rmax;
        return true;
    }

    static public boolean normalizeAmplitude(float[] v, double max)
    {
        if (v == null) return false;
        if (max == 0.0) return false;
        int nValues = v.length;
        if (nValues == 0) return false;
        double rmax = 1.0f / max;
        for (int n = 0; n < nValues; n++)
            v[n] *= rmax;
        return true;
    }

    static public boolean normalizeAmplitude(float[][] v, double max)
    {
        if (v == null) return false;
        int nVectors = v.length;
        for (int n = 0; n < nVectors; n++)
            if (!normalizeAmplitude(v[n], max)) return false;
        return true;
    }

    static public boolean normalizeAmplitude(float[] v)
    {
        if (v == null) return false;
        int nValues = v.length;
        if (nValues == 0) return false;
        float max = 0.0f;
        for (int n = 0; n < nValues; n++)
            max = Math.max(max, Math.abs(v[n]));
        if (max == 0.0f) return false;
        float rmax = 1.0f / max;
        for (int n = 0; n < nValues; n++)
            v[n] *= rmax;
        return true;
    }

    static public boolean normalizeAmplitude(double[] v)
    {
        if (v == null) return false;
        int nValues = v.length;
        if (nValues == 0) return false;
        double max = 0.0;
        for (int n = 0; n < nValues; n++)
            max = Math.max(max, Math.abs(v[n]));
        if (max == 0.0f) return false;
        double rmax = 1.0f / max;
        for (int n = 0; n < nValues; n++)
            v[n] *= rmax;
        return true;
    }

    static public boolean normalizeAmplitude(float[] v, int nValues)
    {
        if (v == null) return false;
        if (nValues == 0) return false;
        float max = 0.0f;
        for (int n = 0; n < nValues; n++)
            max = Math.max(max, Math.abs(v[n]));
        if (max == 0.0f) return false;
        float rmax = 1.0f / max;
        for (int n = 0; n < nValues; n++)
            v[n] *= rmax;
        return true;
    }

    static public long interpolateValue(long a, long b, float f)
    {
        if (a == b) return a;
        return a + (b - a) * (long)f;
    }

    static public float interpolateValue(float a, float b, float f)
    {
        if (a == b) return a;
        return a + (b - a) * f;
    }

    /**
     * given two vectors  a and b which define two variables as a function of a common index and a is monotonically increasing,
     * return an interpolated value of b based on the location of av in vector a.
     */
    static public float interpolateValue(float av, float[] a, float[] b)
    {
        if (a == null || b == null || a.length != b.length)
            return StsParameters.nullValue;
        int length;
        if (a.length > b.length)
        {
            StsException.systemError("StsMath.interpolatedValue error: length of vector A (" + a.length +
                ") is greater than length of vector B (" + b.length + "). Will use shorter length, but check!");
            length = b.length;
        }
        else
            length = a.length;
        float indexF = arrayIndexF(av, a);
        int i = (int) indexF;
        i = minMax(i, 0, length - 2);
        float f = indexF - i;
        return b[i] + f * (b[i + 1] - b[i]);
    }

    /**
     * given two vectors  a and b which define two variables as a function of a common index and a is monotonically increasing,
     * return an interpolated value of b based on the location of av in vector a.
     */
    static public float interpolateValue(float av, float[] a, float[] b, int length)
    {
        int last = length - 1;
        float av0, av1, bv0, bv1;

        if (av < a[0]) // extrapolate before first value
        {
            av1 = a[1] - a[0];
            av0 = a[0] - av;
            bv1 = b[1] - b[0];
            bv0 = bv1 * av0 / av1;
            return b[0] - bv0;
        }

        if (av > a[last]) // extrapolate after last value
        {
            av1 = av - a[last];
            av0 = a[last] - a[last - 1];
            bv0 = b[last] - b[last - 1];
            bv1 = bv0 * av1 / av0;
            return bv1 + b[last];
        }

        for (int i = 0; i <= last; i++)
        {
            if (av == a[i])
            {
                return b[i];
            }
            if (av > a[i])
            {
                continue;
            }
            if (av < a[i])
            {
                av1 = a[i] - a[i - 1];
                av0 = av - a[i - 1];
                bv1 = b[i] - b[i - 1];
                bv0 = bv1 * av0 / av1;
                return b[i - 1] + bv0;
            }
        }
        return StsParameters.nullValue; // shouldn't reach here
    }

    /**
     * given an array of points with parameter a at index ai, points[n].v[ai] and
     * parameter b at index bi, return an interpolated value of b based on the
     * location of av in values of a.
     */
    static public float interpolateValue(StsPoint[] points, float av, int ai, int bi)
    {
        if (points == null)
        {
            return StsParameters.nullValue;
        }
        int length = points.length;
        if (length == 0) return 0.0f;
        if (length == 1) return points[0].v[bi];
        float[] v0 = points[0].v;
        float[] v1 = points[1].v;
        int p = 1;
        while (av > v1[ai] && p < length - 1)
        {
            v0 = v1;
            v1 = points[++p].v;
        }
        return interpolateValue(av, v0, v1, ai, bi);
    }

    /**
     * given a location av in a vector of equally spaced points, where min is min of points, inc is increment,
     * and b are values at each point, return the value at location av.
     *
     * @param av
     * @param min
     * @param inc
     * @param b
     * @return
     */
    static public float interpolateValue(float av, float min, float inc, float[] b)
    {
        if (b == null || b.length < 2)
            return StsParameters.nullValue;
        int length = b.length;

        float indexF = (av - min) / inc;
        int i;
        if (indexF < 0.0f)
            i = 0;
        else if (indexF >= length - 1)
            i = length - 2;
        else
            i = (int) indexF;
        float f = indexF - i;
        return b[i] + f * (b[i + 1] - b[i]);
    }

    /*
    static public float interpolatedValue(StsPoint[] points, float av, int ai, int bi)
    {
     if (points == null)
     {
      return StsParameters.nullValue;
     }
     int length = points.length;
     if(length == 0) return 0.0f;
     if(length == 1) return points[0].v[bi];
     int last = length - 1;
     float av0, av1;
     float[] v0, v1;

     v0 = points[0].v;
     if (av < v0[ai]) // extrapolate before first value
     {
      v1 = points[1].v;
      return interpolateValue(av, v0, v1, ai, bi);
     }
     else if (av == v0[ai])
     {
      return v0[bi];
     }

     v1 = points[last].v;
     if (av > v1[ai]) // extrapolate after last value
     {
      v0 = points[last - 1].v;
      return interpolateValue(av, v0, v1, ai, bi);
     }
     else if (av == v1[ai])
     {
      return v1[bi];
     }

     v1 = points[0].v;
     for (int i = 1; i <= last; i++)
     {
      v0 = v1;
      v1 = points[i].v;

      if (av == v1[ai])
      {
    return v1[bi];
      }
      if (av > v1[ai])
      {
    continue;
      }
      else // av < v1[ai]
      {
    return interpolateValue(av, v0, v1, ai, bi);
      }
     }
     return StsParameters.nullValue;
    }
    */
    static public float[] interpolateValues(StsPoint[] points, float[] av, int ai, int bi)
    {
        if (points == null)
        {
            return new float[]
                {StsParameters.nullValue};
        }
        int length = points.length;
        if (length == 0) return new float[]
            {0.0f};
        if (length == 1) return new float[]
            {points[0].v[bi]};
        int last = length - 1;
        float av0, av1, bv0, bv1;

        float[] v0 = points[0].v;
        float[] v1 = points[1].v;
        int p = 1;
        int nValues = av.length;
        float[] values = new float[nValues];
        for (int n = 0; n < nValues; n++)
        {
            while (av[n] > v1[ai] && p < length - 1)
            {
                v0 = v1;
                v1 = points[++p].v;
            }
            values[n] = interpolateValue(av[n], v0, v1, ai, bi);
        }
        return values;
    }

    static public float[] interpolateValues(StsPoint[] points, int nValues, float zMin, float zInc, int ai, int bi)
    {
        if (points == null)
        {
            return new float[]
                {StsParameters.nullValue};
        }
        int length = points.length;
        if (length == 0) return new float[]
            {0.0f};
        if (length == 1) return new float[]
            {points[0].v[bi]};

        float[] v0 = points[0].v;
        float[] v1 = points[1].v;
        int p = 1;
        float[] values = new float[nValues];
        float z = zMin;
        for (int n = 0; n < nValues; n++, z += zInc)
        {
            while (z > v1[ai] && p < length - 1)
            {
                v0 = v1;
                v1 = points[++p].v;
            }
            values[n] = interpolateValue(z, v0, v1, ai, bi);
        }
        return values;
    }

    public static final float interpolateValue(float av, float[] v0, float[] v1, int ai, int bi)
    {
        float av1 = v1[ai] - v0[ai];
        if (av1 == 0.0f)
        {
            //System.out.println("Bad juju!");
            return v0[bi];
        }
        float av0 = av - v0[ai];
        float bv1 = v1[bi] - v0[bi];
        float bv0 = bv1 * av0 / av1;
        return v0[bi] + bv0;
    }

    /**
     * Given two points in n-dimensional space v0 and v1, and a value to be interpolated av which
     * is defined in dimension ai, return the interpolated value in dimension bi
     *
     * @param av interpolating value
     * @param v0 point 0
     * @param v1 point 1
     * @param ai interpolating dimension
     * @param bi result dimension
     * @return
     */
    public static final float interpolateValueNoExtrapolation(float av, float[] v0, float[] v1, int ai, int bi)
    {
        if (v0[ai] < v1[ai])
        {
            if (av < v0[ai] || av >= v1[ai])
                return StsParameters.nullValue;
        }
        else
        {
            if (av < v1[ai] || av >= v0[ai])
                return StsParameters.nullValue;
        }

        float av1 = v1[ai] - v0[ai];
        if (av1 == 0.0f)
        {
            //System.out.println("Bad juju!");
            return v0[bi];
        }
        float av0 = av - v0[ai];
        float bv1 = v1[bi] - v0[bi];
        float bv0 = bv1 * av0 / av1;
        return v0[bi] + bv0;
    }

    /**
     * given an array of points with parameter a at index ai, points[n].v[ai] and
     * parameter b at index bi, return an interpolated value of b based on the
     * location of av in values of a.
     */
    static public StsPoint interpolatePoint(StsPoint[] points, float av, int ai)
    {
        if (points == null) return null;
        int length = points.length;
        if (length == 0) return null;
        if (length == 1) return points[0];

        StsPoint v0 = points[0];
        StsPoint v1 = points[1];
        int p = 1;
        while (av > v1.v[ai] && p < length - 1)
        {
            v0 = v1;
            v1 = points[++p];
        }
        return interpolatePoint(av, v0, v1, ai);
    }

    static public StsPoint[] interpolatePoints(StsPoint[] points, int nValues, float zMin, float zInc, int ai, int bi)
    {
        if (points == null) return null;
        int length = points.length;
        if (length == 0) return null;
        StsPoint[] values = new StsPoint[nValues];
        if (length == 1)
        {
            float z = zMin;
            for (int n = 0; n < nValues; n++, z += zInc)
                values[n] = new StsPoint(points[0]);
        }
        else
        {
            StsPoint v0 = points[0];
            StsPoint v1 = points[1];
            int p = 1;
            float z = zMin;
            for (int n = 0; n < nValues; n++, z += zInc)
            {
                while (z > v1.v[ai] && p < length - 1)
                {
                    v0 = v1;
                    v1 = points[++p];
                }
                values[n] = interpolatePoint(z, v0, v1, ai);
            }
        }
        return values;
    }

    static public double[][] interpolatePoints(double[][] points0, double[][] points1, float f)
    {
        if (points0 == null || points1 == null) return null;
        int length = Math.min(points0.length, points1.length);
        double[][] points = new double[length][];
        for (int n = 0; n < length; n++)
            points[n] = StsMath.interpolate(points0[n], points1[n], f);
        return points;
    }

    /*
      if (points == null)
      {
       return null;
      }
      int length = points.length;
      int last = length - 1;
      float av0, av1, bv0, bv1;
      StsPoint p0, p1;

      p0 = points[0];
      if (av < p0.v[ai]) // extrapolate before first value
      {
       p1 = points[1];
       return interpolatePoint(av, p0, p1, ai);
      }
      else if (av == p0.v[ai])
      {
       return p0;
      }

      p1 = points[last];
      if (av > p1.v[ai]) // extrapolate after last value
      {
       p0 = points[last - 1];
       return interpolatePoint(av, p0, p1, ai);
      }
      else if (av == p1.v[ai])
      {
       return p1;
      }

      p1 = points[0];
      for (int i = 1; i <= last; i++)
      {
       p0 = p1;
       p1 = points[i];

       if (av == p1.v[ai])
       {
        return p1;
       }
       if (av > p1.v[ai])
       {
        continue;
       }
       else // av < v1[ai]
       {
        return interpolatePoint(av, p0, p1, ai);
       }
      }
      return null;
       }
      */
    private static final StsPoint interpolatePoint(float av, StsPoint p0, StsPoint p1, int ai)
    {
        float f = (av - p0.v[ai]) / (p1.v[ai] - p0.v[ai]);
        return StsPoint.staticInterpolatePoints(p0, p1, f);
    }

    //-----------------------------------
    // Sorting routines
    //-----------------------------------

    // quick sort of a string array (see K&R 2nd ed., p. 87)

    public static final boolean qsort(String[] strings)
    {
        return qsort(strings, null, null);
    }

    // quick sort of a string array with same sorting on an optional array
    public static final boolean qsort(String[] strings, Object[] array)
    {
        return qsort(strings, array, null);
    }

    // quick sort of a string array with same sorting on two optional arrays
    public static final boolean qsort(String[] strings, Object[] array1, Object[] array2)
    {
        if (strings == null || strings.length <= 1)
        {
            return false;
        }
        if (array1 != null)
        {
            if (array1.length != strings.length)
            {
                return false;
            }
        }
        if (array2 != null)
        {
            if (array1.length != strings.length)
            {
                return false;
            }
        }
        return qsort(strings, array1, array2, 0, strings.length - 1);
    }

    private static boolean qsort(String[] strings, Object[] array1, Object[] array2, int left, int right)
    {
        try
        {
            if (left >= right)
            {
                return true;
            }
            int lrd2 = (left + right) / 2;
            swap(strings, left, lrd2);
            swap(array1, left, lrd2);
            swap(array2, left, lrd2);
            int last = left;
            for (int i = left + 1; i <= right; i++)
            {
                if (strings[i].compareTo(strings[left]) < 0)
                {
                    last++;
                    swap(strings, last, i);
                    swap(array1, last, i);
                    swap(array2, last, i);
                }
            }
            swap(strings, left, last);
            swap(array1, left, last);
            swap(array2, left, last);
            if (!qsort(strings, array1, array2, left, last - 1))
            {
                return false;
            }
            if (!qsort(strings, array1, array2, last + 1, right))
            {
                return false;
            }
        }
        catch (ArrayIndexOutOfBoundsException e)
        {
            return false;
        }
        return true; // shouldn't get here
    }

    // quick sort of a string array (see K&R 2nd ed., p. 87)
    public static final boolean qsort(Comparable[] objects)
    {
        if (objects == null) return false;
        int length = objects.length;
        if (length <= 0) return false;
        if (length == 1) return true;
        return qsort(objects, 0, length - 1);
    }

    public static boolean qsort(Comparable[] objects, int left, int right)
    {
        try
        {
            if (left >= right)
            {
                return true;
            }
            int lrd2 = (left + right) / 2;
            swap(objects, left, lrd2);
            int last = left;
            for (int i = left + 1; i <= right; i++)
            {
                if (objects[i].compareTo(objects[left]) < 0)
                {
                    last++;
                    swap(objects, last, i);
                }
            }
            swap(objects, left, last);
            if (!qsort(objects, left, last - 1))
            {
                return false;
            }
            if (!qsort(objects, last + 1, right))
            {
                return false;
            }
        }
        catch (ArrayIndexOutOfBoundsException e)
        {
            return false;
        }
        return true; // shouldn't get here
    }

    public static final void swap(Object[] s, int i1, int i2) throws ArrayIndexOutOfBoundsException
    {
        if (s == null) return;
        if (i1 == i2) return;
        Object temp = s[i1];
        s[i1] = s[i2];
        s[i2] = temp;
    }

    //------------------------------
    //    copy operations
    //------------------------------

    static public float[][] copyFloatArray(float[][] oldArray)
    {
        int n0 = oldArray.length;
        float[][] newArray = new float[n0][];
        for (int n = 0; n < n0; n++)
        {
            int n1 = oldArray[n].length;
            newArray[n] = new float[n1];
            System.arraycopy(oldArray[n], 0, newArray[n], 0, n1);
        }
        return newArray;
    }

    static public double[][] copyDoubleArray(double[][] oldArray)
    {
        int n0 = oldArray.length;
        double[][] newArray = new double[n0][];
        for (int n = 0; n < n0; n++)
        {
            int n1 = oldArray[n].length;
            newArray[n] = new double[n1];
            System.arraycopy(oldArray[n], 0, newArray[n], 0, n1);
        }
        return newArray;
    }

    static public String[] copyStringArray(String[] oldArray)
    {
        if (oldArray == null)
        {
            return null;
        }
        int len = oldArray.length;
        String[] newArray = new String[len];
        for (int n = 0; n < len; n++)
        {
            if (oldArray[n] != null)
            {
                newArray[n] = new String(oldArray[n]);
            }
        }
        return newArray;
    }

    static public byte[][] copyByteArray(byte[][] oldArray)
    {
        int n0 = oldArray.length;
        byte[][] newArray = new byte[n0][];
        for (int n = 0; n < n0; n++)
        {
            int n1 = oldArray[n].length;
            newArray[n] = new byte[n1];
            System.arraycopy(oldArray[n], 0, newArray[n], 0, n1);
        }
        return newArray;
    }

    //------------------------------
    //    copy/conversion operations
    //------------------------------

    static public double[] convertFloatToDoubleArray(float[] oldArray)
    {
        int n0 = oldArray.length;
        double[] newArray = new double[n0];
        for (int n = 0; n < n0; n++)
        {
            newArray[n] = (double) oldArray[n];
        }
        return newArray;
    }

    static public float[] convertDoubleToFloatArray(double[] oldArray)
    {
        int n0 = oldArray.length;
        float[] newArray = new float[n0];
        for (int n = 0; n < n0; n++)
        {
            newArray[n] = (float) oldArray[n];
        }
        return newArray;
    }

    //------------------------------
    //    ordering operations
    //------------------------------

    static public void reverseOrder(Object[] list)
    {
        if (list == null)
        {
            return;
        }

        int n0 = 0;
        int n1 = list.length - 1;

        while (n0 < n1)
        {
            Object item = list[n0];
            list[n0] = list[n1];
            list[n1] = item;
            n0++;
            n1--;
        }
    }

    //------------------------------
    //    search operations
    //------------------------------

    // return index of arrayValue >= value
    // we assume array values are in increasing order and uniformly spaced
    // or approximately uniformly spaced

    static public int arrayMatchIndex(Object[] objects, Object match)
    {
        if (objects == null || match == null) return -1;

        for (int n = 0; n < objects.length; n++)
            if (objects[n].equals(match)) return n;
        return -1;
    }

    static public int arrayIndexAbove(float value, float[] array)
    {
        if (array == null)
        {
            return -1;
        }
        float minValue = array[0];
        if (value <= minValue)
        {
            return 0;
        }

        int length = array.length;
        float maxValue = array[length - 1];
        if (value >= maxValue)
        {
            return length - 1;
        }

        if (minValue == maxValue)
        {
            return -1;
        }

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        float indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > 0 && array[index - 1] >= value)
            {
                index--;

            }
        }
        else ///
        {
            while (index < length - 1 && array[index + 1] < value)
            {
                index++;
            }
        }
        return index + 1;
    }

    static public int arrayIndexAbove(long value, long[] array)
    {
        if (array == null)
        {
            return -1;
        }
        long minValue = array[0];
        if (value <= minValue)
        {
            return 0;
        }

        int length = array.length;
        long maxValue = array[length - 1];
        if (value >= maxValue)
        {
            return length - 1;
        }

        if (minValue == maxValue)
        {
            return -1;
        }

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        long indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > 0 && array[index - 1] >= value)
            {
                index--;
            }
        }
        else ///
        {
            while (index < length - 1 && array[index + 1] < value)
            {
                index++;
            }
        }
        return index + 1;
    }

    static public int arrayIndexAbove(double value, double[] array)
    {
        if (array == null)
        {
            return -1;
        }
        double minValue = array[0];
        if (value <= minValue)
        {
            return 0;
        }

        int length = array.length;
        double maxValue = array[length - 1];
        if (value >= maxValue)
        {
            return length - 1;
        }

        if (minValue == maxValue)
        {
            return -1;
        }

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        double indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > 0 && array[index - 1] >= value)
            {
                index--;
            }
        }
        else ///
        {
            while (index < length - 1 && array[index + 1] < value)
            {
                index++;
            }
        }
        return index + 1;
    }

    static public int arrayIndexBelow(float value, float[] array)
    {
        if (array == null)
        {
            return -1;
        }
        float minValue = array[0];
        if (value <= minValue)
        {
            return 0;
        }

        int length = array.length;
        float maxValue = array[length - 1];
        if (value >= maxValue)
        {
            return length - 1;
        }

        if (minValue == maxValue)
        {
            return -1;
        }

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        float indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > 0 && array[index] >= value)
            {
                index--;
            }
        }
        else ///
        {
            while (index < length - 1 && array[index + 1] < value)
            {
                index++;
            }
        }
        return index;
    }

    static public int arrayIndexBelow(long value, long[] array)
    {
        if (array == null)
        {
            return -1;
        }
        long minValue = array[0];
        if (value <= minValue)
        {
            return 0;
        }

        int length = array.length;
        long maxValue = array[length - 1];
        if (value >= maxValue)
        {
            return length - 1;
        }

        if (minValue == maxValue)
        {
            return -1;
        }

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        long indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > 0 && array[index] >= value)
            {
                index--;
            }
        }
        else ///
        {
            while (index < length - 1 && array[index + 1] < value)
            {
                index++;
            }
        }
        return index;
    }

    static public float arrayIndexF(float value, float[] array)
    {
        return arrayIndexF(value, array, 0, array.length - 1);
    }

    static public float arrayIndexF(float value, float[] array, int min, int max)
    {
        float minValue = array[min];
        if (value <= minValue) return min;
        int length = array.length;
        float maxValue = array[max];
        if (value >= maxValue) return max;

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        float indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > min && array[index] >= value)
            {
                index--;

            }
        }
        else ///
        {
            while (index < max && array[index + 1] < value)
            {
                index++;
            }
        }
        float indexInc = (value - array[index]) / (array[index + 1] - array[index]);
        return index + indexInc;
    }

    static public float arrayIndexF(double value, double[] array)
    {
        return arrayIndexF(value, array, 0, array.length - 1);
    }

    static public float arrayIndexF(double value, double[] array, int min, int max)
    {
        double minValue = array[min];
        if (value <= minValue) return (float) min;
        int length = array.length;
        double maxValue = array[length - 1];
        if (value >= maxValue) return (float) max;

        int index = (int) (max * (value - minValue) / (maxValue - minValue));
        double indexValue = array[index];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > min && array[index] >= value)
            {
                index--;
            }
        }
        else ///
        {
            while (index < max && array[index + 1] < value)
            {
                index++;
            }
        }
        float indexInc = (float) ((value - array[index]) / (array[index + 1] - array[index]));
        return index + indexInc;
    }

    /**
     * We assume that the float vector points[n].v[arrayIndex] is monotonically increasing and may be regularly spaced.
     * We wish to compute an interpolated point value in this array.
     *
     * @param arrayIndex  index of the point coordinate vector used in interpolation
     * @param value       value in the float vector points[n].v[arrayIndex] which defines the interpolation index and factor
     * @param points      array of points in which points[n].v[arrayIndex] is monotonically increasing (and may be equally spaced)
     * @param extrapolate allows extrapolation before the first point and after the last point
     * @return return an StsPoint which interpolated between the two bounding points
     */
    static public StsPoint interpolatePoint(float value, StsPoint[] points, int arrayIndex, boolean extrapolate)
    {
        StsPoint point;
        int i = 0;
        float f;

        if (points == null)
        {
            return null;
        }

        int nPoints = points.length;
        if (nPoints < 2)
        {
            return new StsPoint(points[0]);
        }
        if (points[0].v.length < arrayIndex + 1)
            return null;

        if (value < points[0].v[arrayIndex])
        {
            if (extrapolate)
            {
                float m0 = points[0].v[arrayIndex];
                float m1 = points[1].v[arrayIndex];
                i = 0;
                f = (value - m0) / (m1 - m0);
            }
            else
            {
                return null;
            }
        }
        else if (value > points[nPoints - 1].v[arrayIndex])
        {
            if (extrapolate)
            {
                float m0 = points[nPoints - 2].v[arrayIndex];
                float m1 = points[nPoints - 1].v[arrayIndex];
                i = nPoints - 2;
                f = (value - m0) / (m1 - m0);
            }
            else
            {
                return null;
            }
        }
        else
        {
            float indexF = arrayIndexF(value, points, arrayIndex);
            i = (int) indexF;
            i = StsMath.minMax(i, 0, nPoints - 2);
            f = indexF - i;
        }

        point = StsPoint.staticInterpolatePoints(points[i], points[i + 1], f);
        return point;
    }


    static public float arrayIndexF(float value, StsPoint[] points, int arrayIndex)
    {
        int index = arrayIndexBelow(value, points, arrayIndex);
        if (index < 0 || index > points.length - 1) return (float) index;
        float indexInc = (value - points[index].v[arrayIndex]) / (points[index + 1].v[arrayIndex] - points[index].v[arrayIndex]);
        return index + indexInc;
    }

    /**
     * returns index of first point spanning interval containing value cast as a float and the interpolation factor
     * between the two points
     *
     * @param value      value being interpolated
     * @param points     array of regularly-spaced points
     * @param arrayIndex index of coordinate in points being interpolated
     * @return return index of first point and interpolation factor for value between the two points
     */
/*
    static public float[] arrayIndexF(float value, StsPoint[] points, int arrayIndex)
    {
        int index = arrayIndexBelow(value, points, arrayIndex);
        if(index == -1) return new float[] { -1.0f, 0.0f };
        if(index >= points.length - 1)
            index = points.length - 2;
        float indexInc = (value - points[index].v[arrayIndex]) / (points[index + 1].v[arrayIndex] - points[index].v[arrayIndex]);
        return new float[] { (float)index, indexInc };
    }
*/
    static public int arrayIndexBelow(float value, StsPoint[] points, int arrayIndex)
    {
        if (points == null)
        {
            return -1;
        }
        float minValue = points[0].v[arrayIndex];
        if (value <= minValue)
        {
            return 0;
        }

        int length = points.length;
        float maxValue = points[length - 1].v[arrayIndex];
        if (value >= maxValue)
        {
            return length - 2;
        }

        if (minValue == maxValue)
        {
            return -1;
        }

        int index = (int) ((length - 1) * (value - minValue) / (maxValue - minValue));
        float indexValue = points[index].v[arrayIndex];

        if (indexValue == value)
        {
            return index;
        }
        else if (indexValue > value)
        {
            while (index > 0 && points[index].v[arrayIndex] >= value)
            {
                index--;

            }
        }
        else ///
        {
            while (index < length - 1 && points[index + 1].v[arrayIndex] < value)
            {
                index++;
            }
        }
        return index;
    }
    //------------------------------
    //    geometry operations
    //------------------------------

    // intersect two XY lines extended to infinity

    static public float[] linePVIntersectXYInfinite(float[] p0, float[] v0, float[] p1, float[] v1)
    {
        int length = p0.length;
        float denom, u;
        float[] ipoint = new float[length];

        denom = v0[1] * v1[0] - v0[0] * v1[1];
        if (denom == 0.0f)
        {
            return null;
        }

        u = ((p0[0] - p1[0]) * v1[1] - (p0[1] - p1[1]) * v1[0]) / denom;

        for (int i = 0; i < length; i++)
        {
            ipoint[i] = p0[i] + u * v0[i];

        }
        return ipoint;
    }

    // intersect two lines: line0 is from p0 to p0+v0 and line1 is from p1+v1
    // assume that x and y are coordinates 0 and 1; interpolate the other coordinates

    static public float[] linePVIntersectXYFinite(float[] p0, float[] v0, float[] p1, float[] v1, float[] factors)
    {
        float[] intersection = linePVIntersectXY(p0, v0, p1, v1, factors);
        if (intersection == null)
        {
            return null;
        }
        if (!betweenInclusive(factors[0], 0.0f, 1.0f))
        {
            return null;
        }
        if (!betweenInclusive(factors[1], 0.0f, 1.0f))
        {
            return null;
        }
        return intersection;
    }

    // same as lineIntersectXYInfinite but returns intersection factors
    static public float[] linePVIntersectXY(float[] p0, float[] v0, float[] p1, float[] v1, float[] intersectionFactors)
    {
        float denom = v0[1] * v1[0] - v0[0] * v1[1];
        if (denom == 0.0f)
        {
            return null;
        }

        int length = p0.length;
        float[] ipoint = new float[length];
        float u = ((p0[0] - p1[0]) * v1[1] - (p0[1] - p1[1]) * v1[0]) / denom;

        for (int i = 0; i < length; i++)
        {
            ipoint[i] = p0[i] + u * v0[i];

        }
        float w;
        if (Math.abs(v1[0]) > Math.abs(v1[1]))
        {
            w = (ipoint[0] - p1[0]) / v1[0];
        }
        else
        {
            w = (ipoint[1] - p1[1]) / v1[1];

        }
        intersectionFactors[0] = u;
        intersectionFactors[1] = w;

        return ipoint;
    }
    
    static public double[] lineIntersectXY(double[] p00, double[]p01, double[] p10, double[] p11, double[] intersectionFactors)
    {
        double[] v0 = new double[] { p01[0] - p00[0], p01[1] - p00[1] };
        double[] v1 = new double[] { p11[0] - p10[0], p11[1] - p10[1] };
        return linePVIntersectXY(p00, v0, p10, v1, intersectionFactors);
    }

	static public boolean lineIntersectXY(StsPoint point00, StsPoint point01, StsPoint point10, StsPoint point11)
	{
		double[] p00 = new double[] { point00.v[0], point00.v[1] };
		double[] p01 = new double[] { point01.v[0], point01.v[1] };
		double[] p10 = new double[] { point10.v[0], point10.v[1] };
		double[] p11 = new double[] { point11.v[0], point11.v[1] };
	    double[] v0 = new double[] { p01[0] - p00[0], p01[1] - p00[1] };
        double[] v1 = new double[] { p11[0] - p10[0], p11[1] - p10[1] };
        return linePVIntersectXY(p00, v0, p10, v1);
	}

	static public boolean linePVIntersectXY(double[] p0, double[] v0, double[] p1, double[] v1)
    {
		double[] intersectionFactors = new double[2];
		linePVIntersectXY(p0, v0, p1, v1, intersectionFactors);
		for(double factor : intersectionFactors)
			if(factor < 0.0 || factor > 1.0) return false;
		return true;
	}

    static public double[] linePVIntersectXY(double[] p0, double[] v0, double[] p1, double[] v1, double[] intersectionFactors)
    {
        double denom = v0[1] * v1[0] - v0[0] * v1[1];
        if (denom == 0.0f)
        {
            return null;
        }

        int length = p0.length;
        double[] ipoint = new double[length];
        double u = ((p0[0] - p1[0]) * v1[1] - (p0[1] - p1[1]) * v1[0]) / denom;

        for (int i = 0; i < length; i++)
        {
            ipoint[i] = p0[i] + u * v0[i];

        }
        double w;
        if (Math.abs(v1[0]) > Math.abs(v1[1]))
        {
            w = (ipoint[0] - p1[0]) / v1[0];
        }
        else
        {
            w = (ipoint[1] - p1[1]) / v1[1];

        }
		if(intersectionFactors != null)
		{
			intersectionFactors[0] = u;
			intersectionFactors[1] = w;
		}
        return ipoint;
    }
    /**
     * given a line from between two points, line0 and line1, return the location of the nearest point on this line to the given point.
     * length specifies the dimension of the line and points, 2 for xy and 3 for xyz. The location is return as the interpolating factor
     * between the two line points with 0.0at line0 and 1.0 at line1.  Thus if the nearest point is half way in between, this method
     * will return 0.5.
     *
     * @param line0
     * @param line1
     * @param point
     * @param dimension
     * @return
     */
    static public double nearestPointOnLineInterpolant(double[] line0, double[] line1, double[] point, int dimension)
    {
        double distSq1 = distanceSq(line0, point, dimension);
        double distSq2 = distanceSq(line1, point, dimension);
        double distSq12 = distanceSq(line0, line1, dimension);
        if (distSq12 < StsParameters.smallFloat)
        {
            return 0.5;
        }
        return 0.5 + 0.5 * (distSq1 - distSq2) / distSq12;
    }

    static public boolean isOnLine(float[] pointVector, float[] lineVector, float tol)
    {
        float lengthSq = lengthSq(lineVector);
        if (lengthSq < roundOff)
        {
            return false;
        }
        float dot = dot(pointVector, lineVector); // projection of vp into v
        float f = dot / lengthSq;
        if (f < roundOff || f > 1.0f + roundOff)
        {
            return false;
        }
        float cross = cross2(lineVector, pointVector);
        float ff = cross / lengthSq;
        return ff <= tol;
    }

    static public int getNearestIndexOfPointOnLine(StsPoint point, StsPoint[] linePoints, int dimension)
    {
        if (linePoints == null) return -1;
        int nLinePoints = linePoints.length;
        double minDistSq = StsParameters.largeDouble;
        int minIndex = -1;
        for (int i = 0; i < nLinePoints; i++)
        {
            float distSq = point.distanceSquared(linePoints[i], dimension);
            if (distSq < minDistSq)
            {
                minDistSq = distSq;
                minIndex = i;
            }
        }

        if (minIndex == nLinePoints - 1)
        {
            minIndex--;

        }
        double f = nearestPointOnLineInterpolant(linePoints[minIndex], linePoints[minIndex + 1], point, dimension);
        if (f < 0.0f)
        {
            if (minIndex == 0)
            {
                f = 0.0f;
            }
            else
            {
                minIndex--;
                f = nearestPointOnLineInterpolant(linePoints[minIndex], linePoints[minIndex + 1], point, dimension);
                if (f < 0.0f || f > 1.0f)
                {
                    minIndex++;
                    f = 0.0f;
                }
            }
        }
        else if (f > 1.0f)
        {
            if (minIndex == nLinePoints - 2)
            {
                f = 1.0f;
            }
            else
            {
                minIndex++;
                nearestPointOnLineInterpolant(linePoints[minIndex], linePoints[minIndex + 1], point, dimension);
                if (f < 0.0f || f > 1.0f)
                {
                    minIndex--;
                    f = 1.0f;
                }
            }
        }
        return minIndex;
    }

    static public StsPoint getNearestPointOnLine(StsPoint point, StsPoint[] linePoints, int dimension)
    {
        if (linePoints == null) return null;
        int nLinePoints = linePoints.length;
        double minDistSq = StsParameters.largeDouble;
        int minIndex = -1;
        for (int i = 0; i < nLinePoints; i++)
        {
            float distSq = point.distanceSquared(linePoints[i], dimension);
            if (distSq < minDistSq)
            {
                minDistSq = distSq;
                minIndex = i;
            }
        }

        if (minIndex == nLinePoints - 1)
        {
            minIndex--;

        }
        double f = nearestPointOnLineInterpolant(linePoints[minIndex], linePoints[minIndex + 1], point, dimension);
        if (f < 0.0f)
        {
            if (minIndex == 0)
            {
                f = 0.0f;
            }
            else
            {
                minIndex--;
                f = nearestPointOnLineInterpolant(linePoints[minIndex], linePoints[minIndex + 1], point, dimension);
                if (f < 0.0f || f > 1.0f)
                {
                    minIndex++;
                    f = 0.0f;
                }
            }
        }
        else if (f > 1.0f)
        {
            if (minIndex == nLinePoints - 2)
            {
                f = 1.0f;
            }
            else
            {
                minIndex++;
                nearestPointOnLineInterpolant(linePoints[minIndex], linePoints[minIndex + 1], point, dimension);
                if (f < 0.0f || f > 1.0f)
                {
                    minIndex--;
                    f = 1.0f;
                }
            }
        }
        return StsPoint.staticInterpolatePoints(linePoints[minIndex], linePoints[minIndex + 1], (float) f);
    }

    static public double nearestPointOnLineInterpolant(StsPoint linePoint0, StsPoint linePoint1, StsPoint point, int dimension)
    {
        double distSq1 = point.distanceSquared(linePoint0, dimension);
        double distSq2 = point.distanceSquared(linePoint1, dimension);
        double distSq12 = linePoint0.distanceSquared(linePoint1, dimension);
        if (distSq12 < StsParameters.smallFloat)
        {
            return 0.5;
        }
        return 0.5 + 0.5 * (distSq1 - distSq2) / distSq12;
    }

    /*
    static public boolean isOnLine(float[] pointVector, float[] lineVector, float angleMax)
    {

     float dot = dot(pointVector, lineVector); // projection of vp into v
     if(dot <= 0.0f) return false;
     float lineLength = length(lineVector);
     if(lineLength < roundOff) return false;
     if(dot >= lineLength) return false;
     float cross = cross2(lineVector, pointVector);
     if(cross == StsParameters.nullValue) return false;
     float angle = (float)( DEGperRAD*Math.atan2((double)cross, (double)dot) );
     return Math.abs(angle) < angleMax;
    }
    */
    //------------------------------
    //    Binomial coefficients
    //------------------------------

    /** Construct binomial coefficients using Pascal's pyramid scheme */
    public static final float[] binomCoefs(int noValues)
    {
        if (noValues < 1)
        {
            return null;
        }

        float[] b = new float[noValues];

        for (int n = 0; n < noValues; n++)
        {
            float bp = 1.0f;
            for (int m = 1; m < n; m++)
            {
                float bm = bp;
                bp = b[m];
                b[m] = bm + bp;
            }
            b[n] = 1.0f;
        }

        return b;
    }

    /** Divide binomial coefficients by sum */
    public static final float[] binomWeights(int noValues)
    {
        if (noValues < 1)
        {
            return null;
        }

        float[] b = binomCoefs(noValues);

        float sum = (float) Math.pow(2.0, (double) noValues - 1);
        float rSum = 1.0f / sum;

        for (int n = 0; n < noValues; n++)
        {
            b[n] = b[n] * rSum;

        }
        return b;
    }

    public static final String toString(double[] array)
    {
        int length = array.length;
        if (length == 3)
        {
            return new String(" x: " + array[0] + " y: " + array[1] + " z: " + array[2]);
        }
        else
        {
            return new String(" ");
        }
    }

    public static final String toString(float[] array)
    {
        int length = array.length;
        if (length == 3)
        {
            return new String(" x: " + array[0] + " y: " + array[1] + " z: " + array[2]);
        }
        else
        {
            return new String(" ");
        }
    }

    public static final String toString(String name, double[] array)
    {
        int length = array.length;
        switch(length)
        {
            case 0:
                return name;
            case 1:
                return name + " x: " + array[0];
            case 2:
                return name + " x: " + array[0] + " y: " + array[1];
            default:
                return name + " x: " + array[0] + " y: " + array[1] + " z: " + array[2];
        }
    }

    public static final String toString(String name, float[] array)
    {
        int length = array.length;
        switch(length)
        {
            case 0:
                return name;
            case 1:
                return name + " x: " + array[0];
            case 2:
                return name + " x: " + array[0] + " y: " + array[1];
            default:
                return name + " x: " + array[0] + " y: " + array[1] + " z: " + array[2];
        }
    }

    /**
     * between minValue and maxValue, create rounded down/up min/max at bold ticks.
     * There are approxNumberIntervals of minorTicks with majorTicks every 2, 5, or 10 minorTicks
     *
     * @return return floats: niceMin, niceMax, minorTickInterval, majorTickInterval
     */
    static public double[] niceScaleBold(double minValue, double maxValue, int approxNumberIntervals, boolean isLinear)
    {
        double niceMin, niceMax, niceInc, roundInc, boldInc;
        double logNiceMin, logNiceMax;
        if (approxNumberIntervals < 1)
            approxNumberIntervals = 1;
        if (minValue == maxValue)
        {
            minValue = minValue - (0.10 * minValue);
            maxValue = maxValue + (0.10 * maxValue);
        }
        if (minValue == 0.0 && maxValue == 0.0)
        {
            maxValue = 1.0;
        }
        if (isLinear)
        {
            double dif = Math.abs(maxValue - minValue);
            double baseNiceInc = dif / approxNumberIntervals;
            baseNiceInc = niceNumber(baseNiceInc, true);
            while (true)
            {
                niceInc = baseNiceInc;
                if (dif / niceInc <= approxNumberIntervals)
                {
                    break;
                }
                niceInc = 2 * baseNiceInc;
                if (dif / niceInc <= approxNumberIntervals)
                {
                    break;
                }
                niceInc = 5 * baseNiceInc;
                if (dif / niceInc <= approxNumberIntervals)
                {
                    break;
                }
                baseNiceInc *= 10;
            }
            boldInc = 10 * baseNiceInc;
            if (maxValue < minValue)
            {
                niceMin = intervalRoundUp(minValue, boldInc);
                niceMax = intervalRoundDown(maxValue, boldInc);
                niceInc = -niceInc;
                boldInc = -boldInc;
            }
            else
            {
                niceMin = intervalRoundDown(minValue, boldInc);
                niceMax = intervalRoundUp(maxValue, boldInc);
            }
        }
        else
        {
            double logMin = log10(minValue);
            double logMax = log10(maxValue);
            logNiceMin = intervalRoundDown(logMin, 1.0);
            logNiceMax = intervalRoundUp(logMax, 1.0);
            niceMin = Math.pow(10, logNiceMin);
            niceMax = Math.pow(10, logNiceMax);
            niceInc = 1.0;
            boldInc = 1.0;
        }
        return new double[]
            {niceMin, niceMax, niceInc, boldInc};
    }

    /**
     * between minValue and maxValue, create rounded down/up min/max values.
     * There are approxNumberIntervals of minorTicks with majorTicks every 2, 5, or 10 minorTicks
     *
     * @return return floats: niceMin, niceMax, minorTickInterval, majorTickInterval
     */
    static public double[] niceScale(double minValue, double maxValue, int approxNumberIntervals, boolean isLinear)
    {
        double niceMin, niceMax, niceInc, roundInc, boldInc;
        double logNiceMin, logNiceMax;

        if (!isValueOk(minValue) || !isValueOk(maxValue)) return new double[]{0.0, 1.0};

        if (minValue == maxValue) // we will get infinity and NAN in the return values.
        {
            minValue = minValue - (0.10 * minValue);
            maxValue = maxValue + (0.10 * maxValue);
        }
        if (minValue == 0.0 && maxValue == 0.0)
        {
            maxValue = 1.0;
        }
        // Compute scale factors
        if (isLinear)
        {
            double baseNiceInc = Math.abs((maxValue - minValue) / approxNumberIntervals);
            double logNiceInc = log10(baseNiceInc);
            logNiceInc = intervalRoundDown(logNiceInc, 1.0);
            baseNiceInc = Math.pow(10, logNiceInc);

            double dif = Math.abs(maxValue - minValue);
            while (true)
            {
                niceInc = baseNiceInc;
                if (dif / niceInc <= approxNumberIntervals)
                {
                    break;
                }
                niceInc = 2 * baseNiceInc;
                if (dif / niceInc <= approxNumberIntervals)
                {
                    break;
                }
                niceInc = 5 * baseNiceInc;
                if (dif / niceInc <= approxNumberIntervals)
                {
                    break;
                }
                baseNiceInc *= 10;
            }
            boldInc = 10 * baseNiceInc;

            if (maxValue < minValue)
            {
                niceMin = intervalRoundUp(minValue, niceInc);
                niceMax = intervalRoundDown(maxValue, niceInc);
                niceInc = -niceInc;
                boldInc = -boldInc;
            }
            else
            {
                niceMin = intervalRoundDown(minValue, niceInc);
                niceMax = intervalRoundUp(maxValue, niceInc);
            }
        }
        else
        {
            double logMin = log10(minValue);
            double logMax = log10(maxValue);
            logNiceMin = intervalRoundDown(logMin, 1.0);
            logNiceMax = intervalRoundUp(logMax, 1.0);
            niceMin = Math.pow(10, logNiceMin);
            niceMax = Math.pow(10, logNiceMax);
            niceInc = 1.0;
            boldInc = 1.0;
        }
        return new double[]{niceMin, niceMax, niceInc, boldInc};
    }

    static public boolean isValueOk(double value)
    {
        if (value >= StsParameters.largeDouble) return false;
        if (value <= -StsParameters.largeDouble) return false;
        return true;
    }

    /** find a "nice" number approximately equal to x. Rounded to nearest 1, 2, or 5. */
    static public double niceNumber(double x)
    {
        /* Exponent */
        int exp;
        /* Fractional part of x */
        double f;
        /* Nice, rounded fraction */
        double nf;

        double lx;
        double fx;
        if (x == 0.0)
        {
            return (0.0);
        }
        fx = Math.abs(x);
        lx = log10(fx);
        exp = floor(lx);
        f = fx / Math.pow(10.0, exp);

        /* Between 1 and 10 */

        if (x > 0.0)
        {
            if (f < 1.5)
            {
                nf = 1.0;
            }
            else if (f < 3.5)
            {
                nf = 2.0;
            }
            else if (f < 7.5)
            {
                nf = 5.0;
            }
            else
                nf = 10.0;
        }
        else
        {
            if (f > 7.5)
            {
                nf = 10.0;
            }
            else if (f > 3.5)
            {
                nf = 5.0;
            }
            else if(f > 1.5)
            {
                nf = 2.0;
            }
            else
                nf = 1.0;
        }
        if (x > 0.0)
        {
            return nf * Math.pow(10.0, exp);
        }
        else
        {
            return -nf * Math.pow(10.0, exp);
        }
    }

    /** find a "nice" number approximately equal to x. Rounded up or down. */
    static public double niceNumber(double x, boolean roundDown)
    {
        /* Exponent */
        int exp;
        /* Fractional part of x */
        double f;
        /* Nice, rounded fraction */
        double nf;

        double lx;
        double fx;
        if (x == 0.0)
        {
            return (0.0);
        }
        fx = Math.abs(x);
        lx = log10(fx);
        exp = floor(lx);
        f = fx / Math.pow(10.0, exp);

        /* Between 1 and 10 */

        if (roundDown && x > 0.0 || !roundDown && x < 0.0)
        {
            if (f < 2.0)
            {
                nf = 1.0;
            }
            else if (f < 5.0)
            {
                nf = 2.0;
            }
            else
            {
                nf = 5.0;
            }
        }
        else
        {
            if (f > 5.0)
            {
                nf = 10.0;
            }
            else if (f > 2.0)
            {
                nf = 5.0;
            }
            else
            {
                nf = 2.0;
            }
        }
        if (x > 0.0)
        {
            return nf * Math.pow(10.0, exp);
        }
        else
        {
            return -nf * Math.pow(10.0, exp);
        }
    }

    static public double log10(double x)
    {
        return Math.log(x) / Math.log(10);
    }

    public static final int numberLoopElements(int nPoints, int n0, int n1)
    {
        if (n0 <= n1)
        {
            return n1 - n0 + 1;
        }
        return n1 + nPoints - n0 + 1;
    }

    public static final int nextBaseTwoInt(int value)
    {
        int next = 1;
        while (value > next)
        {
            next *= 2;
        }
        return next;
    }

    /** Not Debugged!  Don't use! */
    static public float IeeeFloatBigEndian(byte[] bytes)
    {
        byte b1 = bytes[0];
        byte b2 = bytes[1];
        byte b3 = bytes[2];
        byte b4 = bytes[3];

        byte S = (byte) ((b1 & 0x80) >> 7);
        int E = ((b1 & 0x7f) << 1) + ((b2 & 0x80) >> 7);
        long F = ((b2 & 0x7f) << 16) + (b3 << 8) + b4;

        double A = 2.0;
        double B = 127.0;
        double C = 1.0;
        double D2 = -126.0;
        double e23 = 8388608.0; // 2^23

        double M = (double) F / e23;

        double F1;
        if (S == 0) F1 = 1.0;
        else F1 = -1.0;

        if (0 < E && E < 255) return (float) (F1 * (C + M) * Math.pow(A, E - B));
        else if (E == 0 && F != 0) return (float) (F1 * M * Math.pow(A, D2));
        else if (E == 255 && F != 0) return -1; // Not a number
        else if (E == 255 && F == 0 && S == 1) return Float.NEGATIVE_INFINITY; // -Infinity
        else if (E == 255 && F == 0 && S == 0) return Float.POSITIVE_INFINITY; // Infinity
        else if (E == 0 && F == 0 && S == 1) return -0;
        else return 0; // ( E == 0 && F == 0 && S == 0 )
    }

    /** Not Debugged!  Don't use! */
    static public float IeeeFloatLittleEndian(byte[] bytes)
    {
        byte b1 = bytes[3];
        byte b2 = bytes[2];
        byte b3 = bytes[1];
        byte b4 = bytes[0];

        byte S = (byte) ((b1 & 0x80) >> 7);
        int E = ((b1 & 0x7f) << 1) + ((b2 & 0x80) >> 7);
        long F = ((b2 & 0x7f) << 16) + (b3 << 8) + b4;

        double A = 2.0;
        double B = 127.0;
        double C = 1.0;
        double D2 = -126.0;
        double e23 = 8388608.0; // 2^23

        double M = (double) F / e23;

        double F1;
        if (S == 0) F1 = 1.0;
        else F1 = -1.0;

        if (0 < E && E < 255) return (float) (F1 * (C + M) * Math.pow(A, E - B));
        else if (E == 0 && F != 0) return (float) (F1 * M * Math.pow(A, D2));
        else if (E == 255 && F != 0) return -1; // Not a number
        else if (E == 255 && F == 0 && S == 1) return Float.NEGATIVE_INFINITY; // -Infinity
        else if (E == 255 && F == 0 && S == 0) return Float.POSITIVE_INFINITY; // Infinity
        else if (E == 0 && F == 0 && S == 1) return -0;
        else return 0; // ( E == 0 && F == 0 && S == 0 )
    }

    /** Not Debugged!  Don't use! */
    static public float IBMFloatBigEndian(byte[] bytes)
    {
        int b1 = bytes[0] & 0xff;
        int b2 = bytes[1] & 0xff;
        int b3 = bytes[2] & 0xff;
        int b4 = bytes[3] & 0xff;

        int S = (b1 & 0x80) >> 7;
        int E = (b1 & 0x7f);
        long F = (b2 << 16) + (b3 << 8) + b4;

        double A = 2.0;
        double B = 127.0;
        double C = 1.0;
        double D2 = -126.0;
        double e23 = 8388608.0; // 2^23

        double M = (double) F / e23;

        double F1;
        if (S == 0) F1 = 1.0;
        else F1 = -1.0;

        if (0 < E && E < 255) return (float) (F1 * (C + M) * Math.pow(A, E - B));
        else if (E == 0 && F != 0) return (float) (F1 * M * Math.pow(A, D2));
        else if (E == 255 && F != 0) return -1; // Not a number
        else if (E == 255 && F == 0 && S == 1) return Float.NEGATIVE_INFINITY; // -Infinity
        else if (E == 255 && F == 0 && S == 0) return Float.POSITIVE_INFINITY; // Infinity
        else if (E == 0 && F == 0 && S == 1) return -0;
        else return 0; // ( E == 0 && F == 0 && S == 0 )
    }

    /** Not Debugged!  Don't use! */
    static public float IBMFloatLittleEndian(byte[] bytes)
    {
        int b1 = bytes[3] & 0xff;
        int b2 = bytes[2] & 0xff;
        int b3 = bytes[1] & 0xff;
        int b4 = bytes[0] & 0xff;

        int S = (b1 & 0x80) >> 7;
        int E = (b1 & 0x7f);
        long F = (b2 << 16) + (b3 << 8) + b4;

        double A = 2.0;
        double B = 127.0;
        double C = 1.0;
        double D2 = -126.0;
        double e23 = 8388608.0; // 2^23

        double M = (double) F / e23;

        double F1;
        if (S == 0) F1 = 1.0;
        else F1 = -1.0;

        if (0 < E && E < 255) return (float) (F1 * (C + M) * Math.pow(A, E - B));
        else if (E == 0 && F != 0) return (float) (F1 * M * Math.pow(A, D2));
        else if (E == 255 && F != 0) return -1; // Not a number
        else if (E == 255 && F == 0 && S == 1) return Float.NEGATIVE_INFINITY; // -Infinity
        else if (E == 255 && F == 0 && S == 0) return Float.POSITIVE_INFINITY; // Infinity
        else if (E == 0 && F == 0 && S == 1) return -0;
        else return 0; // ( E == 0 && F == 0 && S == 0 )
    }

    /*
      static public int SingleByte2SEF(int FloatType, char[] bytes)
      {
       char S;
       int E;
       long F;
       char b1,b2,b3,b4;

       switch (FloatType)
       {
        case IEEE_SINGLE_FLOAT:
      if (MemoryByteOrder() == BigEndian)
      {
       b1 = bytes[0];
       b2 = bytes[1];
       b3 = bytes[2];
       b4 = bytes[3];
      }
      else
      {
       b1 = bytes[3];
       b2 = bytes[2];
       b3 = bytes[1];
       b4 = bytes[0];
      }

      S = (b1 & 0x80) >> 7;
      E = ( (b1 & 0x7f) << 1) + ( (b2 & 0x80) >> 7);
      F = ( (b2 & 0x7f) << 16) + (b3 << 8) + b4;

      break;
        case IBM_SINGLE_FLOAT:
      if (MemoryByteOrder() == BigEndian)
      {
       b1 = bytes[0];
       b2 = bytes[1];
       b3 = bytes[2];
       b4 = bytes[3];
      }
      else
      {
       b1 = bytes[3];
       b2 = bytes[2];
       b3 = bytes[1];
       b4 = bytes[0];
      }

      S = (b1 & 0x80) >> 7;
      E = (b1 & 0x7f);
      F = (b2 << 16) + (b3 << 8) + b4;

      break;
        case VAX_SINGLE_FLOAT:
      b1 = bytes[1];
      b2 = bytes[0];
      b3 = bytes[3];
      b4 = bytes[2];

      S = (b1 & 0x80) >> 7;
      E = ( (b1 & 0x7f) << 1) + ( (b2 & 0x80) >> 7);
      F = ( (b2 & 0x7f) << 16) + (b3 << 8) + b4;
      break;
       }
       return SingleSEF2Float(FloatType, S, E, F);
      }

      static public double DoubleByte2SEF(int FloatType, char[] bytes)
      {
       char S;
       int E;
       long L1;
       long L2;
       char b1,b2,b3,b4,b5,b6,b7,b8;

       switch (FloatType)
       {
        case IEEE_DOUBLE_FLOAT:
      if (MemoryByteOrder() == BigEndian)
      {
       b1 = bytes[0];
       b2 = bytes[1];
       b3 = bytes[2];
       b4 = bytes[3];
       b5 = bytes[4];
       b6 = bytes[5];
       b7 = bytes[6];
       b8 = bytes[7];
      }
      else
      {
       b1 = bytes[7];
       b2 = bytes[6];
       b3 = bytes[5];
       b4 = bytes[4];
       b5 = bytes[3];
       b6 = bytes[2];
       b7 = bytes[1];
       b8 = bytes[0];
      }

      S = (b1 & 0x80) >> 7;
      E = ( (b1 & 0x7f) << 4) + ( (b2 & 0xf0) >> 4);
      L1 = ( (b2 & 0x0f) << 16) + (b3 << 8) + b4;
      L2 = (b5 << 24) + (b6 << 16) + (b7 << 8) + b8;

      break;
        case IBM_DOUBLE_FLOAT:
      if (MemoryByteOrder() == BigEndian)
      {
       b1 = bytes[0];
       b2 = bytes[1];
       b3 = bytes[2];
       b4 = bytes[3];
       b5 = bytes[4];
       b6 = bytes[5];
       b7 = bytes[6];
       b8 = bytes[7];
      }
      else
      {
       b1 = bytes[7];
       b2 = bytes[6];
       b3 = bytes[5];
       b4 = bytes[4];
       b5 = bytes[3];
       b6 = bytes[2];
       b7 = bytes[1];
       b8 = bytes[0];
      }

      S = (b1 & 0x80) >> 7;
      E = (b1 & 0x7f);
      L1 = (b2 << 16) + (b3 << 8) + b4;
      L2 = (b5 << 24) + (b6 << 16) + (b7 << 8) + b8;
      break;
       }

      }

      static public float SingleSEF2Float(int FloatType, char S, int E, long F)
      {
       double M, F1, A, B, C, D, D2, e23, e24;
       float ret;

       switch(FloatType)
       {
        case IEEE_SINGLE_FLOAT:

      A = 2.0;
      B = 127.0;
      C = 1.0;
      D2 = -126.0;
      e23 = 8388608.0;		// 2^23

      M = (double) F / e23;

      if ( S == 0) F1 = 1.0;
      else F1 = -1.0;

      if ( 0 < E && E < 255 ) ret = F1 * ( C + M ) * pow ( A, E - B ) ;
      else if ( E == 0 && F != 0 ) ret = F1 * M * pow ( A, D2 );
      else if ( E == 0 && F == 0 && S == 1 ) ret = -0;
      else if ( E == 0 && F == 0 && S == 0 ) ret = 0;
      else if ( E == 255 && F != 0) return -1; // Not a number
      else if ( E == 255 && F == 0 && S == 1 ) return Float.NEGATIVE_INFINITY; // -Infinity
      else if ( E == 255 && F == 0 && S == 0 ) return Float.POSITIVE_INFINITY; // Infinity

      break;

        case IBM_SINGLE_FLOAT:

      A = 16.0;
      B = 64.0;
      D = 0.69314718055994529;	// log2
      e24 = 16777216.0;		// 2^24

      M = (double) F / e24;

      if ( S == 0) F1 = 1.0;
      else F1 = -1.0;

      if ( S == 0 && E == 0 && F == 0 ) ret = 0;
      else ret = F1 * M * pow ( A, E - B ) ;

      break;

        case VAX_SINGLE_FLOAT:

      A = 2.0;
      B = 128.0;
      C = 0.5;
      e24 = 16777216.0;		// 2^24

      M = (double) F / e24;

      if ( S == 0 ) F1 = 1.0;
      else F1 = -1.0;

      if ( 0 < E ) ret = F1 * ( C + M ) * pow ( A, E - B ) ;
       else if ( E == 0 && S == 0 ) ret = 0;
      else if ( E == 0 && S == 1 ) return -1; // reserved

      break;
       }

       return ret;
      }
      */

    // IBM float bytes to IEEE float

    final public static void convertIBMFloatBytes(byte[] from, float[] to, int nFloatValues, boolean littleEndian)
    {
        convertIBMFloatBytes(from, to, 0, nFloatValues, littleEndian);
    }

    final public static void convertIBMFloatBytes(byte[] from, float[] to, int offset, int nFloatValues, boolean littleEndian)
    {
        int ch1, ch2, ch3, ch4, intBits;
        int i = offset;

        for (int n = 0; n < nFloatValues; n++)
        {
            ch1 = from[i++] & 0xff;
            ch2 = from[i++] & 0xff;
            ch3 = from[i++] & 0xff;
            ch4 = from[i++] & 0xff;
            if (littleEndian)
            {
                intBits = ((ch1 << 0) + (ch2 << 8) + (ch3 << 16) + (ch4 << 24));
            }
            else
            {
                intBits = ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
            }
            to[n] = convertIBMFloatIntBits(intBits);
        }
    }

    final public static float convertIBMFloatBytes(byte[] from, int offset, boolean littleEndian)
    {
        int intBits;

        int ch1 = from[offset++] & 0xff;
        int ch2 = from[offset++] & 0xff;
        int ch3 = from[offset++] & 0xff;
        int ch4 = from[offset++] & 0xff;
        if (littleEndian)
        {
            intBits = ((ch1 << 0) + (ch2 << 8) + (ch3 << 16) + (ch4 << 24));
        }
        else
        {
            intBits = ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
        }
        return convertIBMFloatIntBits(intBits);
    }

    final public static short convertBytesToShort(byte[] from, int offset, boolean littleEndian)
    {
        if (!littleEndian)
        {
            return (short) (((from[offset] & 0xff) << 8) | (from[offset + 1] & 0xff));
        }
        else
        {
            return (short) ((from[offset] & 0xff) | ((from[offset + 1] & 0xff) << 8));
        }
    }

    final public static float[] convertShortBytesToFloats(byte[] from, int offset, int nValues, boolean littleEndian)
    {
        float[] to = new float[nValues];
        for (int n = 0; n < nValues; n++, offset += 2)
        {
            if (!littleEndian)
            {
                to[n] = (short) ((from[offset] & 0xff) << 8) | (from[offset + 1] & 0xff);
            }
            else
            {
                to[n] = (short) (from[offset] & 0xff) | ((from[offset + 1] & 0xff) << 8);
            }
        }
        return to;
    }

    final public static boolean convertShortBytesToFloats(byte[] from, float[] to, int offset, int nValues, boolean littleEndian)
    {
        for (int n = 0; n < nValues; n++, offset += 2)
        {
            if (!littleEndian)
            {
                to[n] = (short) ((from[offset] & 0xff) << 8) | (from[offset + 1] & 0xff);
            }
            else
            {
                to[n] = (short) (from[offset] & 0xff) | ((from[offset + 1] & 0xff) << 8);
            }
        }
        return true;
    }

    final public static float convertIBMFloatBytes(byte[] from, int offset)
    {
        int ch1 = from[offset++] & 0xff;
        int ch2 = from[offset++] & 0xff;
        int ch3 = from[offset++] & 0xff;
        int ch4 = from[offset++] & 0xff;
        int intBits = ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
        return convertIBMFloatIntBits(intBits);
    }

    final public static short convertShortBytes(byte[] from, int offset)
    {
        return (short) (((from[offset] & 0xff) << 8) | ((from[offset + 1] & 0xff)));
    }

    final public static int convertIntBytes(byte[] from, int offset)
    {
        return (int) (((from[offset] & 0xff) << 24) | ((from[offset + 1] & 0xff) << 16) |
            ((from[offset + 2] & 0xff) << 8) | ((from[offset + 3] & 0xff)));
    }

    final public static int convertIntBytes(byte[] from, int offset, boolean littleEndian)
    {
        int intBits;

        int ch1 = from[offset] & 0xff;
        int ch2 = from[offset + 1] & 0xff;
        int ch3 = from[offset + 2] & 0xff;
        int ch4 = from[offset + 3] & 0xff;

        if (littleEndian)
        {
            intBits = ((ch1 << 0) | (ch2 << 8) | (ch3 << 16) | (ch4 << 24));
        }
        else
        {
            intBits = ((ch1 << 24) | (ch2 << 16) | (ch3 << 8) | (ch4 << 0));
        }
        return intBits;
    }

    final public static boolean convertIntBytes(byte[] from, float[] to, int nValues, boolean littleEndian)
    {
        int intBits;

        int offset = 0;
        for (int n = 0; n < nValues; n++, offset += 4)
        {
            int ch1 = from[offset] & 0xff;
            int ch2 = from[offset + 1] & 0xff;
            int ch3 = from[offset + 2] & 0xff;
            int ch4 = from[offset + 3] & 0xff;

            if (littleEndian)
            {
                intBits = ((ch1 << 0) | (ch2 << 8) | (ch3 << 16) | (ch4 << 24));
            }
            else
            {
                intBits = ((ch1 << 24) | (ch2 << 16) | (ch3 << 8) | (ch4 << 0));
            }
            to[n] = intBits;
        }
        return true;
    }

    final public static boolean convertIEEEBytesToFloats(byte[] from, float[] to, int nValues, boolean littleEndian)
    {
        int intBits;

        int offset = 0;
        for (int n = 0; n < nValues; n++, offset += 4)
        {
            int ch1 = from[offset] & 0xff;
            int ch2 = from[offset + 1] & 0xff;
            int ch3 = from[offset + 2] & 0xff;
            int ch4 = from[offset + 3] & 0xff;

            if (littleEndian)
            {
                intBits = ((ch1 << 0) | (ch2 << 8) | (ch3 << 16) | (ch4 << 24));
            }
            else
            {
                intBits = ((ch1 << 24) | (ch2 << 16) | (ch3 << 8) | (ch4 << 0));
            }
            to[n] = Float.intBitsToFloat(intBits);
        }
        return true;
    }

    final public static float convertIBMFloatBytes(byte[] from, boolean littleEndian)
    {
        int intBits;

        int ch1 = from[0] & 0xff;
        int ch2 = from[1] & 0xff;
        int ch3 = from[2] & 0xff;
        int ch4 = from[3] & 0xff;
        if (littleEndian)
        {
            intBits = ((ch1 << 0) + (ch2 << 8) + (ch3 << 16) + (ch4 << 24));
        }
        else
        {
            intBits = ((ch1 << 24) + (ch2 << 16) + (ch3 << 8) + (ch4 << 0));
        }
        return convertIBMFloatIntBits(intBits);
    }

    final public static float convertIBMFloatBytes(byte[] from)
    {
        return convertIBMFloatBytes(from, false);
    }

    // conversion of integer bits to an IBM float
    final public static float convertIBMFloatIntBits(int fconv)
    {
        int fmant, t;
        if (fconv != 0)
        {
            fmant = 0x00ffffff & fconv;
            if (fmant == 0) return 0.0f; // all zero bytes is illegal IBM format, but return it as a 0.0f
            t = (int) ((0x7f000000 & fconv) >> 22) - 130;
            int count = 0;
            while ((fmant & 0x00800000) == 0 && count++ < 32)
            {
                --t;
                fmant <<= 1;
            }
            if (t > 254)
            {
                fconv = (0x80000000 & fconv) | 0x7f7fffff;
            }
            else if (t <= 0)
            {
                fconv = 0;
            }
            else
            {
                fconv = (0x80000000 & fconv) | (t << 23) | (0x007fffff & fmant);
            }
        }
        return Float.intBitsToFloat(fconv);
    } // end of ibm_to_float(int)

    /**
     * *******************************************************************
     * float_to_ibm - convert between 32 bit IBM and IEEE floating numbers
     * **********************************************************************
     * Input:
     * from       input vector
     * n          number of floats in vectors
     * endian     =0 for little endian machine, =1 for big endian machines
     * <p/>
     * Output:
     * to         output vector, can be same as input vector
     * <p/>
     * **********************************************************************
     * Notes:
     * Up to 3 bits lost on IEEE -> IBM
     * <p/>
     * IBM -> IEEE may overflow or underflow, taken care of by
     * substituting large number or zero
     * <p/>
     * Only integer shifting and masking are used.
     * **********************************************************************
     * Credits:     CWP: Brian Sumner
     * *********************************************************************
     */
    static final public void convertIeeeBitsFloatToIbmBits(int[] ieee, int[] ibm, int n, boolean littleEndian)
    {
        int fconv, fmant, t;
        int i;

        for (i = 0; i < n; ++i)
            ibm[i] = convertIeeeBitsToIbmBits(ieee[i], littleEndian);
    }

    static final public int convertIeeeBitsToIbmBits(int ieee, boolean littleEndian)
    {
        int fconv = ieee;
        if (fconv != 0)
        {
            int fmant = (0x007fffff & fconv) | 0x00800000;
            int t = (int) ((0x7f800000 & fconv) >> 23) - 126;
            while ((t & 0x3) != 0)
            {
                ++t;
                fmant >>= 1;
            }
            fconv = (0x80000000 & fconv) | (((t >> 2) + 64) << 24) | fmant;
        }
        if (littleEndian)
            fconv = (fconv << 24) | ((fconv >> 24) & 0xff) |
                ((fconv & 0xff00) << 8) | ((fconv & 0xff0000) >> 8);

        return fconv;
    }

    static final public int convertIeeeFloatToIbmBits(float ieee, boolean littleEndian)
    {
        int ieeeBits = Float.floatToIntBits(ieee);
        return convertIeeeBitsToIbmBits(ieeeBits, false);
    }

    static final public short shortValue(int i, byte[] hdr, boolean isLittleEndian)
    {
//        return (short) (((hdr[i] & 0xff) << 8) | ((hdr[i + 1] & 0xff)));
        return (short) StsMath.convertBytesToShort(hdr, i, isLittleEndian);
    }

    static final public int intValue(int i, byte[] hdr, boolean isLittleEndian)
    {
//        return (int) ( ((hdr[i] & 0xff) << 24) | ((hdr[i + 1] & 0xff) << 16)
//                     | ((hdr[i + 2] & 0xff) << 8) | ((hdr[i + 3] & 0xff)));
        return StsMath.convertIntBytes(hdr, i, isLittleEndian);
    }

    static final public boolean verifyIBMFloatType(byte[] data, boolean isLittleEndian)
    {
        return verifyIBMFloatType(data, 0, isLittleEndian);
    }

    static final public boolean verifyIBMFloatType(byte[] data, int offset, boolean isLittleEndian)
    {
        int start = 0;
        int num;
        String hexString;

        for (int i = offset; i < data.length; i = i + 4)
        {
            num = StsMath.intValue(i, data, isLittleEndian);
            if (num == 0) return true; // bytes are null; illegal IBM format, but call it 0.0f and return ok

            int fmant = 0x00ffffff & num;
            if (fmant == 0)
            {
                if (debug)
                {
                    float value = convertIBMFloatBytes(data, offset, isLittleEndian);
                    System.out.println("Verify failed, fmant == 0. num = " + num + " Bytes are " + data[i] + " " + data[i + 1] + " " + data[i + 2] + " " + data[i + 3] + " value is " + value + " fmant is " + fmant);
                    return true;
                }
                return false;
            }
            /*
                hexString = Integer.toHexString(num);
                if(hexString.length() == 8)
                {
                    if(hexString.charAt(2) == '0')
                    {
                        if(debug)
                        {
                            float value = convertIBMFloatBytes(data, offset, getIsLittleEndian);
                            System.out.println("Verify failed, char[2] == '0'.  Bytes are " + data[i] + " " + data[i+1] + " " + data[i+2] + " " + data[i+3] + " value is " + value + " fmant is " + fmant);
                            return true;
                        }
                        return false;
                    }
                }
            */
        }
        return true;
    }

    /**
     * converts an int in the range 0 to 254 to an unsigned byte 0 to 254.
     * Since java doesn't support an unsigned byte, this will be converted on the cast
     * to 0 to 127, -128 to -1.  Graphics ("C" code) will read this as an unsigned byte.
     * 255 is reserved as the null value.
     */

    static final public byte unsignedIntToUnsignedByte254(int i)
    {
        if (i >= 255) i = 254;
        if (i < 0) i = 0;
        return (byte) i;
    }

    // Version of above which uses ByteBuffers
    // May be done in-place (from==to)
    // Obviously the buffer position(s) must be at a 4 byte word boundary.
    //     private static void ibm_to_float(ByteBuffer from, ByteBuffer to) {
    //         assert from.position()%4 == 0 : from.position();
    //         assert from.limit()%4 == 0 : from.limit();
    //         if ( to != from ) {
    //             assert to.position()%4 == 0 : to.position();
    //             assert to.limit()%4 == 0 : to.limit();
    //         }
    //         int fconv, i;
    //         float fout;
    //         for (i = from.position(); i < from.limit(); i += 4) {
    //             fconv = from.getInt(i); // Byte swapping occurs here if necessary
    //             fout = ibm_to_float(fconv);
    //             to.putFloat(i, fout);
    //         }
    //     } // end of ibm_to_float(ByteBuffer...)

    /** converts float value to 4 bytes */
    final public static void floatToBytes(byte[] byteArray, int offset, float value, boolean littleEndian)
    {
        int i = Float.floatToIntBits(value);
        if (!littleEndian)
        {
            byteArray[offset] = (byte) (i >> 24);
            byteArray[offset + 1] = (byte) (i >> 16);
            byteArray[offset + 2] = (byte) (i >> 8);
            byteArray[offset + 3] = (byte) i;
        }
        else
        {
            byteArray[offset + 3] = (byte) (i >> 24);
            byteArray[offset + 2] = (byte) (i >> 16);
            byteArray[offset + 1] = (byte) (i >> 8);
            byteArray[offset] = (byte) i;
        }
        return;
    }

    /**
     * converts a "C" unsigned byte (0 to 255) to a signed byte (-128 to 127)
     * Since Java doesn't support unsigned bytes, we are converting from signed equivalents.
     * unsigned    JavaSignedEquivalent  signedResult
     * 0              0                -128
     * 127            127              -1
     * 128            -128             0
     * 254            -2               126
     * 255            -1               127
     */

    final public static byte unsignedByteToSignedByte(byte b)
    {
        return (byte) (b ^ 0x80);
    }

    /** converts a "C" unsigned byte (0 to 255) to an int between -127 and 128 */
    final public static int unsignedByteToSignedInt(byte b)
    {
        return (int) (byte) (b ^ 0x80);
    }

    /**
     * converts a Java byte which is actually a "C" unsigned byte to its unsigned int equivalent.
     * signedByteToUnsignedInt((byte)0) = 0
     * signedByteToUnsignedInt((byte)127) = 127
     * signedByteToUnsignedInt((byte)-128) = 128
     * signedByteToUnsignedInt((byte)-1) = 255
     */
    final public static int signedByteToUnsignedInt(byte b)
    {
        return (int) (b & 0xFF);
    }

    final public static int unsignedCByteToUnsignedInt(byte b)
    {
        return (int) (b & 0xFF);
    }

    /** converts an unsigned int back to its original signedByte value */
    final public static byte unsignedIntToUnsignedByte(int i)
    {
        if (i > 255)
        {
            i = 255;
        }
        if (i < 0)
        {
            i = 0;
        }
        return (byte) i;
    }

    final public static float signedByteToFloat(byte signedByte, float dataMin, float dataMax)
    {
        float fValue = StsMath.signedByteToUnsignedInt(signedByte);
        return dataMin + (fValue / 254) * (dataMax - dataMin);
    }

    final public static float signedByteToFloatWithScale(byte signedByte, float scale, float scaleOffset)
    {
        float fValue = StsMath.signedByteToUnsignedInt(signedByte);
        return scale * fValue + scaleOffset;
    }

    final public static boolean signedBytesToFloatsWithScale(byte[] signedBytes, float scale, float scaleOffset, float[] floats)
    {
        if (signedBytes == null || floats == null || signedBytes.length != floats.length)
            return false;
        for (int n = 0; n < signedBytes.length; n++)
        {
            float fValue = StsMath.signedByteToUnsignedInt(signedBytes[n]);
            floats[n] = scale * fValue + scaleOffset;
        }
        return true;
    }

    /** converts an unsigned int back to its original signedByte value */
    final public static int intBits(byte[] bytes)
    {
        return (int) (((bytes[0] & 0xff) << 24) | ((bytes[1] & 0xff) << 16) | ((bytes[2] & 0xff) << 8) |
            ((bytes[3] & 0xff)));
    }

    static final public byte floatToUnsignedByte254(float value, float dataMin, float dataMax)
    {
        float scale = floatToUnsignedByteScale(dataMin, dataMax);
        float scaleOffset = floatToUnsignedByteScaleOffset(scale, dataMin);
        int i = Math.round(scale * value + scaleOffset);
        if (i >= 255) i = 254;
        if (i < 0) i = 0;
        return (byte) i;
    }

    static final public byte[] floatsToUnsignedBytes254(float[] floats, float dataMin, float dataMax)
    {
        int nValues = floats.length;
        byte[] bytes = new byte[nValues];
        if (!floatsToUnsignedBytes254(floats, dataMin, dataMax, bytes)) return null;
        return bytes;
    }

    static final public boolean floatsToUnsignedBytes254(float[] floats, float dataMin, float dataMax, byte[] bytes)
    {
        float scale = floatToUnsignedByteScale(dataMin, dataMax);
        float scaleOffset = floatToUnsignedByteScaleOffset(scale, dataMin);
        return floatsToUnsignedBytes254(floats, bytes, scale, scaleOffset);
    }

    static final public boolean floatsToUnsignedBytes254(float[] floats, byte[] bytes, float scale, float scaleOffset)
    {
        if (floats == null) return false;
        int nValues = floats.length;
        for (int n = 0; n < nValues; n++)
        {
            if (floats[n] == StsParameters.nullValue)
                bytes[n] = StsParameters.nullByte;
            else
            {
                int i = Math.round(scale * floats[n] + scaleOffset);
                if (i >= 255) i = 254;
                if (i < 0) i = 0;
                bytes[n] = (byte) i;
            }
        }
        return true;
    }

    /** If conversion is used repeatedly, compute scale and offset with above methods */
    static final public byte floatToUnsignedByte254WithScale(float value, float scale, float scaleOffset)
    {
        int i = Math.round(scale * value + scaleOffset);
        if (i >= 255) i = 254;
        if (i < 0) i = 0;
        return (byte) i;
    }
    static final public int floatToInt254WithScale(float value, float scale, float scaleOffset)
    {
        int i = Math.round(scale * value + scaleOffset);
        if (i >= 255) return 254;
        if (i < 0) return 0;
        return i;
    }

    static public final float floatToUnsignedByteScale(float dataMin, float dataMax)
    {
        return 254 / (dataMax - dataMin);
    }

    /**
     * if the data is roughly sinusoidal, then the offset is 0 and we scale the larger of -dataMin and dataMax to 1.
     * Otherwise, our data center is the half way between dataMin and dataMax
     *
     * @param dataMin
     * @param dataMax
     * @return
     */
    static public final float[] floatToNormalizedFloatScaleAndOffset(float dataMin, float dataMax)
    {
        float scale, offset = 0.0f;

        boolean hasZeroCrossing = dataMin <= 0 && dataMax >= 0;
        if (hasZeroCrossing)
        {
            if (dataMax >= -dataMin)
                scale = 1 / dataMax;
            else
                scale = -1 / dataMin;
        }
        else
        {
            scale = 2 / (dataMax - dataMin);
            offset = -dataMin * scale - 1.0f;
        }
        return new float[]{scale, offset};
    }

    /**
     * Java doesn't have an unsigned byte, but we can compute one anyway which can be handled by "C" graphics routines
     * which need an unsigned byte for colorscale index.
     */
    static public final float floatToUnsignedByteScale(float[] range)
    {
		if(range[0] == range[1]) return 0.0f;
        return 254 / (range[1] - range[0]);
    }

    static public final float floatToUnsignedByteScaleOffset(float scale, float dataMin)
    {
        return -dataMin * scale;
    }

    static public final float unsignedIntToFloatScale(float dataMin, float dataMax)
    {
        return (dataMax - dataMin) / 254;
    }

    static public final float unsignedIntToFloat(int value, float scale, float scaleOffset)
    {
        return value * scale + scaleOffset;
    }


    static final public byte floatToUnsignedByte254WithScale(float value, float scale, float dataMin, float dataMax)
    {
        if (value == StsParameters.nullValue) return StsParameters.nullByte;
        if (value <= dataMin) return (byte) 0;
        if (value >= dataMax) return (byte) 254;
        return (byte) ((value - dataMin) * scale);
    }

    static final public boolean floatToUnsignedByte254WithScale(FloatBuffer floats, ByteBuffer bytes, int length, float scale, float scaleOffset)
    {
        try
        {
            for (int n = 0; n < length; n++)
            {
                int i = Math.round(scale * floats.get() + scaleOffset);
                if (i >= 255) i = 254;
                if (i < 0) i = 0;
                bytes.put((byte) i);
            }
            return true;
        }
        catch (Exception e)
        {
            StsException.outputWarningException(StsMath.class, "floatToUnsignedByte254WithScale", e);
            return false;
        }
    }

    /** If conversion is used repeatedly, compute scale and offset with above methods */
    static final public byte floatToSignedByte254WithScale(float value, float scale, float scaleOffset)
    {
        int i = Math.round(scale * value + scaleOffset) - 127;
        if (i >= 127) i = 127;
        if (i < -127) i = -127;
        return (byte) i;
    }

    static float[] simpson6 = new float[]
        {3f / 8f, 7f / 6f, 23f / 24f};

    public static float simpsonCoefs(int n, int nPoints)
    {
        if (nPoints > 5)
        {
            if (n < 3) return simpson6[n];
            else if (n >= nPoints - 3) return simpson6[nPoints - n - 1];
            else return 1.0f;
        }
        else if (nPoints == 5)
        {
            if (n == 0 || n == 4) return 14f / 45f;
            else if (n == 1 || n == 3) return 64f / 45f;
            else return 24f / 45f;
        }
        else if (nPoints == 4)
        {
            if (n == 0 || n == 3) return 3f / 8f;
            else return 9f / 8f;
        }
        else if (nPoints == 3)
        {
            if (n == 1) return 4f / 3f;
            else return 1f / 3f;
        }
        else
            return 1f / 2f;
    }

    static public double interpolateGridBilinear(double[][] values, float rowF, float colF)
    {
        int i = 0, j = 0;
        double dx, dy;
        boolean debug = false;

        try
        {
            i = (int) rowF;
            dy = rowF - i;
            j = (int) colF;
            dx = colF - j;

            if (debug)
            {
                System.out.println("\txInt = " + rowF +
                    ", dx = " + dx + ", yInt = " +
                    colF +
                    ", dy = " + dy);
            }

            double weight = 0.0f;
            double z;
            double zWeighted = 0.0f;
            double w;

            w = (1.0f - dy) * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = values[i][j];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i][j] = " + z);
                }
            }
            w = dy * (1.0f - dx);
            if (w > StsParameters.roundOff)
            {
                z = values[i + 1][j];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i+1][j] = " + z);
                }
            }
            w = (1.0f - dy) * dx;
            if (w > StsParameters.roundOff)
            {
                z = values[i][j + 1];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i][j+1] = " + z);
                }
            }
            w = dy * dx;
            if (w > StsParameters.roundOff)
            {
                z = values[i + 1][j + 1];
                weight += w;
                zWeighted += w * z;
                if (debug)
                {
                    System.out.println("\tz[i+1][j+1] = " + z);
                }
            }
            if (weight > 0.0f)
            {
                return zWeighted / weight;
            }
            else
            {
                return doubleNullValue;
            }
        }
        catch (Exception e)
        {
            return doubleNullValue;
        }
    }

    static public String formatNumber(double number, int nSigDigits, int nColumns)
    {
        if (number == 0) return new String("0.0");

        int nIntegerDigits = (int) StsMath.log10(Math.abs(number)) + 1;
        boolean exponentialPattern = nIntegerDigits > nSigDigits || nIntegerDigits < -nSigDigits;

        DecimalFormat format = new DecimalFormat();

        if (exponentialPattern)
        {
            format.applyPattern("0.E0");
            format.setMinimumIntegerDigits(1);
            format.setMaximumFractionDigits(nSigDigits - 1);
        }
        else if (nIntegerDigits > 0)
        {
            format.setMinimumFractionDigits(0);
            format.setMaximumFractionDigits(Math.max(0, nSigDigits - nIntegerDigits));
        }
        else
        {
            format.setMinimumFractionDigits(0);
            format.setMaximumFractionDigits(nSigDigits - nIntegerDigits);
        }
        String numberString = format.format(number);
        int nChars = numberString.length();
        if (nChars > nColumns && !exponentialPattern)
        {
//			System.out.println("string " + numberString + " larger than " + nColumns + " columns. Revising...");
            format.applyPattern("0.E0");
            format.setMinimumIntegerDigits(1);
            format.setMaximumFractionDigits(Math.max(0, nColumns - 3));
        }
//		System.out.println(format.toPattern() + " number: " + number + " string: " + numberString);
        return numberString;
    }

    static public void convertInstantToAvgVelocity(float[] velocities, double scaleMultiplier, float dt)
    {
        double zSum = 0.0;
        double tSum = 0.0;
        int nValues = velocities.length;
        double f = scaleMultiplier * 0.5 * dt;
        double z1 = f * velocities[0];
        velocities[0] *= (float) scaleMultiplier;
        for (int n = 1; n < nValues; n++)
        {
            double z0 = z1;
            z1 = f * velocities[n];
            zSum += z0 + z1;
            tSum += dt;
            velocities[n] = (float) (zSum / tSum);
        }
    }

    static public void convertAvgToInstantVelocity(float[] velocities, double scaleMultiplier, float dt)
    {
        int nValues = velocities.length;
        double t = -dt;
        double tp1 = 0.0;
        double z = -dt*velocities[0];
        double zp1 = 0.0;
        double f = scaleMultiplier * 0.5 / dt;
        for (int n = 0; n < nValues-1; n++)
        {
            t = tp1;
            tp1 += dt;
            double zm1 = z;
            z = zp1;
            zp1 = tp1*velocities[n+1];
            velocities[n] = (float)(f*(zp1 - zm1));
        }
        velocities[nValues-1] = velocities[nValues-2];
    }

    public static boolean windowAverages(double[] values, int windowWidth, int offset, int nValues)
    {
        if (values == null) return false;
        int i;
        double sum = 0.0;
        double value;
        double[] averages = new double[nValues];
        // compute initial window for first point
        for (i = 0; i < windowWidth + 1; i++)
        {
            value = values[offset + i];
            //System.out.println("i " + i + " value "+ value);
            sum += value;
        }
        int nwin = windowWidth + 1;
        averages[0] = sum / nwin;

        // ramping on
        for (i = 1; i <= windowWidth; i++)
        {
            value = values[offset + i + windowWidth];
            sum += value;
            ++nwin;
            averages[i] = sum / nwin;
        }

        // middle range -- full avg window
        for (i = windowWidth + 1; i <= nValues - 1 - windowWidth; i++)
        {
            value = values[offset + i + windowWidth];
            sum += value;
            value = values[offset + i - windowWidth - 1];
            sum -= value;
            averages[i] = sum / nwin;
        }

        // ramping off
        for (i = nValues - windowWidth; i <= nValues - 1; i++)
        {
            value = values[offset + i - windowWidth];
            sum -= value; // rounding could make sum negative!
            --nwin;
            averages[i] = sum / nwin;
        }
        System.arraycopy(averages, 0, values, offset, nValues);
        return true;
    }

    static class TestComparator implements Comparator
    {
        TestComparator()
        {
        }

        public int compare(Object o1, Object o2)
        {
            int n1 = (int) ((StsPoint) o1).v[0];
            int n2 = (int) ((StsPoint) o2).v[0];
            return n2 - n1;
        }

        public boolean equals(Object o1)
        {
            return false;
        }
    }

    /*
      public static final void main(String[] args)
      {
      try
       {
        double v0 = -1;
        double s0 = 0;
        double v1 = 1;
        double s1 = 0;
        double[] coef = computeCoefs(v0, v1, s0, s1);
        double root = findRoot(coef, 0, 1, 0.005, true);
        System.out.println("root found " + root + " should be 0.5.");
        v0 = -0.9;
        coef = computeCoefs(v0, v1, s0, s1);
        root = findRoot(coef, 0, 1, 0.005, true);
        double[] checkRoot = evalPoly(coef, 1, root);
        System.out.println("root found " + root + " checked solution (should be zero) " + checkRoot[0]);
       }
       catch (Exception e)
       {
        e.printStackTrace();
       }
      }
      static final double[] computeCoefs(double v0, double v1, double s0, double s1)
      {
       double a = (s1 + s0) - 2*(v1 - v0);
       double b = v1 - s0 - v0 - a;
       double c = s0;
       double d = v0;
       return new double[] { d, c, b, a };
      }
      */

    public static final void main(String[] args)
    {
        double f = 0.99;
        double fn = niceNumber(f);
        System.out.println(f + " nice is " + fn);
        f = -.99;
        fn = niceNumber(f);
        System.out.println(f + " nice is " + fn);
    }
/*
    public static final void main(String[] args)
    {
        double v0 = 8;
        double v1 = -4;
        double s0 = -12;
        double s1 = -12;
        double a = 2 * (v0 - v1) + s0 + s1;
        double b = 3 * (v1 - v0) - 2 * s0 - s1;
        double c = s0;
        double d = v0;
        double[] coefs = new double[]{d, c, b, a};
        double f = StsMath.findRoot(coefs, 0, 1, 0.001, false);
        if (f != StsParameters.doubleNullValue)
        {
            double[] values = StsMath.evalPoly(coefs, 1, f);
            if (Math.abs(values[0]) > 2 * 0.001)
                StsException.systemError(StsMath.class, "main", "zero crossing calc failed. Value should be 0, but is " + values[0] +
                    " amp value: " + values[1]);
        }
    }
*/
    /*
public static final void main(String[] args)
{
int[] array = new int[]{0, 1, 2, 3, 4, 5, 6, 7};
int[] remaining = (int[]) StsMath.arrayDeleteElementRange(array, 0, 2);
remaining = (int[]) StsMath.arrayDeleteElementRange(array, 1, 5);
remaining = (int[]) StsMath.arrayDeleteElementRange(array, 0, 7);
remaining = (int[]) StsMath.arrayDeleteElementRange(array, 0, 8);
remaining = (int[]) StsMath.arrayDeleteElementRange(array, -3, -1);

remaining = (int[]) StsMath.arrayDeleteElements(array, new int[]{2, 4, 6});
}
    */
    /*
       public static final void main(String[] args)
       {
           computeShortFromInt((int)65000);
       }
    */

    /*
        static private void computeShort(int b1, int b2)
        {
            byte[] bytes = new byte[] { b1, b2 };
            short result = StsMath.convertShortBytes(bytes, 0);
            System.out.println("bytes: " + b1 + " " + b2 + " = " + result);
        }
    */

    static private void computeShortFromInt(int i)
    {
        float f = (short) i;
        System.out.println("i: " + i + " f: " + f);
    }

    /*
      public static final void main(String[] args)
      {
         float result;
         byte[] bytes = new byte[] { 0, 0, 0, 0 };
         result = StsMath.convertIBMFloatBytes(bytes, false);
         System.out.println("result old " + result);
         result = StsMath.convertIBMFloatBytes(bytes, true);
         System.out.println("result new " + result);

         bytes = new byte[] {24, -31, 0, 0};
         result = StsMath.convertIBMFloatBytes(bytes, false);
         System.out.println("result old " + result);
         result = StsMath.convertIBMFloatBytes(bytes, false);
         System.out.println("result new " + result);
      }
    */
    /*
        public static final void main(String[] args)
        {
            testConvert(12345e6f);
            testConvert(-12345e6f);
            testConvert(12345e-6f);
            testConvert(-12345e-6f);
        }
    */
    private static void testConvert(float ieee)
    {
        int ibmBits = convertIeeeFloatToIbmBits(ieee, false);
        float ieeeFloat = convertIBMFloatIntBits(ibmBits);
        System.out.println(ieee + " converted to " + ieeeFloat);
    }

    /**
     * divides all values in floatArray by absolute max value in floatArray
     *
     * @param floatArray
     * @return floatArray normalized by max value or null if floatArray is null
     */
    public static float[][][] normalizeAmplitude(float[][][] floatArray)
    {
        if (floatArray == null) return null;
        float maxAmp = StsMath.absMax(floatArray);
        float[][][] normalizedArray = floatArray.clone();
        for (int i = 0; i < floatArray.length; i++)
        {
            for (int j = 0; j < floatArray[0].length; j++)
            {
                StsMath.normalizeAmplitude(normalizedArray[i][j], maxAmp);
            }
        }
        return normalizedArray;
    }

    /**
     * finds absolute maximum value in floatArray
     *
     * @param floatArray
     * @return maximum amplitude or zero if floatArray is null
     */
    private static float absMax(float[][][] floatArray)
    {
        float maxVal = 0.0f;
        if (floatArray == null) return maxVal;
        for (int i = 0; i < floatArray.length; i++)
        {
            for (int j = 0; j < floatArray[0].length; j++)
            {
                for (int k = 0; k < floatArray[0][0].length; k++)
                {
                    maxVal = Math.max(Math.abs(floatArray[i][j][k]), maxVal);
                }
            }
        }
        return maxVal;
    }

    /**
     * divides all values in float array by the RMS amplitude of that array.
     * Ignores zeroed amplitudes (mute zone)
     *
     * @return normalized by RMS value or null if floatArray is null
     */
    public static float[][][] normalizeAmplitudeRMS(float[][][] array)
    {
        if (array == null || array[0] == null || array[0][0] == null) return null;
        double rmsAmp = StsMath.rmsIgnoreZero(array);
        for (int i = 0; i < array.length; i++)
        {
            for (int j = 0; j < array[0].length; j++)
            {
                StsMath.normalizeAmplitude(array[i][j], rmsAmp);
            }
        }
        return array;
    }

    /**
     * returns root-mean-square of an array.
     * Ignores zeroes.
     *
     * @param array
     * @return rms of array or zero for null array
     */
    public static final double rmsIgnoreZero(float[][][] array)
    {
        if (array == null || array[0] == null || array[0][0] == null) return 0;
        float[][][] squareArray = StsMath.square(array);
        double mean = StsMath.meanIgnoreZero(squareArray);
        return Math.sqrt(mean);
    }

    public static final double meanIgnoreZero(float[][][] array)
    {
        if (array == null || array[0] == null || array[0][0] == null) return 0;
        double mean = 0;
        int size = 0;
        for (int i = 0; i < array.length; i++)
            for (int j = 0; j < array[i].length; j++)
                for (int k = 0; k < array[i][j].length; k++)
                    if (Math.abs(array[i][j][k]) > 0.0)
                    {
                        mean += array[i][j][k];
                        size++;
                    }
        return mean / size;
    }

    public static final double rmsIgnoreZero(float[] array)
    {
        if (array == null) return 0;
        double mean = 0;
        int size = 0;
        for (int n = 0; n < array.length; n++)
        {
            float value = array[n];
            if (value == 0.0f) continue;
            mean += value * value;
            size++;
        }
        if (size == 0) return 0.0;
        return Math.sqrt(mean / size);
    }


    /**
     * returns the input array squared
     *
     * @param array
     * @return
     */
    public static float[][][] square(float[][][] array)
    {
        float[][][] output = new float[array.length][array[0].length][array[0][0].length];
        for (int i = 0; i < array.length; i++)
            for (int j = 0; j < array[i].length; j++)
                for (int k = 0; k < array[i][j].length; k++)
                    output[i][j][k] = array[i][j][k] * array[i][j][k];
        return output;
    }

    /**
     * divides all amplitudes in floatArray by maxAmp
     * <p/>
     * if maxAmp == 0, does nothing.
     *
     * @param floatArray
     * @param maxAmp
     * @return
     */
    public static float[][][] normalizeAmplitude(float[][][] floatArray, float maxAmp)
    {
        if (maxAmp == 0) return floatArray;
        if (floatArray == null) return null;
        for (int i = 0; i < floatArray.length; i++)
        {
            for (int j = 0; j < floatArray[0].length; j++)
            {
                StsMath.normalizeAmplitude(floatArray[i][j], maxAmp);
            }
        }
        return floatArray;
    }

    public static double computeNormalizedRMSAmplitude(float[][] traces)
    {
        int nTraces = traces.length;
        double sum = 0;
        int nValues = 0;
        for (int t = 0; t < nTraces; t++)
        {
            if (traces[t] == null) continue;
            int nSamples = traces[t].length;
            for (int n = 0; n < nSamples; n++)
            {
                float value = traces[t][n];
                if (value != 0.0f)
                {
                    sum += value * value;
                    nValues++;
                }
            }
        }
        return Math.sqrt(sum / nValues);
    }

    public static float[] normalizeAmplitudeRMS(float[] values)
    {
        if (values == null) return null;
        int nSamples = values.length;
        double sum = 0;
        for (int n = 0; n < nSamples; n++)
        {
            double value = values[n];
            sum += value * value;
        }
        double recipAmp = Math.sqrt(sum / nSamples);
        for (int n = 0; n < nSamples; n++)
        {
            values[n] *= recipAmp;
        }
        return values;
    }

    public static float[] arrayMult(float[] array, double factor)
    {
        if (array == null) return null;
        float[] output = new float[array.length];
        for (int i = 0; i < array.length; i++)
        {
            output[i] = (float) (array[i] * factor);
        }
        return output;
    }

    public static double min(double[] vals)
    {
        if (vals == null || vals.length == 0) return 0.0;
        double min = vals[0];
        for (int i = 0; i < vals.length; i++)
        {
            min = Math.min(min, vals[i]);
        }
        return min;
    }

    public static double[] abs(double[] signedArray)
    {
        if (signedArray == null) return null;
        double[] absArray = StsMath.copy(signedArray);
        for (int i = 0; i < absArray.length; i++)
            absArray[i] = Math.abs(absArray[i]);
        return absArray;
    }

    /**
     * finds minimum and maximum of array
     * returns {0,0} for zero--length array
     *
     * @param array
     * @return float array {min, max} or null
     */
    public static float[] minMax(float[] array)
    {
        if (array == null) return null;
        if (array.length == 0) return new float[]{0, 0};
        float min = array[0];
        float max = array[0];
        for (int i = 1; i < array.length; i++)
        {
            min = Math.min(min, array[i]);
            max = Math.max(max, array[i]);
        }
        return new float[]{min, max};
    }

    public static float[] minMax(float[][][] array)
    {
        if (array == null) return null;
        if (array.length == 0) return new float[]{0, 0};
        float min = Float.MAX_VALUE;
        float max = 0.0f - min;
        for (int i = 0; i < array.length; i++)
        {
            if (array[i] == null) continue;
            for (int j = 0; j < array[i].length; j++)
            {
                float[] minMax = StsMath.minMax(array[i][j]);
                if (minMax != null)
                {
                    min = Math.min(min, minMax[0]);
                    max = Math.max(max, minMax[1]);
                }
            }
        }
        return new float[]{min, max};
    }
    public static final double length(double x, double y)
    {
        return Math.sqrt(x*x + y*y);
    }

    public static final double lengthSq(double x, double y)
    {
        return x*x + y*y;
    }
    static final public byte floatMinusOneToOneToUnsignedByte254(float value)
    {
        int i = Math.round(127*(value + 1.0f));
        if (i >= 255) i = 254;
        if (i < 0) i = 0;
        return (byte) i;
    }


}
/*
  public static final void main(String[] args)
  {
 float result;
 byte[] bytes = new byte[] { 0, 0, -96, 65};
 result = StsMath.convertIBMFloatBytes(bytes, false);
 System.out.println("result old " + result);
 result = StsMath.IBMFloatBigEndian(bytes);
 System.out.println("result new " + result);

 bytes = new byte[] {7, -76, -128, -7};
 result = StsMath.convertIBMFloatBytes(bytes, false);
 System.out.println("result old " + result);
 result = StsMath.IBMFloatBigEndian(bytes);
 System.out.println("result new " + result);
  }
 */

/*
 public static final void main(String[] args)
 {
  try
  {
   Comparator comparator = new TestComparator();
   StsPoint[] points = new StsPoint[5];
   for (int n = 0; n < 5; n++)
   {
 points[n] = new StsPoint(n, n);
   }
   StsPoint point;
   point = new StsPoint( -1, -1);
   points = (StsPoint[]) StsMath.arrayAddSortedElement(points, point, comparator);
   point = new StsPoint(.5f, .5f);
   points = (StsPoint[]) StsMath.arrayAddSortedElement(points, point, comparator);
   point = new StsPoint(3.5f, 3.5f);
   points = (StsPoint[]) StsMath.arrayAddSortedElement(points, point, comparator);
   point = new StsPoint(5.5f, 5.5f);
   points = (StsPoint[]) StsMath.arrayAddSortedElement(points, point, comparator);
   point = new StsPoint(.5f, .5f);
   points = (StsPoint[]) StsMath.arrayAddSortedElement(points, point, comparator);
  }
  catch (Exception e)
  {
   e.printStackTrace();
  }
 }
 */
/*
 float f1 = 1.0f;
 float f2 = 1.00001f;
 boolean same = sameAsTol(f1, f2);
 f2 = 1.00001f;
 same = sameAsTol(f1, f2);
 System.out.println(f1 + " " + f2 + " " + same);
 f2 = 1.000001f;
 same = sameAsTol(f1, f2);
 System.out.println(f1 + " " + f2 + " " + same);
 f2 = 1.0000001f;
 same = sameAsTol(f1, f2);
 System.out.println(f1 + " " + f2 + " " + same);
 f1 = 1000000f;
 f2 = 1000100f;
 same = sameAsTol(f1, f2);
 System.out.println(f1 + " " + f2 + " " + same);
 f2 = 1000010f;
 same = sameAsTol(f1, f2);
 System.out.println(f1 + " " + f2 + " " + same);
 f2 = 1000001f;
 same = sameAsTol(f1, f2);
 System.out.println(f1 + " " + f2 + " " + same);

 double d1 = 1.0;
 double d2 = 1.00000000001;
 same = sameAsTol(d1, d2);
 d2 = 1.00000000001;
 same = sameAsTol(d1, d2);
 System.out.println(d1 + " " + d2 + " " + same);
 d2 = 1.000000000001;
 same = sameAsTol(d1, d2);
 System.out.println(d1 + " " + d2 + " " + same);
 d2 = 1.0000000000001;
 same = sameAsTol(d1, d2);
 System.out.println(d1 + " " + d2 + " " + same);
 d1 = 1000000000000000.0;
 d2 = 1000000000000100.0;
 same = sameAsTol(d1, d2);
 System.out.println(d1 + " " + d2 + " " + same);
 d2 = 1000000000000010.0;
 same = sameAsTol(d1, d2);
 System.out.println(d1 + " " + d2 + " " + same);
 d2 = 1000000000000001.0;
 same = sameAsTol(d1, d2);
 System.out.println(d1 + " " + d2 + " " + same);
 */

/*
   float x;
   System.out.println(simpsonCoefs(0, 6));
   System.out.println(simpsonCoefs(1, 6));
   System.out.println(simpsonCoefs(2, 6));
   System.out.println(simpsonCoefs(3, 6));
   System.out.println(simpsonCoefs(4, 6));
   System.out.println(simpsonCoefs(5, 6));

   System.out.println(simpsonCoefs(0, 7));
   System.out.println(simpsonCoefs(1, 7));
   System.out.println(simpsonCoefs(2, 7));
   System.out.println(simpsonCoefs(3, 7));
   System.out.println(simpsonCoefs(4, 7));
   System.out.println(simpsonCoefs(5, 7));
   System.out.println(simpsonCoefs(6, 7));
 */
/*
 System.out.println(unsignedByteToSignedByte( (byte) 0)); // 0 -> 0 -> -128
 System.out.println(unsignedByteToSignedByte( (byte) 127)); // 127 -> 127 -> -1
 System.out.println(unsignedByteToSignedByte( (byte) - 128)); // -128 -> 128 -> 0
 System.out.println(unsignedByteToSignedByte( (byte) - 1)); // -1 -> 255 -> 127
 */
/*
   System.out.println(unsignedByteToSignedInt((byte)0));  // 0 -> 0 -> -128
   System.out.println(unsignedByteToSignedInt((byte)127));  // 127 -> 127 -> -1
   System.out.println(unsignedByteToSignedInt((byte)-128)); // -128 -> 128 -> 0
   System.out.println(unsignedByteToSignedInt((byte)-1));   // -1 -> 255 -> 127
 */
/*
   byte[] bytes = new byte[4];
   bytes[3] = 64;
   bytes[2] = 31;
   bytes[0] = 0;
   bytes[1] = 0;
   int i = intBits(bytes);
   System.out.println(i);
 */
/*
  float value = 10560.044f;
  ByteArrayOutputStream baos = new ByteArrayOutputStream(4);
  DataOutputStream dos = new DataOutputStream(baos);
  dos.writeFloat(value);
  byte[] bytes = baos.toByteArray();
 */
/*
 byte[] bytes = new byte[4];
 bytes[0] = newBytes[3];
 bytes[1] = newBytes[2];
 bytes[2] = newBytes[1];
 bytes[3] = newBytes[0];
 */
/*
  int intBits =  intBits(bytes);
  float newValue = Float.intBitsToFloat(intBits);
  System.out.println("value in: " + value + " value out: " + newValue);
 */

/*
  float max;

  max = findQuadraticMinMax(0.0f, 1.0f, 0.0f, true);
  max = findQuadraticMinMax(0.0f, 1.0f, -1.0f, true);
 */
/*
  float max;
  float v0, v1, s0, s1;
  boolean mainDebug = false;

  s0 = 0.0f;

   s1 = 0.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


   s1 = 1.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


   s1 = -1.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

  s0 = 1.0f;

   s1 = 0.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


   s1 = 1.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


   s1 = -1.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);
  s0 = -1.0f;

   s1 = 0.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


   s1 = 1.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


   s1 = -1.0f;

 v0 = -1.0f;
 v1 = 1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);

 v0 = 1.0f;
 v1 = -1.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);


 v0 = 0.0f;
 v1 = 0.0f;
 max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 if(mainDebug) System.out.print("v0: " + v0 + " v1: " + v1 + " s0: " + s0 + " s1: " + s1);
 if(max == StsParameters.nullValue)
  System.out.println("No max found.");
 else
  System.out.println("Max at: " + max);
 */
/*
  v0 = 0f;
  v1 = 1f;
  s0 = 0f;
  s1 = 4f;
  max = StsMath.findCubicMax(v0, v1, s0, s1, true);

  v0 = 0f;
  v1 = 1f;
  s0 = 1f;
  s1 = 0f;
  max = StsMath.findCubicMax(v0, v1, s0, s1, true);

  v0 = 1f;
  v1 = 0f;
  s0 = 0f;
  s1 = -1f;
  max = StsMath.findCubicMax(v0, v1, s0, s1, true);

  v0 = 17.6f;
  v1 = 20.6f;
  s0 = (20.6f-9.6f)/2;
  s1 = (19.6f-17.6f)/2;
  max = StsMath.findCubicMax(v0, v1, s0, s1, true);

  v0 = v1;
  v1 = 19.6f;
  s0 = s1;
  s1 = (10.6f-20.6f)/2;
  max = StsMath.findCubicMax(v0, v1, s0, s1, true);

  v0 = 0f;
  v1 = 1f;
  s0 = 1f;
  s1 = 1f;
  max = StsMath.findCubicMax(v0, v1, s0, s1, true);
 */

//        max = StsMath.findCubicMax(0.0f, 0.0f, 1.0f, -1.0f, true);
//        max = StsMath.findCubicMax(0.0f, 0.0f, 1.0f, 1.0f, true);
//        max = StsMath.findCubicMax(0.0f, 0.0f, -1.0f, -1.0f, true);

//    niceScaleBold(3508.0, 3492.0, 78, true);
/*
 System.out.println((int)(byte)127);  // 127
 System.out.println((int)(byte)-128); // 128
 System.out.println((int)(byte)-1);   // 255

 System.out.println(signedByteToUnsignedInt((byte)127));  // 127
 System.out.println(signedByteToUnsignedInt((byte)-128)); // 128
 System.out.println(signedByteToUnsignedInt((byte)-1));   // 255

 System.out.println(UnsignedIntToSignedByte(127));  // 127
 System.out.println(UnsignedIntToSignedByte(128)); // -128
 System.out.println(UnsignedIntToSignedByte(255));   // -1
 */
/*
  double[] scale = niceScale2(21050.0, 22175.0, 10, true);
  System.out.println("scale: " + scale[0] + " " + scale[1] + " " + scale[2] + " " + scale[3]);
 */
/*
 com.Sts.Types.StsGridLink[] gridLinks = new com.Sts.Types.StsGridLink[4];
 gridLinks[0] = new com.Sts.Types.StsGridLink();
 gridLinks[1] = new com.Sts.Types.StsGridLink();
 gridLinks[2] = new com.Sts.Types.StsGridLink();
 gridLinks[3] = new com.Sts.Types.StsGridLink();

 gridLinks = (com.Sts.Types.StsGridLink[])arraycopy(gridLinks, 10);
 for(int n = 0; n < 1; n++)
 {
  com.Sts.Types.StsGridLink gridLink = gridLinks[n];
 }
 */
//        double[] scale = niceScale(0.0, 200.0, true);
/*
   float[] pointVector, lineVector;

   lineVector = new float[] { 1.0f, 1.0f };
   pointVector = new float[] { .5f, 1.0f };
   boolean isOnLine = isOnLine(pointVector, lineVector, 0.01f);
 */
/*
 float value = 12347.1f;
 boolean isIt = isExactIntegral(value);
 System.out.println(isIt);
 */

/*
 float[] floats = new float[] { 106.999985f, 106.9999985f };
 int nFloats = floats.length;
 for(int n = 0; n < nFloats; n++)
 {
  Float f = new Float(floats[n]);
  System.out.println("float: " + f.intValue());
  double d = (double)floats[n];
  BigDecimal bd = new BigDecimal(d);
  bd.setScale(0, BigDecimal.ROUND_UP);
  System.out.println("bigDecimal: " + bd.intValue());
 }
 */
/*
  int index;
  float[] values = new float[] { 1, 2, 3, 4, 5, 6, 7 };

  index = arrayIndexAbove(2.5f, values);
  index = arrayIndexBelow(2.5f, values);

  index = arrayIndexAbove(5.0f, values);
  index = arrayIndexBelow(5.0f, values);
 */
/*
  int i;
  boolean isEqual;
  float v[] = new float[] { -4.000001f, -3.999999f, 4.000001f, 3.999999f, 3.99999999f,
 -4.00001f, -3.9999f, 4.0001f, 3.9999f, 3.99999f };

  for(int n = 0; n < v.length; n++)
  {
   i = StsMath.ceilingRoundOff(v[n]);
   isEqual = StsMath.sameAsTol((float)i, v[n]);
   System.out.println("input: " + v[n] + " ceiling: " + i + " are equal: " + isEqual);

   i = StsMath.floorRoundOff(v[n]);
   isEqual = StsMath.sameAsTol((float)i, v[n]);
   System.out.println("input: " + v[n] + " floor: " + i + " are equal: " + isEqual);
  }
 */
/*
 float[] p0 = new float[] { 0.0f, 0.0f, 0.0f };
 float[] v0 = new float[] { 1.0f, 1.0f, 1.0f };
 float[] p1 = new float[] { 0.0f, 1.0f, 1.0f };
 float[] v1 = new float[] { 1.0f, -1.0f, -1.0f };

 float[] i = StsMath.lineIntersectXYFinite(p0, v0, p1, v1);
 System.out.println(" x: " + i[0] + " y: " + i[1] + " z: " + i[2]);
 */

/*
  int i;
  boolean isEqual;
  float v[] = new float[] { -4.000001f, -3.999999f, 4.000001f, 3.999999f, 3.99999999f };

  for(int n = 0; n < v.length; n++)
  {
   i = StsMath.ceiling(v[n]);
   isEqual = StsMath.equals(i, v[n]);
   System.out.println("input: " + v[n] + " ceiling: " + i + " are equal: " + isEqual);
   isEqual = i == v[n];
   System.out.println("input: " + v[n] + " ceiling: " + i + " == : " + isEqual);

   i = StsMath.floor(v[n]);
   isEqual = StsMath.equals(i, v[n]);
   System.out.println("input: " + v[n] + " floor: " + i + " are equal: " + isEqual);
   isEqual = i == v[n];
   System.out.println("input: " + v[n] + " ceiling: " + i + " == : " + isEqual);
  }
 */
/*
 boolean isIntegral;
 float v[] = new float[] { -4.0001f, -3.9999f, 4.0001f, 3.9999f, 3.9999999f };

 for(int n = 0; n < v.length; n++)
 {
  isIntegral = isIntegral(v[n]);
  System.out.println("input: " + v[n] + " output: " + isIntegral);
 }
 */

//        float[] b = binomWeights(5);

//        System.out.println(b[0] + " "+ b[1] + " " + b[2] + " " + b[3] + " " + b[4]+ " " + b[5]);
/*
  float v = -4.0f;
  float d  = 5.0f;
  float r = remainder(v, d);
  System.out.println("input: " + v + " " + d + " output: " + r);
  v = -9.0f;
  d  = 5.0f;
  r = remainder(v, d);
  System.out.println("input: " + v + " " + d + " output: " + r);
 */
