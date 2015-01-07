//  JAVA ALGORITHMS
//  ---------------
//  Copyright 1997 Scott Robert Ladd
//  All rights reserved
//
//  This software source code is sold as a component of the
//  book JAVA ALGORITHMS, written by Scott Robert Ladd and
//  published by McGraw-Hill, Inc. Please read the LICENSE
//  AGREEMENT and DISCLAIMER OF WARRANTY printed in the book.
//
//  You may freely compile and link this source code into your 
//  non-commercial software programs, providing that you do
//  not redistribute the source code or object code derived
//  therefrom. If you want to use this source code in a
//  commercial application, you must obtain written permission
//  by contacting:
//
//      Scott Robert Ladd
//      P.O. Box 617
//      Silverton, Colorado
//      81433-0617 USA
//
//  This software is sold "as is" without warranty of any kind.

package com.Sts.Utilities.Sort;

//
//
// Translator
//
//
public class Translator 
{

    public static Short [] toShort(short [] a)
    {
        Short [] result = new Short [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = new Short(a[n]);

        return result;
    }

    public static short [] fromShort(Short [] a)
    {
        short [] result = new short [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = a[n].shortValue();

        return result;
    }

    public static Integer [] toInteger(int [] a)
    {
        Integer [] result = new Integer [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = new Integer(a[n]);

        return result;
    }

    public static int [] fromInteger(Integer [] a)
    {
        int [] result = new int [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = a[n].intValue();

        return result;
    }

    public static Long [] toLong(long [] a)
    {
        Long [] result = new Long [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = new Long(a[n]);

        return result;
    }

    public static long [] fromLong(Long [] a)
    {
        long [] result = new long [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = a[n].longValue();

        return result;
    }

    public static Float [] toFloat(float [] a)
    {
        Float [] result = new Float [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = new Float(a[n]);

        return result;
    }

    public static float [] fromFloat(Float [] a)
    {
        float [] result = new float [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = a[n].floatValue();

        return result;
    }

    public static Double [] toDouble(double [] a)
    {
        Double [] result = new Double [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = new Double(a[n]);

        return result;
    }

    public static double [] fromDouble(Double [] a)
    {
        double [] result = new double [a.length];

        for (int n = 0; n < a.length; ++n)
            result[n] = a[n].doubleValue();

        return result;
    }
}

