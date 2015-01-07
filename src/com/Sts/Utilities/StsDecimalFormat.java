
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Utilities;

import java.text.*;

public class StsDecimalFormat extends DecimalFormat
{
    int nSignificantDigits = 0;
    int nIntegerDigits = 0;
    int nFractionalDigits = 0;

    public StsDecimalFormat()
    {
        setGroupingUsed(false);
    }

    public StsDecimalFormat(String pattern)
    {
        super(pattern);
        setGroupingUsed(false);
    }

    public StsDecimalFormat(int nDigits)
    {
        setSignificantDigits(nDigits);
        setGroupingUsed(false);
    }
    public void setSignificantDigits(int digits)
    {
        this.nSignificantDigits = digits;
        if(digits == 0) return;
		setMaximumIntegerDigits(digits);
        setMaximumFractionDigits(digits);
    }

    public void formatValueRange(double startValue, double endValue)
    {
        formatValue(startValue);
        int nFractionalDigitsStart = nFractionalDigits;
        int nIntegerDigitsStart = nIntegerDigits;
        formatValue(endValue);
        setMaximumIntegerDigits(Math.max(nIntegerDigits, nIntegerDigitsStart));
        setMaximumFractionDigits(Math.max(nFractionalDigits, nFractionalDigitsStart));
        setMinimumFractionDigits(Math.max(nFractionalDigits, nFractionalDigitsStart));    
    }

    public String formatValue(double value)
    {
        double absValue = Math.abs(value);
        double log10 = Math.log10(absValue);
        if(log10 >= 0)
        {
            nIntegerDigits = (int)log10 + 1;
            if(nIntegerDigits >= nSignificantDigits)
                nFractionalDigits = 0;
            else
                nFractionalDigits = nSignificantDigits - nIntegerDigits;
        }
        else
        {
            nIntegerDigits = 0;
            nFractionalDigits = StsMath.ceiling(-log10) + 1;
            if(nFractionalDigits < nSignificantDigits)
                nFractionalDigits = nSignificantDigits;
        }
        setMaximumIntegerDigits(nIntegerDigits);
        setMaximumFractionDigits(nFractionalDigits);
        return format(value);
    }

    public String stripLeadingZeroes(float value)
    {
        return stripLeadingZeroes(format(value));
    }

    public String stripLeadingZeroes(String fmtString)
    {
        if (fmtString == null) return null;
        int pos = 0;
        StringBuffer fmt = new StringBuffer(fmtString);
        char sign = ' ';
        if (fmt.charAt(0) == '+' || fmt.charAt(0) == '-')
        {
            sign = fmt.charAt(0);
            fmt.setCharAt(0, ' ');
            pos++;
        }
        int length = fmt.length();
        for (int i=pos; i<length-1; i++)
        {
            if (fmt.charAt(i) == '0' && fmt.charAt(i+1) != '.' && i != length-2)
            {
                fmt.setCharAt(i, ' ');
            }
            else
            {
                if (i > 0) fmt.setCharAt(i-1, sign);
                break;
            }
        }
        return fmt.toString();
    }

    static public void main(String[] args)
    {
        StsDecimalFormat format = new StsDecimalFormat(8);
        testFormat(format, 0.000011f);
        testFormat(format, 1111.1111f);
        testFormat(format, 0.111111f);
        testFormat(format, 111111.11f);

        testFormat(format, 0.000011);
        testFormat(format, 1111.1111);
        testFormat(format, 0.111111);
        testFormat(format, 111111.11);
    }

    static private void testFormat(StsDecimalFormat format, float value)
    {
        String valueString = format.formatValue(value);
        System.out.println("float value: " + value + " string: " + valueString);
    }


    static private void testFormat(StsDecimalFormat format, double value)
    {
        String valueString = format.formatValue(value);
        System.out.println("double value: " + value + " string: " + valueString );
    }
}