
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;

import java.util.*;

public class StsHistogram
{
    private int dataCnt[] = new int[255];
    private float[] histVals = null;
    private int ttlHistogramSamples, minIdx, maxIdx;
    private float dataMin = StsParameters.largeFloat;
    private float dataMax = StsParameters.smallFloat;
    private float dataRangeMin = StsParameters.largeFloat;
    private float dataRangeMax = StsParameters.smallFloat;
    private int indexOfMaxPercent = -1;

    public StsHistogram()
    {
        ttlHistogramSamples = 0;
        minIdx = 0;
        maxIdx = 254;
    }

    public void clearHistogram()
    {
        dataCnt = null;
        histVals = null;
        ttlHistogramSamples = 0;
        dataCnt = new int[255];
        histVals = new float[255];
    }

    public int getTotalSamples() { return ttlHistogramSamples; }
    public float[] getValues() { return histVals; }
    public void setDataMin(float min)
    {
        dataMin = min;
        if(dataMin < dataRangeMin)
            dataMin = dataRangeMin;
    }
    public void setDataMax(float max)
    {
        dataMax = max;
        if(dataMax > dataRangeMax)
            dataMax = dataRangeMax;
    }
    public float getDataMin() { return dataMin; }
    public float getDataMax() { return dataMax; }

    public void setDataRangeMin(float min) { dataRangeMin = min; }
    public void setDataRangeMax(float max) { dataRangeMax = max; }
    public float getDataRangeMin() { return dataRangeMin; }
    public float getDataRangeMax() { return dataRangeMax; }
    public void resetLimits()
    {
        dataMin = dataRangeMin;
        dataMax = dataRangeMax;
        dataRangeMin = StsParameters.largeFloat;
        dataRangeMax = StsParameters.smallFloat;
    }

    public void computeMinMax(float min, float max)
    {
        if((dataMin < min) || (dataMin == StsParameters.largeFloat))
           dataMin = min;
       if((dataMax > max) || (dataMax == StsParameters.smallFloat))
           dataMax = max;
       if(dataRangeMin > min)
           dataRangeMin = min;
       if(dataRangeMax < max)
           dataRangeMax = max;
    }

    public void calculateHistogram()
    {
        for(int i=0; i<255; i++)
            histVals[i] = (float)((float)dataCnt[i]/(float)ttlHistogramSamples)*100.0f;
        return;
    }

    public float getMaxPercentage()
    {
        float maxPercent = 0.0f;
        indexOfMaxPercent = -1;

        if(histVals != null)
        {
            for (int i = 0; i < 255; i++)
            {
                if(histVals[i] > maxPercent)
                {
                    maxPercent = histVals[i];
                    indexOfMaxPercent = i;
                }
            }
        }
        return maxPercent;
    }

    public float getRecommendedVerticalScale()
    {
        if(histVals != null)
        {
            double[] sortVals = StsMath.copyDouble(histVals);
            Arrays.sort(sortVals);
            for(int i=sortVals.length-1; i>0; i--)
            {
                if((sortVals[i]/10.0f) < sortVals[i-1])
                    return (float)sortVals[i-1];
            }
        }
        return 1.0f;
    }

    public int indexOfNull()
    {
        getMaxPercentage();
        return indexOfMaxPercent;
    }

    public void accumulateHistogram(float bindex, int count)
    {
        int index = computeIndex(dataMin, dataMax, bindex);
        dataCnt[index] = dataCnt[index] + count;
        ttlHistogramSamples += count;
    }

    public int computeIndex(float min, float max, float bindex)
    {
        byte bsamp = 0;
        float scale = 254 / (max - min);
        float scaledFloat = (bindex - min)*scale;
        int scaledInt = Math.round(scaledFloat);
        bsamp = unsignedIntToSignedByte254(scaledInt);

        int index = StsMath.signedByteToUnsignedInt(bsamp);
        if(index > 254) index = 254;
        if(index < 0) index = 0;
        return index;
    }

    /** converts an unsigned int to a signedByte value between/including 0 to 254 */
    public static final byte unsignedIntToSignedByte254(int i)
    {
        if(i >= 255) i = 254;
        if(i < 0) i = 0;
        return (byte)i;
    }
}





