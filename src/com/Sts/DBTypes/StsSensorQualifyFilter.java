//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.Utilities.StsMath;

public class StsSensorQualifyFilter extends StsSensorFilter implements StsTreeObjectI
{
    String[] curveNames = new String[0];
    float[] lows = new float[0];
    float[] highs = new float[0];

    /** default constructor */
    public StsSensorQualifyFilter()
    {
        super(QUALIFY, true);
    }
    public StsSensorQualifyFilter(String name)
    {
        super(QUALIFY, true);
        this.name = name + getIndex();
    }

     /**
     * Get the filter class containing this object
     * @return
     */
    static public StsSensorQualifyFilterClass getSensorQualifyFilterClass()
    {
        return (StsSensorQualifyFilterClass)currentModel.getCreateStsClass(StsSensorQualifyFilter.class);
    }

    public boolean addRange(String curveName, float low, float high)
    {
        if(curveName == null) return false;
        if(low > high) return false;

        curveNames = (String[]) StsMath.arrayAddElement(curveNames, curveName);
        float[] newLows = new float[lows.length+1];
        float[] newHighs = new float[highs.length+1];
        for(int i=0; i<lows.length; i++)
        {
            newLows[i] = lows[i];
            newHighs[i] = highs[i];
        }
        newLows[lows.length] = low;
        newHighs[highs.length] = high;

        lows = newLows;
        highs = newHighs;
        return true;
    }
    
    public boolean filter(StsSensor sensor)
    {
        if(!enable) return true;

        sensor.setClustering(true);
        for(int i=0; i<curveNames.length; i++)
        {
            StsTimeCurve curve = sensor.getTimeCurve(curveNames[i]);
            if(curve != null)
                sensor.setAttributeLimits(curve, lows[i], highs[i]);
        }
        currentModel.viewObjectRepaint(this, sensor);        
        return true;
    }
}
