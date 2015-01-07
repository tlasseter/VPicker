//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.Types.StsPoint;
import com.Sts.Utilities.StsMath;

public class StsSensorProximityFilter extends StsSensorFilter implements StsTreeObjectI
{
    StsWell well = null;
    float distance = 0.0f;
    boolean insideLimits = true;
    boolean ignoreVertical = false;
    
    /** default constructor */
    public StsSensorProximityFilter()
    {
        super(PROXIMITY, true);
    }
    public StsSensorProximityFilter(String name)
    {
        super(PROXIMITY, true);
        this.name = name + getIndex();
    }

    /**
     * Set well to compute distance from
     */
    public void setWell(StsWell well)
    {
        this.well = well;
    }

    /**
     * Set distance criteria
     * @param distance from well
     * @param inside distance from well or outside it
     */
    public void setDistanceCriteria(float distance, boolean inside, boolean ignoreVert)
    {
        this.distance = distance;
        this.insideLimits = inside;
        this.ignoreVertical = ignoreVert;
    }

    /**
     * Set distance criteria
     * @param distance from well
     * @param inside distance from well or outside it
     */
    public void setDistanceCriteria(float distance, boolean inside)
    {
        this.distance = distance;
        this.insideLimits = inside;
        this.ignoreVertical = false;
    }

    /**
     * Set distance criteria
     * @param distance from well
     */
    public void setDistanceCriteria(float distance)
    {
        this.distance = distance;
        this.insideLimits = true;
        this.ignoreVertical = false;
    }
     /**
     * Get the filter class containing this object
     * @return
     */
    static public StsSensorProximityFilterClass getSensorQualifyFilterClass()
    {
        return (StsSensorProximityFilterClass)currentModel.getCreateStsClass(StsSensorProximityFilter.class);
    }


    public boolean filter(StsSensor sensor)
    {
        if(!enable) return true;

        // If the sensor has no positional curves it cannot be filtered
        if(sensor.getTimeCurve(StsLogVector.types[StsLogVector.X]) == null)
            return false;

        sensor.setClustering(true);
        long[] times = sensor.getTimeCurve(StsLogVector.types[StsLogVector.X]).getTimeVectorLongs();
        int[] passed = new int[times.length];
        StsPoint sensorPt = null;
        float dist;
        for(int i=0; i<times.length; i++)
        {
  			passed[i] = -1;
			dist = -1.0f;
			sensorPt = ((StsDynamicSensor)sensor).getZPoint(i);
			StsPoint pt = null;
			if (!ignoreVertical)
				pt = StsMath.getNearestPointOnLine(sensorPt, well.getRotatedPoints(), 3);
			else
				pt = StsMath.getNearestPointOnLine(sensorPt, well.getRotatedPoints(), 2);

			dist = pt.distance(sensorPt);
			if (((dist < distance) && insideLimits)|| ((dist > distance) && !insideLimits))
				passed[i] = 99;
        }
        sensor.setClusters(passed);        
        currentModel.viewObjectRepaint(this, sensor);
        return true;
    }

}