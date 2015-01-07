//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTreeObjectI;

public class StsSensorPartitionFilter extends StsSensorFilter implements StsTreeObjectI
{
    long[] starts = null;
    long[] ends = null;
    int[] colorIdxs = null;

    /** default constructor */
    public StsSensorPartitionFilter()
    {
        super(PARTITION, true);
    }
    public StsSensorPartitionFilter(String name)
    {
        super(PARTITION, true);
        this.name = name + getIndex();
    }

    /**
     * Get the filter class containing this object
     * @return
     */
    static public StsSensorPartitionFilterClass getSensorQualifyFilterClass()
    {
        return (StsSensorPartitionFilterClass)currentModel.getCreateStsClass(StsSensorPartitionFilter.class);
    }

    public void addTimeRanges(long[] startTimes, long[] endTimes, int[] colorIdx)
    {
        starts = startTimes;
        ends = endTimes;
        colorIdxs = colorIdx;
    }

    public boolean filter(StsSensor sensor)
    {
        long[] times = sensor.getTimeCurve(0).getTimeVectorLongs();
        int[] passed = new int[times.length];
        sensor.setClustering(true);
        for(int i=0; i<times.length; i++)
        {
            passed[i] = -1;
            for(int j=0; j<ends.length; j++)
            {
                if((times[i] >= starts[j]) && (times[i] <= ends[j]))
                {
                    passed[i] = colorIdxs[j];
                    break;
                }
            }
        }
        sensor.setClusters(passed);
        currentModel.viewObjectRepaint(this, sensor);
        return true;
    }

}