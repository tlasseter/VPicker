
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types.PreStack;

import java.io.Serializable;

import com.Sts.DB.StsSerializable;
import com.Sts.DBTypes.StsSerialize;
import com.Sts.DBTypes.*;

public class StsMicroseismicGatherFile extends StsSerialize implements Cloneable, StsSerializable, Serializable, Comparable
{
	String segyFilename = null;
	long startTime = 0l;
	long endTime = 0l;
	
	transient StsMicroseismicGather gather = null;
	
    public StsMicroseismicGatherFile()
    {
    }

    public StsMicroseismicGatherFile(String file, long start, long end)
    {
        segyFilename = file;
        startTime = start;
        endTime = end;
    }
    
    public int compareTo(Object otherObj)
    {
    	StsMicroseismicGatherFile otherGather = (StsMicroseismicGatherFile)otherObj;
    	if(otherGather.segyFilename.equals(segyFilename))
    		return 0;
    	else
    		return -1;
    }
    
    public long getStartTime() { return startTime; }
    public long getEndTime() { return endTime; }
    public String getSegyFilename() { return segyFilename; }
    
    public StsMicroseismicSuperGather buildSuperGather(StsPreStackMicroseismicSet lineSet)
    {
    	StsMicroseismicSuperGather gather = StsMicroseismicSuperGather.constructor(currentModel, lineSet, this);
    	return gather;
    }
}





