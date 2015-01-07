package com.Sts.Interfaces;

/**
 * Title:        Workflow development
 * Description:  Interface for any class that reacts to real-time data.
 * Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author S.A. Jackson
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.Types.*;

public interface StsMonitorable
{
	// Supplied the source name and location
	public int addNewData(String sourceLocation, byte sourceType, long lastPollTime, boolean compute, boolean reload, boolean replace);
	
	// attValues and attNames have all attributes including x, y, z, ...
	public int addNewData(double[] attValues, long time, String[] attNames);
	
	// object must be same type as class containing method.
	public int addNewData(StsObject object);
	
	// point with all attributes. Names are in attNames array including x,y,z...
	public int addNewData(StsPoint point, long time, String[] attNames);

	// get the alarms associated with the monitor objects
	public boolean hasAlarms();

    // get alarms
    public StsAlarm[] getAlarms();

    // add alarm
    public boolean addAlarm(StsAlarm alarm);
}
