package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.StsException;

public class StsMonitorClass extends StsClass implements StsSerializable
{
   boolean enable = true;
   int[] numPolls;
   int[] acceptedValues;
   long[] lastPollTime;
   long[] lastFileSize;
   
   public StsMonitorClass()
   {
       userName = "Realtime Data Monitors";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]  
       {  
    	   new StsBooleanFieldBean(this, "enable", "Enable Monitoring:"), 
       };
       defaultFields = new StsFieldBean[]  {   };
   }

   public boolean setCurrentObject(StsObject object)
   {
        if(object == null) return false;
		if (currentObject == object)
			return false;
		currentObject = object;
        ((StsMonitor)object).updateProperties();
		return super.setCurrentObject(object);
   }

   public void setEnable(boolean enable)
   {
       if(this.enable == enable) return;
       this.enable = enable;
//       setDisplayField("enable", enable);
       currentModel.win3dDisplayAll();
   }
   public boolean getEnable() {	return enable; }

    public void close()
    {
        int nValues = getSize();
        numPolls = new int[nValues];
        acceptedValues = new int[nValues];
        lastPollTime = new long[nValues];
        lastFileSize = new long[nValues];
        for(int n = 0; n < nValues; n++)
        {
            StsMonitor monitor = (StsMonitor)getElement(n);
            if(monitor != null)
            {
                numPolls[n] = monitor.numPolls;
                acceptedValues[n] = monitor.acceptedChanges;
                lastPollTime[n] = monitor.lastPollTime;
                lastFileSize[n] = monitor.lastFileSize;
            }
        }
    }

    public void finalInitialize()
    {
        int nValues = getSize();
        if(numPolls == null) return;
        int length = numPolls.length;
        if(length != nValues)
            StsException.systemError(this, "finalInitialize", "Saved numPolls.length " + length + " doesn't agree with number of StsMonitor instances " + nValues);

        nValues = Math.min(nValues, length);
        for(int n = 0; n < nValues; n++)
        {
            StsMonitor monitor = (StsMonitor)getElement(n);
            if(monitor != null)
            {
                monitor.numPolls = numPolls[n];
                monitor.acceptedChanges = acceptedValues[n];
                monitor.lastPollTime = lastPollTime[n];
                monitor.lastFileSize = lastFileSize[n];
            }
        }
    }

}