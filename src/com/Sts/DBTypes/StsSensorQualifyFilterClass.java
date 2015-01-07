package com.Sts.DBTypes;

import com.Sts.Interfaces.StsClassTimeDisplayable;
import com.Sts.Interfaces.StsClassTimeSeriesDisplayable;
import com.Sts.UI.Beans.*;

import java.util.Iterator;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSensorQualifyFilterClass extends StsSensorFilterClass
{
    boolean enable = true;

   public StsSensorQualifyFilterClass()
   {
       userName = "Sensor Qualify Filters";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
           new StsBooleanFieldBean(this, "enable", "Enable"),
       };

       defaultFields = new StsFieldBean[] {  };
   }

   public void setEnable(boolean val)
   {
       if(enable == val) return;
       enable = val;
       Iterator iter = getObjectIterator();
       while(iter.hasNext())
       {
           StsSensorQualifyFilter filter = (StsSensorQualifyFilter)iter.next();
           filter.setIsVisible(val);
       }
   }

   public boolean getEnable() { return enable; }

   public boolean filterSensor(StsSensor sensor)
   {
       if(!enable) return true;

       sensor.setClustering(true);

       Iterator iter = getObjectIterator();
       while(iter.hasNext())
       {
           StsSensorQualifyFilter sensorFilter = (StsSensorQualifyFilter)iter.next();
           sensorFilter.filter(sensor);
       }
       return true;
   }
}