package com.Sts.DBTypes;

import com.Sts.UI.Beans.StsBooleanFieldBean;
import com.Sts.UI.Beans.StsFieldBean;

import java.util.Iterator;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSensorPartitionFilterClass extends StsSensorFilterClass
{
   boolean enable = true;

   public StsSensorPartitionFilterClass()
   {
       userName = "Sensor Parition Filters";
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
           StsSensorPartitionFilter filter = (StsSensorPartitionFilter)iter.next();
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
           StsSensorPartitionFilter sensorFilter = (StsSensorPartitionFilter)iter.next();
           sensorFilter.filter(sensor);
       }
       return true;   }
}