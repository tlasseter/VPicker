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

public class StsSensorProximityFilterClass extends StsSensorFilterClass
{
   boolean enable = true;

   public StsSensorProximityFilterClass()
   {
       userName = "Sensor Proximity Filters";
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
           StsSensorProximityFilter filter = (StsSensorProximityFilter)iter.next();
           filter.setIsVisible(val);
       }
   }

   public boolean getEnable() { return enable; }

   public boolean filterSensor(StsSensor sensor)
   {
       return true;
   }
}