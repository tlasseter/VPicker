package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.StsSerializable;
import com.Sts.UI.Beans.StsBooleanFieldBean;
import com.Sts.UI.Beans.StsFieldBean;
import com.Sts.Utilities.StsMath;

public class StsSensorFilterClass extends StsClass implements StsSerializable
{
   public StsSensorFilterClass()
   {
       userName = "Sensor Filters";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
    	   new StsBooleanFieldBean(this, "isVisible", "Enable Filters"),
       };
   }

   public boolean filterSensor(StsSensor sensor)
   {
       if(!isVisible) return true;

       if(!((StsSensorQualifyFilterClass)currentModel.getStsClass("com.Sts.DBTypes.StsSensorQualifyFilter")).filterSensor(sensor))
            return false;
       if(!((StsSensorPartitionFilterClass)currentModel.getStsClass("com.Sts.DBTypes.StsSensorPartitionFilter")).filterSensor(sensor))
            return false;
       if(!((StsSensorProximityFilterClass)currentModel.getStsClass("com.Sts.DBTypes.StsSensorProximityFilter")).filterSensor(sensor))
            return false;
       return true;

   }

    public StsObject[] getFilters()
    {
       return getObjectList();
    }

    public StsObject[] getAllFilters()
    {
       StsObject[] objects = ((StsSensorQualifyFilterClass)currentModel.getCreateStsClass("com.Sts.DBTypes.StsSensorQualifyFilter")).getFilters();
       objects = (StsObject[])StsMath.arrayAddArray(objects, ((StsSensorPartitionFilterClass)currentModel.getCreateStsClass("com.Sts.DBTypes.StsSensorPartitionFilter")).getFilters());
       objects = (StsObject[])StsMath.arrayAddArray(objects, ((StsSensorProximityFilterClass)currentModel.getCreateStsClass("com.Sts.DBTypes.StsSensorProximityFilter")).getFilters());

       return objects;
    }
}