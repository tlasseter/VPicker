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

public class StsAlarmClass extends StsClass implements StsSerializable
{
   public StsAlarmClass()
   {
       userName = "Alarms";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
    	   new StsBooleanFieldBean(this, "isVisible", "Enable Alarms"),
       };
   }
}