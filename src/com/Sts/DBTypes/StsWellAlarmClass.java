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

public class StsWellAlarmClass extends StsAlarmClass implements StsSerializable
{
   public StsWellAlarmClass()
   {
       userName = "Well Proximity Alarms";
   }

}