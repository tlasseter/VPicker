package com.Sts.DBTypes;

import com.Sts.Interfaces.StsClassTimeDisplayable;
import com.Sts.Interfaces.StsClassTimeSeriesDisplayable;
import com.Sts.UI.Beans.*;

import java.util.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsStaticSensorClass extends StsSensorClass implements StsClassTimeDisplayable, StsClassTimeSeriesDisplayable
{
   boolean enable = true;

   public StsStaticSensorClass()
   {
       userName = "Fixed Location Sensors";
   }

   public void initializeFields()
   {
       displayFields = new StsFieldBean[]
       {
           new StsBooleanFieldBean(this, "enable", "Enable"),
           new StsBooleanFieldBean(this, "displayNames", "Names"),
           new StsBooleanFieldBean(this, "enableTime", "Enable Time"),
           new StsBooleanFieldBean(this, "enableSound", "Sound Accent"),
           new StsBooleanFieldBean(this, "enableAccent", "Graphic Accent"),
           new StsBooleanFieldBean(this, "enableGoTo", "Go To New Data"),
           new StsIntFieldBean(this, "goToOffset", 1, 100, "Go To Offset (*scale):", true),
       };

       defaultFields = new StsFieldBean[]
       {
           new StsIntFieldBean(this, "defaultSize", 1, 500, "Size:", true),
           new StsColorComboBoxFieldBean(this, "defaultSensorColor", "Color:", StsColor.colors32),
           new StsColorComboBoxFieldBean(this, "accentColor", "Accent Color:", StsColor.colors32),
           new StsComboBoxFieldBean(this, "defaultSymbolString", "Symbol:", StsSensor.SYMBOL_TYPE_STRINGS),
           new StsComboBoxFieldBean(this, "defaultSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums)
       };
   }

   public void setEnableStatic(boolean val)
   {
       if(enable == val) return;
       enable = val;
       Iterator iter = getObjectIterator();
       while(iter.hasNext())
       {
           StsStaticSensor sensor = (StsStaticSensor)iter.next();
           sensor.setIsVisible(val);
       }
   }

   public boolean getEnableStatic()
   {
       return enable;
   }

   public StsObject[] getSensors()
   {
       return getObjectList();
   }

   public void setEnable(boolean val)
   {
       if(enable == val) return;
       enable = val;
       Iterator iter = getObjectIterator();
       while(iter.hasNext())
       {
           StsStaticSensor sensor = (StsStaticSensor)iter.next();
           sensor.setIsVisible(val);
       }
   }

   public boolean getEnable() { return enable; }
}