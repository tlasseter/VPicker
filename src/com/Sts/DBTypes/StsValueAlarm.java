//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.UI.Beans.*;
import com.Sts.UI.Sounds.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;
import com.Sts.Interfaces.*;

public class StsValueAlarm extends StsAlarm implements StsTreeObjectI
 {
     // Event symbols
     float value = StsParameters.nullValue;
     float threshold = 0.0f;  // five percent threshold
     boolean inside = true;   // Inside equals true, value inside range triggers alarm

     static public StsFieldBean[] displayFields = null;
     static public StsFieldBean[] propertyFields = null;
     static protected StsObjectPanel objectPanel = null;
     /** default constructor */
     public StsValueAlarm()
     {
         super();
     }

     public StsValueAlarm(String alarmName, String sound, float value, float threshold, boolean inside)
     {
         super(sound);
         this.value = value;
         this.threshold = threshold;
         this.inside = inside;
         this.name = alarmName;
     }

     static public StsValueAlarm nullAlarmConstructor(String name)
     {
         return new StsValueAlarm("Value_Alarm", null, StsParameters.nullValue, 0.0f, true);
     }

     public StsFieldBean[] getDisplayFields()
     {
         try
         {
             if (displayFields == null)
             {
                 displayFields = new StsFieldBean[]
                 {
                     new StsBooleanFieldBean(this, "Enable", "Enable Alarm:"),
                     new StsComboBoxFieldBean(this, "Sound", "Sound:", StsSound.sounds),
                     new StsComboBoxFieldBean(this, "QueueString", "Type:", QUEUE_TYPE_STRINGS)
                 };
             }
             return displayFields;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "getDisplayFields", e);
             return null;
         }
     }

     public StsFieldBean[] getPropertyFields()
     {
         try
         {
             if (propertyFields == null)
             {
                 propertyFields = new StsFieldBean[]
                 {
                     new StsBooleanFieldBean(StsValueAlarm.class, "Inside", "Inside Range", true),
                     new StsFloatFieldBean(StsValueAlarm.class, "Threshold", true, "Percentage:"),
                     new StsFloatFieldBean(StsValueAlarm.class, "Value", true, "Value:")
                 };
             }
             return propertyFields;
         }
         catch (Exception e)
         {
             StsException.outputWarningException(this, "getPropertyFields", e);
             return null;
         }
     }

     public boolean getInside()
     {
         return inside;
     }

     public void setInside(boolean val)
     {
         inside = val;
     }

     public boolean checkAlarm(double[] vals)
     {
         if (!enable) return false;

         boolean insideRange = false;
         float halfRange = value + (value * (threshold / 100.0f));
         if ((vals[3] < (value - halfRange)) || (vals[3] > (value + halfRange)))
             insideRange = true;

         if ((insideRange && inside) || (!insideRange && !inside))
         {
             soundAlarm("Value range alarm triggered with value of " + vals[3]);
             return true;
         }
         return false;
     }

     public float getThreshold() { return threshold; }

     public float getValue() { return value; }

     public void setThreshold(float threshold)
     {
         this.threshold = threshold;
     }

     public void setValue(float value)
     {
         this.value = value;
     }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }
 }