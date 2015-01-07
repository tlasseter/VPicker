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
import com.Sts.Types.StsPoint;
import com.Sts.Interfaces.*;

public class StsWellAlarm extends StsAlarm implements StsTreeObjectI
 {
     // Event symbols
     StsWell well = null;
     float distanceToWell = 0.0f;
     boolean inside = false;  // default to trigger alarm outside range.

     static public StsFieldBean[] displayFields = null;
     static public StsFieldBean[] propertyFields = null;
     static protected StsObjectPanel objectPanel = null;

     /** default constructor */
     public StsWellAlarm()
     {
         super();
     }

     public StsWellAlarm(String sound, StsWell well, float offset, boolean inside)
     {
         super(sound);
         this.well = well;
         this.inside = inside;
         this.distanceToWell = offset;
         this.name = well.getName() + "_WellAlarm";
     }
     
     public StsFieldBean[] getDisplayFields()
     {
         try
         {
             if (displayFields == null)
             {
                 displayFields = new StsFieldBean[]
                 {
                     new StsBooleanFieldBean(this, "enable", "Enable Alarm:"),
                     new StsComboBoxFieldBean(this, "Sound", "Sound:", StsSound.sounds),
                     new StsComboBoxFieldBean(this, "QueueString", "Type:", QUEUE_TYPE_STRINGS)
                 };
             }
             return displayFields;
         }
         catch (Exception e)
         {
             StsException.outputException("StsWellAlarm.getDisplayFields() failed.", e, StsException.WARNING);
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
                     new StsBooleanFieldBean(StsWellAlarm.class, "Inside", "Inside Range", true),
                     new StsFloatFieldBean(StsWellAlarm.class, "distanceToWell", true, "Distance to Well:")
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

     public float getDistanceToWell() { return distanceToWell; }

     public void setDistanceToWell(float distance)
     {
         this.distanceToWell = distance;
     }

     public boolean checkAlarm(double[] xyz)
     {
         if (!enable) return false;

 // Check Alarm
         float[] xy = currentModel.getProject().getRelativeXY(xyz[0], xyz[1]);
         StsPoint point = new StsPoint(xy[0],xy[1],xyz[2],xyz[3]);
         StsPoint wellPoint = StsMath.getNearestPointOnLine(point, well.getRotatedPoints(), 3);
         Float distance = StsMath.distance(point.getXYZ(), wellPoint.getXYZ(), 3);
         if (((distance > distanceToWell) && (!inside)) || ((distance < distanceToWell) && (inside)))
         {
             soundAlarm("Well proximity alarm triggered for " + well.getName());
             return true;
         }
         return true;
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