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

public class StsSurfaceAlarm extends StsAlarm implements StsTreeObjectI
{
    // Need reference direction
    static public final byte ABOVE = 0;
    static public final byte BELOW = 1;
    static public final String[] REF_DIRECTION = new String[]{"Above", "Below"};

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;
    static protected StsObjectPanel objectPanel = null;

    // Event symbols
    StsSurface surface = null;
    float distanceToSurface = 0.0f;
    byte direction = BELOW;     // Direction from surface

    /** default constructor */
    public StsSurfaceAlarm()
    {
        super();
    }

    public StsSurfaceAlarm(String sound, StsSurface surface, float distance, byte direction)
    {
        super(sound);
        this.surface = surface;
        this.distanceToSurface = Math.abs(distance);
        this.direction = direction;
        this.name = surface.getName() + "_Alarm";
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
            StsException.outputException("StsSurfaceAlarm.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsFieldBean[] getPropertyFields()
    {
        if (propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
            {
                new StsComboBoxFieldBean(StsSurfaceAlarm.class, "referenceDirection", "Direction:", REF_DIRECTION),
                new StsFloatFieldBean(StsSurfaceAlarm.class, "distanceToSurface", true, "Distance to Surface:")
            };
        }
        return propertyFields;
    }

    public String getReferenceDirection() { return REF_DIRECTION[direction]; }

    public void setReferenceDirection(String dir)
    {
        if (dir.equalsIgnoreCase(REF_DIRECTION[BELOW]))
            direction = BELOW;
        else
            direction = ABOVE;
    }

    public boolean checkAlarm(double[] xyz)
    {
        if (!enable) return false;

        // Check Alarm
        StsPoint surfacePoint = surface.getInterpolatedPoint((float)xyz[0], (float)xyz[1]);
        Float distance = (float)xyz[2] - surfacePoint.getZ();
        if ((direction == ABOVE) && (xyz[2] < surfacePoint.getZ()))
        {
            if (Math.abs(distance) < distanceToSurface)
            {
                soundAlarm("Surface proximity alarm triggered for " + surface.getName());
                return true;
            }
        }
        else if ((direction == BELOW) && (xyz[2] > surfacePoint.getZ()))
        {
            if (Math.abs(distance) < distanceToSurface)
            {
                soundAlarm("Surface proximity alarm triggered for " + surface.getName());
                return true;
            }
        }
        return false;
    }

    public float getDistanceToSurface() { return distanceToSurface; }

    public void setDistanceToSurface(float distance)
    {
        this.distanceToSurface = Math.abs(distance);
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