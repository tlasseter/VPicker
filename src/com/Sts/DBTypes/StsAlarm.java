//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.UI.Beans.*;
import com.Sts.UI.Sounds.*;
import com.Sts.UI.StsMessage;
import com.Sts.UI.StsSelectable;
import com.Sts.UI.ObjectPanel.StsObjectPanel;
import com.Sts.Utilities.*;
import com.Sts.MVC.StsMessageFiles;
import com.Sts.MVC.StsModel;
import com.Sts.Interfaces.StsTreeObjectI;

abstract public class StsAlarm extends StsMainObject
{
    // Event symbols
    transient public static final byte WELL = 0;
    transient public static final byte SURFACE = 1;
    transient public static final byte VALUE = 2;
    transient static public final String[] ALARM_TYPE_STRINGS = new String[] { "Well", "Surface", "Value"};
    transient static public final byte[] ALARM_TYPES = new byte[] { WELL, SURFACE, VALUE };
    // Event symbols
    transient public static final byte AUDIO = 0;
    transient public static final byte VISUAL = 1;
    transient public static final byte MESSAGE = 2;
    transient public static final byte AUDIO_VISUAL = 3;
    transient public static final byte ALL = 4;
    transient static public final String[] QUEUE_TYPE_STRINGS = new String[] { "Audio", "Visual", "Message", "Audio & Visual", "All"};
    transient static public final byte[] QUEUE_TYPES = new byte[] { AUDIO, VISUAL, MESSAGE, AUDIO_VISUAL, ALL };

    String sound = null;
    byte queue = AUDIO;
    protected boolean enable = true;

    abstract public StsFieldBean[] getDisplayFields();
    abstract public StsObjectPanel getObjectPanel(); 
    /** default constructor */
    public StsAlarm()
    {
    	super();
    }

    public StsAlarm(String sound)
    {
        this.sound = sound;
    }

    public boolean initialize(StsModel model) { return true; }

    public boolean getEnable()
    {
    	return enable;
    }    
    public void setEnable(boolean val)
    {
    	enable = val;
    }
    public String getSound()
    {
    	return sound;
    }
    public void setSound(String sound)
    {
    	this.sound = sound;
        StsSound.play(sound);
    }
    public void soundAlarm()
    {
        soundAlarm("An alarm has been triggered from unknown source.");
    }
    public void soundAlarm(String message)
    {
        // Sound Alarm
        if((sound != null) && ((queue == AUDIO) || (queue == AUDIO_VISUAL) || (queue == ALL)))
            StsSound.play(sound);

        // Visual Alarm
        if((queue == VISUAL) || (queue == AUDIO_VISUAL) || (queue == ALL))
            ;

        // Message Alarm
        if((queue == MESSAGE) || (queue == ALL))
            new StsMessage(currentModel.win3d, StsMessage.WARNING,message);             
    }
    public byte getQueue() { return queue; }
    public String getQueueString()
    {
    	return QUEUE_TYPE_STRINGS[queue];
    }
    public void setQueueString(String queueString)
    {
        getQueueType(queueString);
    }
    static public byte getAlarmType(String alarm)
    {
        for(int i=0; i<=ALARM_TYPE_STRINGS.length; i++)
        {
            if(alarm.equalsIgnoreCase(ALARM_TYPE_STRINGS[i]))
            {
                return (byte)i;
            }
        }
    	return WELL;
    }
    static public byte getQueueType(String queue)
    {
        for(int i=0; i<=QUEUE_TYPE_STRINGS.length; i++)
        {
            if(queue.equalsIgnoreCase(QUEUE_TYPE_STRINGS[i]))
            {
                return (byte)i;
            }
        }
    	return AUDIO;
    }
    public boolean checkAlarm(double[] xyz)
    {
        System.err.println("checkAlarm(float[] xyz) not implemented in class " + getClassname());
        return false;
    }

    public StsFieldBean[] getPropertyFields() { return null; }

    public Object[] getChildren() { return new Object[0]; }

    public boolean anyDependencies() { return false; }

    public void updateProperties() {   }

    public void treeObjectSelected() {  currentModel.getCreateStsClass(this).selected(this); }
}