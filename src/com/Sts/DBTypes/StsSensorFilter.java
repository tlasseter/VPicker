//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.StsTreeObjectI;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.StsObjectPanel;
import com.Sts.UI.Sounds.StsSound;
import com.Sts.UI.StsMessage;
import com.Sts.Utilities.StsException;

public class StsSensorFilter extends StsMainObject implements StsTreeObjectI
{
    // Event symbols
    transient public static final byte QUALIFY = 0;
    transient public static final byte PARTITION = 1;
    transient public static final byte PROXIMITY = 2;
    transient static public final String[] FILTER_TYPE_STRINGS = new String[] { "Qualify", "Partition", "Proximity"};
    transient static public final byte[] FILTER_TYPES = new byte[] { QUALIFY, PARTITION, PROXIMITY};

    protected boolean enable = true;
    protected boolean oneOfMany = false;    // If false, it cannot be applied with other filters.

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;
    static protected StsObjectPanel objectPanel = null;

    /** default constructor */
    public StsSensorFilter()
    {
    	super();
    }

    public StsSensorFilter(boolean persistent)
    {
        super(persistent);
    }

    public StsSensorFilter(byte type, boolean oneOfMany)
    {
        this.type = type;
        this.oneOfMany = oneOfMany;
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
    /**
     * Get an array of the user available field beans
     * @return
     */
    public StsFieldBean[] getDisplayFields()
    {
       try
       {
           if (displayFields == null)
           {
               displayFields = new StsFieldBean[]
               {
                    new StsBooleanFieldBean(StsSensorFilter.class, "isVisible", "Enable"),
               };
           }
            return displayFields;
        }
        catch (Exception e)
        {
            StsException.outputException("StsSensorFilter.getDisplayFields() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
            {
                    new StsStringFieldBean(StsSensorFilter.class, "name", true, "Name:"),
                    new StsStringFieldBean(StsSensorFilter.class, "typeString", false, "Type:"),
                    new StsBooleanFieldBean(StsSensorFilter.class, "oneOfMany", "Multiple Allowed:", false, true),
            };
        }
        return propertyFields;
    }

    public String getTypeString()
    {
        return FILTER_TYPE_STRINGS[type];
    }
    public void setTypeString(String typeS)
    {
        for(int i=0; i<FILTER_TYPE_STRINGS.length; i++)
        {
            if(FILTER_TYPE_STRINGS[i].equalsIgnoreCase(typeS))
            {
                setType((byte)i);
                return;
            }
        }
    }
    public boolean getOneOfMany() { return oneOfMany;}
    public void setOneOfMany(boolean val) { oneOfMany = val; }

    public boolean delete()
    {
        boolean val = super.delete();

        // Must update the static beans list of filters
        StsObject[] sensors = ((StsSensorClass)currentModel.getCreateStsClass("com.Sts.DBTypes.StsSensor")).getSensors();
        for(int i=0; i<sensors.length; i++)
        {
            ((StsSensor)sensors[i]).updateFilters();
        }
        return val;
    }

    public boolean filter(StsSensor sensor) { return false; }  // Override in each filter

    public Object[] getChildren() { return new Object[0]; }

    public boolean anyDependencies() { return false; }

    public void updateProperties() {   }

    public void treeObjectSelected() {  currentModel.getCreateStsClass(this).selected(this); }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }
}