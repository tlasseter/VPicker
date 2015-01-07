package com.Sts.DBTypes;

import com.Sts.UI.Beans.*;
import com.Sts.DB.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Jul 15, 2010
 * Time: 7:14:45 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsBuiltModelClass extends StsClass implements StsSerializable
{
    protected byte displayMode = 0;

    static public final byte displayZones = 0;
    static public final byte displayBlocks = 1;
    static public final byte displayLayers = 2;
    static public final byte displayProperties = 3;
    static public final String displayZonesString = "Zones";
    static public final String displayBlocksString = "Blocks";
    static public final String displayLayersString = "Layers";
    static public final String displayPropertiesString = "Properties";
    static public final String[] displayModeStrings = new String[] {displayZonesString, displayBlocksString, displayLayersString, displayPropertiesString};

    public StsBuiltModelClass()
    {
        userName = "GeoCellular Models";
    }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
        {
            new StsButtonListFieldBean(this, "displayModeString", "Select Display Mode", displayModeStrings, false)
        };
    }

    public void setDisplayModeString(String modeString)
    {
        for(byte n = 0; n < displayModeStrings.length; n++)
        {
            if(modeString == displayModeStrings[n])
            {
                displayMode = n;
                currentModel.viewObjectRepaint(this, this);
                break;
            }
        }
    }

    static public String staticGetDisplayModeString()
    {
        StsBuiltModelClass builtModelClass = (StsBuiltModelClass)currentModel.getStsClass(StsBuiltModel.class);
        if(builtModelClass == null) return StsBuiltModelClass.displayZonesString;
        return builtModelClass.getDisplayModeString();
    }

    public String getDisplayModeString()
    {
        return displayModeStrings[displayMode];
    }

    static public void staticSetDisplayMode(byte mode)
    {
        StsBuiltModelClass builtModelClass = (StsBuiltModelClass)currentModel.getStsClass(StsBuiltModel.class);
        builtModelClass.displayMode = mode;    
    }
}
