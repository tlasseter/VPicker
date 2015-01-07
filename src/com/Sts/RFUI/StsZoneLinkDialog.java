
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.RFUI;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

public class StsZoneLinkDialog extends StsSelectStsObjects
{
	static public StsZoneLinkDialog constructor(StsModel model, String title)
	{
		try
		{
			return (StsZoneLinkDialog)StsSelectStsObjects.constructor(model, StsWellZoneSet.class, title,
				"Select a well zone (links shown):", true, false);
		}
		catch(Exception e)
		{
		    StsException.outputException("StsZoneLinkDialog constructor() failed.",
				e, StsException.WARNING);
			return null;
		}
	}

    // overrides StsSelectStsObjects method
    protected void setSelector(StsObject[] objects)
    {
        String[] names;
        StsColor[] colors;
        if (objects!=null)
        {
            names = new String[objects.length];
            colors = new StsColor[objects.length];
            for (int i=0; i<objects.length; i++)
            {
                StsWellZoneSet wellZone = (StsWellZoneSet)objects[i];
                names[i] = wellZone.getName();
                StsObject[] zoneLinks = { wellZone.getParentZone() };
                if (zoneLinks[0]==null) names[i] += "   ( )";
                else names[i] += "   (" + zoneLinks[0].getName() + ")";
                colors[i] = wellZone.getStsColor();
            }
        }
        else  // empty list
        {
            colors = StsListStsObjects.NULL_COLORS;
            names = StsListStsObjects.NULL_NAMES;
        }

        setSelector(colors, names);
    }
}
