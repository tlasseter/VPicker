
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

public class StsWellZoneLinkDialog extends StsSelectStsObjects
{
	static public StsWellZoneLinkDialog constructor(StsModel model)
	{
		try
		{
			return (StsWellZoneLinkDialog)StsSelectStsObjects.constructor(model, StsZone.class, model.getName(),
				"Select a zone (links shown):", true, false);
		}
		catch(Exception e)
		{
		    StsException.outputException("StsWellZoneLinkDialog.constructor() failed.",
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
                StsZone zone = (StsZone)objects[i];
                names[i] = zone.getName();
                StsObject[] wellZoneLinks = { zone.getWellZoneSet() };
                if (wellZoneLinks[0]==null) names[i] += "   ( )";
                else names[i] += "   (" + wellZoneLinks[0].getName() + ")";
                colors[i] = zone.getStsColor();
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
