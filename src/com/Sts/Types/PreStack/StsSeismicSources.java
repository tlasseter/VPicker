package com.Sts.Types.PreStack;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 11:24:34 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSeismicSources extends StsSeismicSourcesReceivers
{
    private StsSeismicSources(StsPreStackLineSet lineSet, StsModel model)
    {
        super(lineSet, model);
        normalColor = StsColor.RED;
	    xString = StsSEGYFormat.SHT_X;
		yString = StsSEGYFormat.SHT_Y;
    }

    public static StsSeismicSources constructor(StsPreStackLineSet lineSet, StsModel model)
    {
        if(lineSet.getAttributeIndex("SELEV") == -1) return null;
        if(lineSet.getAttributeIndex(StsSEGYFormat.SHT_X) == -1) return null;
        if(lineSet.getAttributeIndex(StsSEGYFormat.SHT_Y) == -1) return null;
        StsSeismicSources sources = new StsSeismicSources(lineSet, model);
        return sources;
    }

    public String getName() { return StsPreStackLineSet.ATTRIBUTE_SELEV_STRING; }

    public int getIndex() { return StsPreStackLineSet.ATTRIBUTE_SELEV; }

    public boolean computeBytes()
    {
        StsPreStackLine currentLine = lineSet.currentLine;
        if( currentLine == null) return false ;
        if ( (currentLine.getAttributeIndex("SELEV") == -1) || (currentLine.getAttributeIndex("SHT-X") == -1) || (currentLine.getAttributeIndex("SHT-Y") == -1))
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "No Source Elevations and/or XY's were selected for storage during processing");
            return false;
        }
        return computeBytes("SELEV", StsSEGYFormat.SHT_X, StsSEGYFormat.SHT_Y);
    }
}
