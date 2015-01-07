package com.Sts.Types.PreStack;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 11:27:05 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsSeismicReceivers extends StsSeismicSourcesReceivers
{
    private StsSeismicReceivers(StsPreStackLineSet lineSet, StsModel model)
    {
        super(lineSet, model);
        normalColor = StsColor.GREEN;
	    xString = StsSEGYFormat.REC_X;
		yString = StsSEGYFormat.REC_Y;
    }

   public static StsSeismicReceivers constructor(StsPreStackLineSet lineSet, StsModel model)
   {
        if(lineSet.getAttributeIndex("RELEV") == -1) return null;
        if(lineSet.getAttributeIndex(StsSEGYFormat.REC_X) == -1) return null;
        if(lineSet.getAttributeIndex(StsSEGYFormat.REC_Y) == -1) return null;
        StsSeismicReceivers receivers = new StsSeismicReceivers(lineSet, model);
        return receivers;
   }

    public boolean computeBytes()
    {
        StsPreStackLine currentLine = lineSet.currentLine;
        if( currentLine == null) return false;

        if ( (currentLine.getAttributeIndex("RELEV") == -1) || (currentLine.getAttributeIndex(StsSEGYFormat.REC_X) == -1) || (currentLine.getAttributeIndex(StsSEGYFormat.REC_Y) == -1))
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "No Receiver Elevations and/or XY's were selected for storage during processing");
            return false;
        }
        return computeBytes("RELEV", StsSEGYFormat.REC_X, StsSEGYFormat.REC_Y);
    }

    public String getName() { return StsPreStackLineSet.ATTRIBUTE_RELEV_STRING; }

    public int getIndex() { return StsPreStackLineSet.ATTRIBUTE_RELEV; }
}
