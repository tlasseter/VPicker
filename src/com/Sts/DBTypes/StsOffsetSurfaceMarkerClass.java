package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.StsSerializable;
import com.Sts.MVC.View3d.StsGLPanel3d;
import com.Sts.Utilities.StsException;

import java.util.ArrayList;
import java.util.Iterator;

public class StsOffsetSurfaceMarkerClass extends StsWellMarkerClass implements StsSerializable
{
    public StsOffsetSurfaceMarkerClass()
    {
        userName = "Offset Surface Markers";
//        initColors(StsWellMarker.displayFields);
    }

    public StsOffsetSurfaceMarker[] getOffsetSurfaceMarkerSet(StsMarker marker)
    {
        userName = "Offset Surface Markers";

        if (marker == null) return null;
        int nMarkers = getSize();
        ArrayList tempSet = new ArrayList();
        for (int i = 0; i < nMarkers; i++)
        {
            StsOffsetSurfaceMarker wellMarker = (StsOffsetSurfaceMarker)getElement(i);
            if (wellMarker.getMarker() == marker) tempSet.add(tempSet);
        }
        int nInSet = tempSet.size();
        if (nInSet == 0) return null;
        StsOffsetSurfaceMarker[] wellMarkerSet = new StsOffsetSurfaceMarker[nInSet];
        for (int i = 0; i < nInSet; i++)
            wellMarkerSet[i] = (StsOffsetSurfaceMarker) tempSet.get(i);
        return wellMarkerSet;
    }

    /** set/get the orderedMarkerNames list */
    public void setOrderedMarkerNames(String[] names)
    {
        orderedMarkerNames = names;
    }

    public String[] getOrderedMarkerNames()
    {
        return getOrderedMarkerNames(true);
    }

    public String[] getOrderedMarkerNames(boolean ascending)
    {
        try
        {
            orderedMarkerNames = StsOrderedList.getOrderedNames(this, orderedMarkerNames, ascending);
            return orderedMarkerNames;
        }
        catch (Exception e)
        {
            StsException.outputException(e, StsException.WARNING);
            return null;
        }
    }

	public void displayClass(StsGLPanel3d glPanel3d)
	{
		Iterator iter = getObjectIterator();
		while(iter.hasNext())
		{
			StsOffsetSurfaceMarker marker = (StsOffsetSurfaceMarker)iter.next();
			marker.display(glPanel3d);
		}
	}

}