package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import java.util.*;

public class StsWellMarkerClass extends StsClass implements StsSerializable
{
    public String[] orderedMarkerNames = null; // bottom to top ordering
    
    public StsWellMarkerClass()
    {
        userName = "Geologic Markers";                                                                      

//        initColors(StsWellMarker.displayFields);
    }

    public StsWellMarker[] getWellMarkerSet(StsMarker marker)
    {
        if (marker == null) return null;
        int nMarkers = getSize();
        ArrayList tempSet = new ArrayList();
        for (int i = 0; i < nMarkers; i++)
        {
            StsWellMarker wellMarker = (StsWellMarker)getElement(i);
            if (wellMarker.getMarker() == marker) tempSet.add(tempSet);
        }
        int nInSet = tempSet.size();
        if (nInSet == 0) return null;
        StsWellMarker[] wellMarkerSet = new StsWellMarker[nInSet];
        for (int i = 0; i < nInSet; i++)
            wellMarkerSet[i] = (StsWellMarker) tempSet.get(i);
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
			StsWellMarker marker = (StsWellMarker)iter.next();
			marker.display(glPanel3d);
		}
	}

}
