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

public class StsFMIMarkerClass extends StsWellMarkerClass implements StsSerializable
{
    public StsFMIMarkerClass()
    {
        userName = "Formation Image Scan Markers";                                                              
//        initColors(StsWellMarker.displayFields);
    }

    public StsFMIMarker[] getFMIMarkerSet(StsMarker marker)
    {
        if (marker == null) 
        	return null;
        int nMarkers = getSize();
        ArrayList tempSet = new ArrayList();
        for (int i = 0; i < nMarkers; i++)
        {
            StsFMIMarker wellMarker = (StsFMIMarker)getElement(i);
            if (wellMarker.getMarker() == marker) tempSet.add(tempSet);
        }
        int nInSet = tempSet.size();
        if (nInSet == 0) return null;
        StsFMIMarker[] wellMarkerSet = new StsFMIMarker[nInSet];
        for (int i = 0; i < nInSet; i++)
            wellMarkerSet[i] = (StsFMIMarker) tempSet.get(i);
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
			StsFMIMarker marker = (StsFMIMarker)iter.next();
			marker.display(glPanel3d);
		}
	}

}
