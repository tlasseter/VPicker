
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Types.*;
import com.Sts.Utilities.*;

public class StsOrderedList
{
    /** return list of names in a user-selected ordering */
    static public String[] getOrderedNames(StsClass instList,
            String[] prevOrderedNames, boolean ascending) throws StsException
    {
        if (instList==null) return null;
        int nItems = instList.getSize();
        if (nItems<1) return null;

        StsList list = null;
        if (prevOrderedNames==null) // first time
        {
            list = new StsList(nItems, 1);
            list.add(instList.getFirst());
            for (int i=1; i<nItems; i++)
            {
                addOrderedItem(list, instList.getElement(i), ascending);
            }
        }
        else
        {
            // check for deletions from instance list
            list = instList.getObjectList(prevOrderedNames);
            int nOrderedItems = list.getSize();
            for (int i=nOrderedItems-1; i>=0; i--)
            {
                StsObject obj = (StsObject)list.getElement(i);
                if (instList.getElement(obj.getName())==null) list.delete(obj);
            }

            // check for additions in instance list
            for (int i=0; i<nItems; i++)
            {
                StsObject obj = instList.getElement(i);
                if (!list.contains(obj)) addOrderedItem(list, obj, ascending);
            }
        }
        return listToStrings(list);
    }

    /** return list of names in a user-selected ordering */
    static public String[] getOrderedNames(StsObject[] items, String[] prevOrderedNames,
            boolean ascending)
        throws StsException
    {
        if (items==null || items.length==0) return null;
        int nItems = items.length;

        StsList list = null;
        if (prevOrderedNames==null) // first time
        {
            list = new StsList(nItems, 1);
            list.add(items[0]);
            for (int i=1; i<nItems; i++)
            {
                addOrderedItem(list, items[i], ascending);
            }
        }
        else
        {
            boolean[] itemsFound = new boolean[nItems];
            for (int i=0; i<nItems; i++) itemsFound[i] = false;

            // build list while checking for deletions
            list = new StsList(prevOrderedNames.length);
            for (int i=0; i<prevOrderedNames.length; i++)
            {
                for (int j=0; j<nItems; j++)
                {
                    if (items[j].getName().equals(prevOrderedNames[i]))
                    {
                        list.add(items[j]);
                        itemsFound[j] = true;
                        break;
                    }
                }
            }

            // check for additions in items array
            for (int i=0; i<nItems; i++)
            {
                if (!itemsFound[i]) addOrderedItem(list, items[i], ascending);
            }
        }
        return listToStrings(list);

    }

    static private boolean addOrderedItem(StsList list, StsObject obj,
            boolean ascending) throws StsException
    {
        if (list==null || obj==null) return false;
        float objValue = obj.getOrderingValue();
        if (objValue==StsParameters.nullValue)
        {
            list.add(obj);
            return true;
        }

        int nItems = list.getSize();
        if (ascending)
        {
            for (int i=0; i<nItems; i++)
            {
                StsObject objI = (StsObject)list.getElement(i);
                float objIValue = objI.getOrderingValue();
                if (objValue<objIValue || objIValue==StsParameters.nullValue)
                {
                    list.insertBefore(objI, obj);
                    return true;
                }
            }
            list.add(obj);  // put at the end
        }
        else  // descending
        {
            for (int i=0; i<nItems; i++)
            {
                StsObject objI = (StsObject)list.getElement(i);
                float objIValue = objI.getOrderingValue();
                if (objValue>objIValue || objIValue==StsParameters.nullValue)
                {
                    list.insertBefore(objI, obj);
                    return true;
                }
            }
            list.add(obj);  // put at the end
        }
        return true;
    }

    static public String[] listToStrings(StsList list)
    {
        int nItems = list.getSize();
        String[] names = new String[nItems];
        for (int i=0; i<nItems; i++)
        {
            StsObject obj = (StsObject)list.getElement(i);
            names[i] = obj.getName();
        }
        return names;
    }

}


