
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Utilities.*;

/** A highlightedList is a temporary list of StsObjects which can be highlighted
  * When an item is added to the list, its highlight property is toggled.
  * The entire list can then be cleared: all objects are toggled off and the
  * list items are deleted.
  */

public class StsHighlightedList extends StsObjectList
{

	public StsHighlightedList(int initSize, int incSize)
	{
    	super(initSize, incSize);
	}

    public void clear()
    {
    	for(int n = 0; n < getSize(); n++)
        {
        	try
            {
        		StsObject object = (StsObject)getElement(n);
            	object.setHighlight(false);
            }
            catch(Exception e)
            {
            	StsException.outputException(e, StsException.WARNING);
            }
        }

        deleteAll();
    }

    public int add(StsObject object)
    {
    	try
        {
        	int index = super.add(object);
    		object.setHighlight(true);
        	return index;
        }
        catch(Exception e)
        {
        	StsException.outputException(e, StsException.WARNING);
            return -1;
        }
    }

    /** Called when we want to clear the list and add the first object */
    public void clearAdd(StsObject object)
	{
    	clear();
        add(object);
    }

    public void display(StsGLPanel glPanel)
    {
        for(int n = 0; n < getSize(); n++)
        {
        	try
            {
            	StsObject object = (StsObject)getElement(n);
                if(object != null) object.display(glPanel);
            }
            catch(Exception e)
            {
            	StsException.outputException(e, StsException.WARNING);
            }
        }
    }
}
