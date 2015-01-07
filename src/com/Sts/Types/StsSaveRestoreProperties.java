
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;


/** This class defines String keys paired with integer values *
  * These need to be saved in the database as they define the state of
  * any flags of interest.
  *
  * Typically call save() before an operation (e.g., turn off sections),
  * change view properties for operation, and then call restore() after.
  */

public class StsSaveRestoreProperties extends StsProperties
{
    transient int[] saveValues = null; /* A temporary list of saved int values */
    protected StsIntList currentValues = null; /* A vector of ints containing current values */

	public StsSaveRestoreProperties()
	{
        super(10);
 	}

	public StsSaveRestoreProperties(int initSize)
	{
    	super(initSize);
 	}

    /** Save all the key values to a temporary int list */
    public void save()
    {
		if(currentValues == null)
		{
			saveValues = null;
			return;
		}

    	int size = currentValues.getSize();
    	saveValues = new int[size];
        for(int n = 0; n < size; n++)
        	saveValues[n] = currentValues.list[n];
    }

    /** Restore the value list from the temp list */
    public void restore()
    {
		if(saveValues == null) return;
    	int size = saveValues.length;

        for(int n = 0; n < size; n++)
        	currentValues.list[n] = saveValues[n];
    }

   	public static void main(String[] args)
    {
        StsSaveRestoreProperties properties = new StsSaveRestoreProperties(10);

        properties.set("Sections", 1);
        properties.set("Edges", 3);
        properties.set("Patches", 5);

        properties.print();

        properties.save();

        properties.set("Edges", -3);
        properties.print();

        properties.restore();
        properties.print();
    }
}

