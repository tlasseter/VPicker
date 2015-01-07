
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;

public class StsNote extends StsObject
{
	public String note = null;
	public long timeStamp = 0;
    public String author = null;

	public StsNote()
	{
	}

    public StsNote(String note)
    {
        super(false);
        this.note = note;
        this.author = new String("unknown");
        this.timeStamp = System.currentTimeMillis();
		addToModel();
     }

    public boolean initialize(StsModel model)
    {
        return true;
    }
}




