
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.UI;

import javax.swing.filechooser.FileFilter;
import java.io.*;

public class StsFileFilter extends FileFilter
{
    String filter;
    boolean isPrefix;

    public StsFileFilter(String filter, boolean isPrefix)
    {
        this.filter = filter;
        this.isPrefix = isPrefix;
    }

    public boolean accept(File file)
    {
        if (file.isDirectory()) return true;
        if (!file.isFile()) return false;
        if( filter == null ) return true;

        if (isPrefix) return file.getName().startsWith(filter);
        else return file.getName().endsWith(filter);
    }

    public boolean accept(String filename)
    {
        if( filter == null ) return true;
        if (isPrefix) return filename.startsWith(filter);
        else return filename.endsWith(filter);
    }

	public String getDescription()
	{
		return filter;
	}
}

