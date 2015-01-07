package com.Sts.IO;

import java.io.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsFilenameEndingFilter extends StsFilenameFilter
{
	String[] filter = null;
	int length;

    public StsFilenameEndingFilter() {}

	public StsFilenameEndingFilter(String[] filter)
	{
        super();
		this.filter = new String[filter.length];
		for (int i = 0; i < filter.length; i++)
		{
			this.filter[i] = filter[i].toLowerCase();
		}
	}

	public boolean accept(File dir, String name)
	{
		for (int i = 0; i < filter.length; i++)
		{
			if (name.toLowerCase().endsWith(filter[i]))
				return true;
		}
		return false;
	}

	public String getFilenameStemName(String filename)
	{
		return getFilenameName(filename);
	}

	public String getFilenameName(String filename)
	{
		int filterStart = 1;
		for (int i = 0; i < filter.length; i++)
		{
			if (filename.toLowerCase().indexOf(filter[i]) > 0)
			{
				filterStart = filename.toLowerCase().indexOf(filter[i]);
				break;
			}
		}
		return filename.substring(0, filterStart - 1);
	}

    public String getFilenameEnding(String filename)
    {
        int dotIndex = filename.lastIndexOf(".");
        return filename.substring(dotIndex + 1);
    }
}
