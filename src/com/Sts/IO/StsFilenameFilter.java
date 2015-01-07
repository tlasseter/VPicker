package com.Sts.IO;

import com.Sts.Actions.Import.*;

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
public class StsFilenameFilter implements StsFilenameFilterFace
{
	public String format = null;
	public String group = null;

    public StsFilenameFilter() {}
	public StsFilenameFilter(String grp, String fmt)
    {
        this.format = fmt;
        this.group = grp;
    }

	public boolean accept(File dir, String filename)
	{
	  return StsKeywordIO.fileMatchesGroupAndFormat(filename, group, format);
	}

    public String getFilenameName(String filename)
    {
        return StsKeywordIO.getFileEndString(filename);
    }
    public String getFilenameStemName(String filename)
    {
        return StsKeywordIO.getFileStemName(filename);
    }
    public String getFilenameEnding(String filename)
    {
        return StsKeywordIO.getFileEndString(filename);
    }
}
