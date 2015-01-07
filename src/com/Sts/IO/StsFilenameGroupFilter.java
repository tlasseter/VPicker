
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.IO;

import com.Sts.Actions.Import.*;

import java.io.*;

public class StsFilenameGroupFilter implements FilenameFilter
{
    String group;
    String format;
    String name;

    public StsFilenameGroupFilter(String group, String format, String name)
    {
        this.group = group;
        this.format = format;
        this.name = name;
    }

    public boolean accept(File file)
    {
        if (file.isDirectory()) return true;
        if (!file.isFile()) return false;
        return accept(null, file.getName());
    }


    public boolean accept(String filename)
    {
        return accept(null, filename);
    }

    public boolean accept(File dir, String filename)
    {
        StsKeywordIO.parseAsciiFilename(filename);
        if(group != null && !group.equals(StsKeywordIO.group)) return false;
        if(format != null && !format.equals(StsKeywordIO.format)) return false;
        if(name != null && !name.equals(StsKeywordIO.name)) return false;
        return true;
   }
}
