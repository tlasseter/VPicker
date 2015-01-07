package com.Sts.Actions.Import;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.MVC.*;
import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

public class StsKeywordIO
{
    static StsProject project;

    static public String filename = null;
    static public String group = null; // well, well-log, etc
    static public byte GROUP = 0;
    static public String format = null; // txt (ascii), bin (binary)
    static public byte FORMAT = 1;
    static public String name = null; // name of object (example: wellname)
    static public byte NAME = 2;
    static public String subname = null; // name of subobject (example: curvename)
    static public byte SUBNAME = 3;
    static public int version = 0;  // version numbers if there identical filenames
	static public long fileCreationDate = 0; // date this file was created
    static public int[] tokenOrder = { GROUP, FORMAT, NAME, SUBNAME };
    static public String delimiter = ".";
    static public int numberFields = 3;
    static public boolean hasVersion = false;

    static public void initialize(StsModel model)
	{
       delimiter = ".";
       numberFields = 3;
       tokenOrder = new int[] { GROUP, FORMAT, NAME, SUBNAME };
    }

    static public boolean setParseFilename(String _filename, String _group, String _format)
    {
        filename = _filename;
        if (!parseAsciiFilename(filename))
            return false;

		if(!group.equals(_group) || !format.equals(_format))
		{
			StsMessageFiles.logMessage(" File: " + filename + " is not of type: " + _group + "." + _format);
		    return false;
        }
		return true;
    }

    static public void setTokenOrder(int[] order) { tokenOrder = order; }
    static public int getTokenOrder(int type) { return tokenOrder[type]; }

    static public void setDelimiter(String delimit) { delimiter = delimit; }
    static public String getDelimiter() { return delimiter; }

    static public void setNumberFields(int numFields) { numberFields = numFields; }
    static public int getNumberFields() { return numberFields; }

    static public boolean isFileAscii(String filename)
    {
        if(!parseAsciiFilename(filename)) return false;
        return format.equals("txt");
    }
    
    /** these are multi-column files with name: prefix.format.name or prefix.format.name.version where version is an integer */
    static public boolean parseAsciiFilename(String filename)
    {
		File file = new File(filename);
		fileCreationDate = file.lastModified();
        int tokenCount = 0;

        String[] fields = new String[numberFields];
        for(int i=0; i<numberFields; i++) fields[i] = "";
        group = ""; format = ""; name = ""; version = 0;
        StringTokenizer sTokens = new StringTokenizer(filename, delimiter);
        while(sTokens.hasMoreTokens())
        {
            fields[tokenCount] = sTokens.nextToken();
            tokenCount++;
            if(tokenCount == numberFields)
                break;
        }
        if(tokenCount != numberFields)
            return false;
        if(sTokens.hasMoreTokens())
        {
            try
            {
                version = Integer.parseInt(sTokens.nextToken());
                hasVersion = true;
            }
            catch(Exception e) { }
        }

        if(tokenCount > tokenOrder[GROUP]) group = fields[tokenOrder[GROUP]];
        if(tokenCount > tokenOrder[FORMAT]) format = fields[tokenOrder[FORMAT]];
        if(tokenCount > tokenOrder[NAME]) name = fields[tokenOrder[NAME]];
/*        if(tokenCount > tokenOrder[NAME])
        {
            if(hasVersion)
                name = fields[tokenOrder[NAME]];
            else
            {
                int length = new String(group + delimiter + format + delimiter).length();
                name = filename.substring(length);
                if(name.indexOf(delimiter) > 0)
                    name = name.substring(0,name.indexOf(delimiter));
            }
        }*/
        return true;
    }

    // ascii multicolumn files are written out as single-column binary files: prefix.format.name.subname.version
    static public boolean parseBinaryFilename(String filename)
    {
        group = ""; format = ""; name = ""; subname = ""; version = -1;

        StringTokenizer sTokens = new StringTokenizer(filename, ".");
        if(!sTokens.hasMoreTokens()) return false;
        group = sTokens.nextToken();
        if(!sTokens.hasMoreTokens()) return false;
        format = sTokens.nextToken();
        if(!sTokens.hasMoreTokens()) return false;
        name = sTokens.nextToken();
        if(!sTokens.hasMoreTokens()) return false;
        subname = sTokens.nextToken();

        if(sTokens.hasMoreTokens())
        {
            try {  version = Integer.parseInt(sTokens.nextToken()); }
            catch(Exception e) { version = 0; }
        }
        else
            version = 0;

        return true;
    }

    static public boolean fileMatchesGroupAndFormat(String filename, String _group, String _format)
    {
        if(!parseAsciiFilename(filename)) return false;
        return group.equals(_group) && format.equals(_format);
    }

    static public String getFileEndString(String filename)
    {
        if(!parseAsciiFilename(filename)) return null;
        int prefixLength = new String(group + delimiter + format + delimiter).length();
        return filename.substring(prefixLength);
    }

    static public String getFileStemName(String filename)
    {
        if(!parseAsciiFilename(filename)) return null;
        return name;
    }

    static protected boolean parsesOK(String _group, String _format)
    {
        return group.equals(_group) && format.equals(_format);
    }

	static public boolean lineHasKeyword(String line, String keyword)
	{
		return line.indexOf(keyword) >= 0;
	}

	static protected boolean findKeyword(BufferedReader bufRdr, String keyword)
	{
		String line;
		try
		{
			while( (line = bufRdr.readLine().trim()) != null)
				if(line.indexOf(keyword) >= 0) return true;
			return false;
		}
		catch(Exception e)
		{
			StsException.outputException("StsKeywordIO.findKeyword() failed.",
				e, StsException.WARNING);
			return false;
		}
	}

    static protected void closeBufRdr(BufferedReader bufRdr)
    {
        try
        {
            if(bufRdr == null) return;
            bufRdr.close();
        }
        catch(Exception e) { }
    }

    static public String constructFilename(String group, String format, String name, String subname, int version)
    {
        return group + "." + format + "." + name + "." + subname + "." + version;
    }
}
