package com.Sts.IO;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.MVC.Main;
import com.Sts.UI.*;
import com.Sts.Utilities.*;
import com.Sts.Actions.Import.*;

import java.awt.*;
import java.io.*;
import java.net.*;

/** This class handles data files processed by S2S data loaders: wells, surfaces, seismic, etc.
 *  All data files have a standard name structure:  group.format.name.subname.version where subname and version are optional.
 *  "group" identifies the type of file, e.g. "seismic-grid" is the group name for surface grids in seismic format.
 *  "format" for these files is either "txt" for ASCII files or "bin" for binary files.  Other files may be in Java object format "obj", for example.
 *  "name" is the name of the object being read in and instantiated.
 *  "subname" is used for objects which are part of another object identified by name.
 *  "version" is an integer version number which is incremented as output files with changes for this same object are written so as not to overwrite
 *  the original file.
 *  The fileType is either an ascii file whose type is the same as group, or a binary file which is either an ordinary binary file, or an entry in a jar
 *  file set or a webstart jar file set.
 */
abstract public class StsAbstractFile implements Comparable
{
    public String filename;
    public URL url;
    public String group;
    public String format;
    public String name;
    public String subname;
    public int version;
    /** fileType is set by the wizard and is used by data loaders to specify the type of loader to use. */
    public String fileType;

    static public final String asciiFormat = "txt";
    static public final String binaryFormat = "bin";

    static public final String WEBJAR = "webJar";
    static public final String JAR = "jar";
    
    static public final String[] invalidCharacters = { "\\", "/", "*", "?", "\"", "<", ">", "|", "#" };

    abstract public InputStream getInputStream() throws IOException;
    abstract public OutputStream getOutputStream(boolean append) throws FileNotFoundException;
    abstract public OutputStream getOutputStreamAndPosition(long position) throws FileNotFoundException;
    abstract public InputStream getMonitoredInputStream(Component component);
	abstract public boolean isWritable();

    public StsAbstractFile()
    {
    }

    public String getFilename() { return filename; }
    public URL getURL() { return url; }

    public int compareTo(Object other)
    {
        String otherFilename = ((StsAbstractFile)other).getFilename();
        return filename.compareTo(otherFilename);
    }

    public String getURLDirectory()
    {
        try
        {
            String path = url.toURI().getPath();
            int last = path.lastIndexOf("/");
            if(last == -1) return new String("./");
            return new String(path.toCharArray(), 0, last+1);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getURLDirectory", e);
            return null;
        }
    }

    public String getPathname() { return getURLPathname(); }

    public String getURLPathname()
    {
        try
        {
            if(Main.isJarDB)
                return url.getPath();
            else
                return url.toURI().getPath();
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "getURLDirectory", e);
            return null;
        }
    }

    public long length()
    {
        File file = getFile();
        if(file == null) return 0;
        return file.length();
    }

    public File getFile()
    {
        try
        {
            if(this instanceof StsJarEntryFile)
            {
                return new File(url.toExternalForm());
            }
            else
            {
                URI uri = url.toURI();
                return new File(uri);
            }
        }
        catch(Exception e)
        {
            StsException.systemError(this, "getFile", "Failed for URL: " + url.getPath());
            return null;
        }
    }

    public boolean createNewFile()
    {
        File file = getFile();
        if(file == null) return false;
        try
        {
            file.createNewFile();
            return true;
        }
        catch(Exception e)
        {
            return false;
        }
    }

    public boolean isAFile()
    {
        File file = getFile();
        if(file == null) return false;
        return file.isFile();
    }

    public boolean isAFileAndNotDirectory()
    {
        File file = getFile();
        if(file == null) return false;
        return !file.isDirectory() && file.isFile();
    }
    
	public boolean exists()
	{
         File file = getFile();
         if(file == null) return false;
         return file.exists();
    }

    public String toString()
    {
        if(format == null || format.equals("")) return filename;
        if(format.equals(StsAbstractFile.binaryFormat))
            return name + " (binary)";
        else
            return name;
    }

    public static boolean isFilenameValid(String fullFilename)
    {
        String fname = null;
        for(int i=0; i<invalidCharacters.length; i++)
        {
            // Need to just check the filename since some characters are allowed in directory names that are not in filenames
            // ex. URL (which we use in StsFile) does not handle spaces in filenames but is okay if in directory name
            try {fname = StsFile.getFilenameFromPathname(fullFilename);}
            catch(Exception e) { System.out.println("Error parsing filename from fullname.");}
            if (fname.indexOf(invalidCharacters[i]) > 0)
            {
                new StsMessage(null, StsMessage.ERROR, "StsFile: Invalid character (" + invalidCharacters[i] + ") in name " + fullFilename);
                return false;
            }
        }
        return true;
    }

    public void parseFilename()
    {
        StsKeywordIO.parseAsciiFilename(filename);
        this.group = StsKeywordIO.group;
        this.format = StsKeywordIO.format;
        this.name = StsKeywordIO.name;
        this.subname = StsKeywordIO.subname;
        this.version = StsKeywordIO.version;
        if(this.format.equals(binaryFormat))
            fileType = binaryFormat;
        else
            fileType = group;
    }

    boolean nameSame(StsAbstractFile otherFilename)
    {
        return name.equals(otherFilename.name);
    }
}
