package com.Sts.IO;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Utilities.*;

import javax.jnlp.*;
import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.net.*;

public class StsWebStartFile extends StsAbstractFile
{
    private int tag = NOT_CACHED; // file is NOT_CACHED, CACHED, DIRTY, or TEMPORARY
    private FileContents fileContents;

    // first 3 tag types are JavaWebStart definitions, so don't change.
    // DIRTY and TEMPORARY are not currently used by S2S
    static final int CACHED = PersistenceService.CACHED;
    static final int DIRTY = PersistenceService.DIRTY;
    static final int TEMPORARY = PersistenceService.TEMPORARY;
    static final int NOT_CACHED = 3;
    static final int DONT_CACHE = 4;
    static public final String[] tagLabels = new String[] { "CACHED", "DIRTY", "TEMPORARY", "NOT_CACHED", "DONT_CACHE" };

    static final int bufSize = 2048;

    static private final boolean debug = false;

    private StsWebStartFile(String directory, String filename, boolean cacheFiles) throws MalformedURLException, IOException
    {
        if(cacheFiles) tag = NOT_CACHED;
        else           tag = DONT_CACHE;

        this.filename = filename;
        url = new URL(directory + filename);
        parseFilename();
        fileType = WEBJAR;
    }

    static public StsWebStartFile constructor(String directory, String filename, boolean cacheFiles)
    {
        try
        {
            return new StsWebStartFile(directory, filename, cacheFiles);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFile.constructor() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void setTag(int tag) { this.tag = tag; }
    public int getTag() { return tag; }

    public void setFileContents(FileContents fileContents) { this.fileContents = fileContents; }

    public InputStream getInputStream()
    {
        try
        {
            if(tag == CACHED)
            {
                if(debug) System.out.println("FileContents: " + fileContents.getName());
                return fileContents.getInputStream();
            }
            else if(tag == DONT_CACHE)
            {
                return url.openStream();
            }
            else
                return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFile.getInputStream() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public InputStream getMonitoredInputStream(Component parentComponent)
    {
        InputStream is;

        try
        {
            String filename = fileContents.getName();
            if(tag == CACHED)
            {
                if(debug) System.out.println("FileContents: " + filename);
                is = fileContents.getInputStream();
            }
            else if(tag == DONT_CACHE)
                is = url.openStream();
            else
                return null;

            return new BufferedInputStream(
                    new ProgressMonitorInputStream(parentComponent, "Reading " + filename, is));
        }
        catch(Exception e)
        {
//            StsException.systemError("StsFile.getInputStream() failed, file not found: " + filename);
            return null;
        }
    }

    public OutputStream getOutputStream(boolean append)
    {
        try
        {
            if(tag == CACHED)
            {
                return fileContents.getOutputStream(true);
            }
            else if(tag == DONT_CACHE)
            {
                URLConnection connection = url.openConnection();
                return connection.getOutputStream();
            }
            else
                return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFile.getOutputStream() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public OutputStream getOutputStreamAndPosition(long position)
    {
        return getOutputStream(true);
    }

    public boolean checkCache()
    {
        try
        {
            if(tag == DONT_CACHE) return true;

            if(debug)
            {
                System.out.print("DEBUG: StsWebStartFileSet.getURL() URL: " + url.toString());
                System.out.println("    tag: " + tagLabels[tag]);
            }

            PersistenceService persistenceService = JNLPUtilities.getPersistenceService();

            if(tag == NOT_CACHED)
            {
                URLConnection connection = url.openConnection();

                 // create the cached URL
                int length = connection.getContentLength();
                persistenceService.create(url, length);
                fileContents = persistenceService.get(url);

                // copy data from server URL to cache
                OutputStream os = fileContents.getOutputStream(true);

               InputStream is = connection.getInputStream();
                byte[] buf = new byte[bufSize];
                BufferedInputStream bis = new BufferedInputStream(is);
                int len = 0;
                while( (len = bis.read(buf, 0, bufSize)) != -1)
                {
                    os.write(buf, 0, len);
                }
                is.close();
                os.flush();
                os.close();

                tag = persistenceService.getTag(url);
                setTag(tag);
                if(debug) System.out.print("  " + length + " bytes has been cached; tag set to " + tagLabels[tag]);
                return true;
            }
            else if(tag == CACHED)
            {
                if(debug) System.out.print("  " + " currently cached.");
                if(fileContents == null) fileContents = persistenceService.get(url); // doublecheck: should have been set in StsWebStartFileSet
                return true;
            }
            else
                return false;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFileSet.getURL(index) failed." +
                "filename: " + filename,
                e, StsException.WARNING);
            return false;
        }
    }
	public boolean isWritable() { return false; }
/*
    public JNLPRandomAccessFile getRandomAccessFile(String mode)
    {
        try
        {
            if(tag == CACHED)
                return fileContents.getRandomAccessFile(mode);
            else if(tag == DONT_CACHE)
            {
                String filePath = url.toString();
                if(debug) System.out.println("StsWebStartFile.getRandomAccessFile() file: " + filePath);
                return new StsRandomAccessFile(filePath, mode);
            }
            else
                return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFile.getRandomAccessFile() failed.",
                e, StsException.WARNING);
            return null;
        }
    }
*/
}
