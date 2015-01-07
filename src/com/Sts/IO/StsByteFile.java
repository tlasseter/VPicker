package com.Sts.IO;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

public class StsByteFile
{
    StsAbstractFile file = null;
    String filename = null;
    InputStream is = null;
    OutputStream os = null;

    public StsByteFile()
    {
    }

    public StsByteFile(StsAbstractFile file)
    {
        this.file = file;
        filename = file.getFilename();
    }


    public boolean openRead(Component parentComponent)
    {
        try
        {
            if(is != null) close();
            if(parentComponent == null)
                is = file.getInputStream();
            else
                is = file.getMonitoredInputStream(parentComponent);

            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsByteFile.openReadAndCheck() failed." +
                "Can't read: " + file.getFilename(), e, StsException.WARNING);
            return false;
        }

    }

    public boolean openWrite()
    {
        try
        {
            if(os != null) close();
            os = file.getOutputStream(true); // true: append write to end of file
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsByteFile.openWrite() failed." +
                "Can't write: " + filename, e, StsException.WARNING);
            return false;
        }
    }

    public boolean openReadWrite()
    {
        return openRead(null) && openWrite();
    }

    /** close this binary file */
    public boolean close()
    {
        try
        {
            if (os != null)
            {
                os.flush();
                os.close();
                os = null;
            }
            if (is != null)
            {
                is.close();
                is = null;
            }
			return true;
        }
        catch (Exception e)
        {
            StsException.outputException("StsByteFile.close() failed."
                    + "Unable to close file " + filename, e, StsException.WARNING);
			return false;
        }
    }

    public byte[] getBytes(int size) throws IOException
    {
        byte[] bytes = new byte[size];
        is.read(bytes);
        return bytes;
    }

    // returns number of bytes read; returns -1 if EOF
    public int read(byte[] bytes) throws IOException
    {
        return is.read(bytes);
    }

    public long skip(int size) throws IOException
    {
        try {  return is.skip(size); }
        catch(Exception e) { return 0; }
    }

    public byte[] getMonitoredBytes(int size, Frame frame) throws IOException
    {
        ProgressMonitorInputStream pmis;

        byte[] bytes = new byte[size];
        pmis = new ProgressMonitorInputStream(frame, "Reading " + file.getPathname(), is);
        int nRead = pmis.read(bytes);
        if(nRead != size)
        {
            new StsMessage(frame, StsMessage.WARNING, "For file: " + file.getPathname() +
                ", only " + nRead + " bytes read. Expected to read " + size);
            return null;
        }
        return bytes;
    }

    public void write(byte[] bytes) throws IOException
    {
        os.write(bytes);
    }
}