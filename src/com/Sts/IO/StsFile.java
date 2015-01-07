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

import javax.swing.*;
import java.awt.*;
import java.io.*;
import java.net.*;
import java.nio.channels.*;
import java.util.*;

public class StsFile extends StsAbstractFile
{
	StsFilenameEndingFilter filenameFilter = null;
    RandomAccessFile raf;

    protected StsFile(String directory, String filename) throws URISyntaxException, MalformedURLException, StsException
	{
        this.filename = filename;
		String pathname = new String(directory + filename);
        constructURL(pathname);
        parseFilename();
    }

    protected StsFile(String directory, String filename, StsFilenameEndingFilter filenameFilter) throws URISyntaxException, MalformedURLException, StsException
	{
        this.filename = filename;
		this.filenameFilter = filenameFilter;
        String pathname = filename;
        if(!filename.contains(directory.substring(0,directory.lastIndexOf(File.separator))))
		    pathname = new String(directory + filename);
        constructURL(pathname);
        parseFilename();
    }

    private StsFile(String pathname) throws URISyntaxException, MalformedURLException, StsException
    {
        filename = getFilenameFromPathname(pathname);
        constructURL(pathname);
        parseFilename();
    }

    private void constructURL(String pathname) throws URISyntaxException, MalformedURLException,  StsException
    {
        if(pathname.startsWith("file:"))
            pathname = StsStringUtils.trimPrefix(pathname, "file:");
        File file = new File(pathname);
        url = file.toURI().toURL();
//        url = new URL("file", "", pathname);
    }

	static public StsFile constructor(String directory, String filename)
	{
		return constructor(directory, filename, null);
	}

    static public StsFile constructor(String directory, String filename, StsFilenameEndingFilter filenameFilter)
    {
        try
        {
            if(filename == null) return null;
            if(!isFilenameValid(filename))
                return null;

            if(!directory.endsWith(File.separator))
                return new StsFile(directory + File.separator, filename, filenameFilter);
            else
                return new StsFile(directory, filename, filenameFilter);
        }
        catch(Exception e)
        {
            StsException.outputException("StsFile.constructor(directory, filename) failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public String getDirectory()
    {
        try
        {
            return StsFile.getDirectoryFromPathname(this.getPathname());
        }
        catch(Exception e)
        {
            StsException.systemError("StsFile.getParent() failed.");
            return null;
        }
    }

    static public StsFile constructor(String pathname)
    {
        try
        {
            String filename = getFilenameFromPathname(pathname);
            if(!isFilenameValid(filename))
                return null;
            return new StsFile(pathname);
        }
        catch(Exception e)
        {
            StsException.outputException("StsFile.constructor(pathname) failed.",
                e, StsException.WARNING);
            return null;
        }
    }

	static public DataInputStream constructDIS(String pathname)
	{
		try
		{
            String filename = getFilenameFromPathname(pathname);
            if(!isFilenameValid(filename)) return null;
            StsFile file = new StsFile(pathname);
			return file.getDataInputStream();
		}
		catch(Exception e)
		{
			StsException.outputException("StsFile.constructor(pathname) failed.",
				e, StsException.WARNING);
			return null;
		}
	}

	static public DataOutputStream constructDOS(String pathname)
	{
		try
		{
            String filename = getFilenameFromPathname(pathname);
            if(!isFilenameValid(filename)) return null;
            StsFile file = new StsFile(pathname);
			return file.getDataOutputStream();
		}
		catch(Exception e)
		{
			StsException.outputException("StsFile.constructor(pathname) failed.",
				e, StsException.WARNING);
			return null;
		}
    }
    static public String getFilenameFromPathname(String pathname) throws StsException
    {
        String fileSeparator = File.separator;
        int separatorIndex = pathname.lastIndexOf(fileSeparator);
        if(separatorIndex == -1 && !fileSeparator.equals("/"))
            separatorIndex = pathname.lastIndexOf("/");
        if(separatorIndex == -1)
        {
            //throw new StsException(StsException.WARNING, "StsFile(pathname) failed. Didn't find separator.");
            return pathname; // May just be a filename with no path.
        }
        int length = pathname.length();
        return pathname.substring(separatorIndex+1, length);
    }

    static public String getParentDirectoryFromPathname(String pathname)
    {
        File file = new File(pathname);
        return file.getParent();
    }

    // Find last fileSeparator, but ignore the one at the end of the string if it exists.
    // Return string up to but not including this separator.
    static public String getDirectoryFromPathname(String pathname) throws StsException
    {
        return getDirectoryFromPathname(pathname, false);
    }

    public InputStream getInputStream() throws IOException
    {
		return new FileInputStream(getFile());
    }

      // Find last fileSeparator, optionally ignore the one at the end of the string if it exists.
    // Return string up to but not including this separator.
    static public String getDirectoryFromPathname(String pathname, boolean includeLast) throws StsException
    {
        if(pathname == null) return null;
        int length = pathname.length();
        String fileSeparator = File.separator;
        int separatorIndex = 0;
        if(!includeLast)
        {
            separatorIndex = pathname.lastIndexOf(fileSeparator, length-2);
            if(separatorIndex == -1 && !fileSeparator.equals("/"))
                separatorIndex = pathname.lastIndexOf("/", length-2);
        }
        else
        {
            separatorIndex = pathname.lastIndexOf(fileSeparator, length);
            if(separatorIndex == -1 && !fileSeparator.equals("/"))
                separatorIndex = pathname.lastIndexOf("/", length); 
        }
        if(separatorIndex == -1) throw new StsException(StsException.WARNING, "StsFile(pathname) failed. Didn't find separator.");
        return pathname.substring(0, separatorIndex+1);
    }

	public FileInputStream getFileInputStream()
	{
		try
		{
			return (FileInputStream)url.openStream();
		}
		catch(Exception e)
		{
            StsException.systemError(this, "getFileInputStream", "file not found: " + filename);
			return null;
		}
    }
    public InputStream getMonitoredInputStream(Component parentComponent)
    {
        try
        {
            InputStream is = url.openStream();
            return new BufferedInputStream(
                    new ProgressMonitorInputStream(parentComponent, "Reading " + url.getPath(), is));
        }
        catch(Exception e)
        {
//            StsException.systemError("StsFile.getInputStream() failed, file not found: " + filename);
            return null;
        }
    }


	public DataInputStream getDataInputStream() throws IOException
	{
		InputStream is = getInputStream(); // true: append write to end of file
		BufferedInputStream bis = new BufferedInputStream(is);
		return new DataInputStream(bis);
	}


    public OutputStream getOutputStream() throws FileNotFoundException
    {
        return getOutputStream(false);
    }

    public OutputStream getOutputStream(boolean append) throws FileNotFoundException
    {
        File file = getFile();
        if(file == null) return null;
        return new FileOutputStream(file, append);
    }

    public OutputStream getOutputStreamAndPosition(long position) throws FileNotFoundException
    {
        File file = getFile();
        try
        {
            FileOutputStream outputStream = new FileOutputStream(file, true);
            FileChannel channel = outputStream.getChannel();
            channel.position(position);
            long filePosition = channel.position();
            if(filePosition != position)
            {
                StsException.systemError(this, "openWrite(position", "Failed to set file position: " + filePosition + " to desired position: " + position);
            }
            System.out.println("FILE DEBUG: Open output file and position " + filename +
								   " position: " + filePosition + " length: " + file.length() +
								   " thread: " + Thread.currentThread().getName() +
								   " time: " + System.currentTimeMillis());
            return outputStream;
        }
        catch(Exception e)
        {
            return getOutputStream(true);
        }
    }

    public DataOutputStream getDataOutputStream() throws FileNotFoundException
	{
		OutputStream os = getOutputStream(true); // true: append write to end of file
		BufferedOutputStream bos = new BufferedOutputStream(os);
		return new DataOutputStream(bos);
	}

    static public boolean checkCreateDirectory(String directory)
    {
        StsFile file = StsFile.constructor(directory);   
        if(file.exists()) return true;
        try
        {
            file.getFile().mkdirs();
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(StsFile.class, "checkCreateDirectory", e);
            return false;
        }
    }

    public boolean delete()
	{
         File file = getFile();
         if(file == null) return false;
         return file.delete();
	 }

	 public void writeStringsToFile(String[] strings)
	 {
		 StsAsciiFile asciiFile = null;

		 try
		 {
			 if(strings == null) return;
			 asciiFile = new StsAsciiFile(this);
			 if(!asciiFile.openWrite()) return;
			 for(int n = 0; n < strings.length; n++)
				 asciiFile.writeLine(strings[n]);
			 asciiFile.close();
		 }
		 catch(Exception e)
		 {
			 StsException.outputException("StsProject.writeStringsToFile() failed.",
					 e, StsException.WARNING);
		 }
		 finally
		 {
			 if(asciiFile != null) asciiFile.close();
		 }
	 }

	 public String[] readStringsFromFile() throws IOException
	 {
		 String[] strings = new String[0];
		 String string;

		 StsAsciiFile asciiFile = new StsAsciiFile(this);
		 if(!asciiFile.openRead()) return strings;
		 try
		 {
			 while ((string = asciiFile.readLine()) != null)
			 {
				 strings = (String[])StsMath.arrayAddElement(strings, string);
			 }
			 return strings;
		 }
		 catch(Exception e)
		 {
			 StsException.systemError("Failed to read " + getPathname());
			 return strings;
		 }
		 finally
		 {
			 if (asciiFile != null) asciiFile.close();
		 }
	 }

	public void setFilenameFilter(StsFilenameEndingFilter filenameFilter)
	{
		this.filenameFilter = filenameFilter;
	}

	public String getFilenameStem()
	{
		if(filenameFilter == null)
			return filename;
		else
			return filenameFilter.getFilenameName(filename);
	}

    static public boolean copy(String source, String destination) throws IOException
    {
        File testFile = new File(source);
        if(!testFile.exists())
            return false;

        FileChannel srcChannel = new FileInputStream(source).getChannel();

        // Create channel on the destination
        FileChannel dstChannel = new FileOutputStream(destination).getChannel();

        // Copy file contents from source to destination
        dstChannel.transferFrom(srcChannel, 0, srcChannel.size());

        // Close the channels
        srcChannel.close();
        dstChannel.close();

        return true;
    }

    public boolean copy(String destination) throws IOException
    {
        FileChannel srcChannel = new FileInputStream(this.getPathname()).getChannel();

        // Create channel on the destination
        FileChannel dstChannel = new FileOutputStream(destination).getChannel();

        // Copy file contents from source to destination
        dstChannel.transferFrom(srcChannel, 0, srcChannel.size());

        // Close the channels
        srcChannel.close();
        dstChannel.close();

        return true;
    }

    void removeLastRepeat(ArrayList<StsFile> files)
    {
        for (StsAbstractFile file : files)
            if (nameSame(file)) return;
        files.add(this);
    }

    void removeFirstRepeat(ArrayList<StsFile> files)
    {
        for (StsAbstractFile file : files)
        {
            if (nameSame(file))
            {
                files.remove(file);
                files.add(this);
                return;
            }
        }
        files.add(this);
    }

    static public void clearDirectoryAndFiles(String pathname)
    {
        clearDirectoryAndFiles(new File(pathname));
    }

    /** recursively removes all files in this directory tree */
    static public void clearDirectoryAndFiles(File parentFile)
    {
        if(!parentFile.exists()) return;
        if(parentFile.isDirectory())
        {
            File[] files = parentFile.listFiles();
            for(File file : files)
                clearDirectoryAndFiles(file);
        }
        parentFile.delete();
    }

	public boolean isWritable() { return true; }
    static public void main(String[] args)
    {
        try
        {
            String s = File.separator;
            String dir = "c:";
            String one = "data";
            String two = "Q" + s + "WellsTest";
            String filename = "well-dev.txt.1G#03";

            String name, parent, file;

            name = dir + s + one + s + two + s + filename;
            StsFile ff = StsFile.constructor(name);
        /*
            file = getFilenameFromPathname(name);
            parent = getDirectoryFromPathname(name);
            name = dir + s + one + s + two + s;
            parent = getDirectoryFromPathname(name);
            name = dir + s + one + s + two;
            parent = getDirectoryFromPathname(name);
            name = dir + s;
        */
        }
        catch(Exception e) { }
    }
}