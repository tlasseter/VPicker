package com.Sts.IO;

import com.Sts.Utilities.*;

import java.util.*;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

/** Defines a set of files of a certain type and origin.  Subclassed by StsFileSet, StsJar, and StsWebStartFileSet. */
abstract public class StsAbstractFileSet
{
    protected String directory;
    protected StsAbstractFile[] files = new StsAbstractFile[0];

    abstract public StsAbstractFile[] initializeAvailableFiles();
    abstract public String getDescription();

    public StsAbstractFileSet()
    {
    }

    public int size()
    {
        if(files == null) return 0;
        else return files.length;
    }

    public void setDirectory(String directory)
    {
        this.directory = directory;
    }

    public StsAbstractFile[] getFiles()
    {
        return files;
    }

    public String[] getFilenames()
    {
        int nFilenames = files.length;
        String[] filenames = new String[nFilenames];
        for(int n = 0; n < nFilenames; n++)
            filenames[n] = files[n].getFilename();
        return filenames;
    }


    public StsAbstractFile getFileEndingWith(String suffix)
    {
        String[] filenames = getFilenames();
        for(int i=0; i<filenames.length; i++)
        {
            if(filenames[i] == null) continue;
            if(filenames[i].endsWith(suffix))
                return getFile(i);
        }
        return null;
    }

/*
    public URL getURL(int index)
    {
        try
        {
            return files[index].getURL();
        }
        catch(Exception e)
        {
            StsException.outputException("StsFileSet.getURL(index) failed." +
                "directory: " + directory + " index: " + index,
                e, StsException.WARNING);
            return null;
        }
    }
*/
    public int getIndex(String filename)
    {
        if(files == null)
        {
            StsException.systemError("StsAbstractFileSet.index() failed. No files are available.");
            return -1;
        }
        int nFiles = files.length;
        for(int n = 0; n < nFiles; n++)
            if(filename.equals(files[n].getFilename())) return n;

        StsException.systemError("StsAbstractFileSet.index() failed. No file found matching: " + filename);
        return -1;
    }

	 public int getIndexFromName(String name)
	 {
		  if(files == null)
		  {
				StsException.systemError("StsAbstractFileSet.getIndexFromName() failed. No files are available.");
				return -1;
		  }
		  int nFiles = files.length;
		  for(int n = 0; n < nFiles; n++)
		  {
			  if (files[n].getFilename().indexOf(name) != -1)  return n;
		  }
		  StsException.systemError("StsAbstractFileSet.getIndexFromName() failed. No file found matching: " + name);
		  return -1;
	 }

	 public int getIndexFromName(String name, StsFilenameFilterFace filenameFilter)
	 {
		  if(files == null)
		  {
				StsException.systemError("StsAbstractFileSet.getIndexFromName() failed. No files are available.");
				return -1;
		  }
		  int nFiles = files.length;
		  for(int n = 0; n < nFiles; n++)
		  {
			  String stemName = filenameFilter.getFilenameName(files[n].getFilename());
			  if (stemName.equals(name))
			  {
				  return n;
			  }
			  //if (files[n].createFilename().indexOf(name) != -1)  return n;
		  }
		  StsException.systemError("StsAbstractFileSet.getIndexFromName() failed. No file found matching: " + name);
		  return -1;
	 }

    public StsAbstractFile getFile(String filename)
    {
        int index = getIndex(filename);
        if(index < 0) return null;
        return getFile(index);
    }

	 public StsAbstractFile getFileFromName(String name)
	 {
		  int index = getIndexFromName(name);
		  if(index < 0) return null;
		  return getFile(index);
	 }

	 public StsAbstractFile getFileFromName(String name, StsFilenameFilterFace filenameFilter)
	 {
		  int index = getIndexFromName(name, filenameFilter);
		  if(index < 0) return null;
		  return getFile(index);
	 }

    public StsAbstractFile getFile(int index)
    {
        if(files == null || index < 0 || index >= files.length) return null;
        return files[index];
    }

    public void sort()
    {
        if(files == null) return;
        int nFiles = files.length;
        if(nFiles < 2) return;
        String[] filenames = new String[nFiles];
        for(int n = 0; n < nFiles; n++)
            filenames[n] = files[n].getFilename();
        StsMath.qsort(filenames, files);
    }

    static ArrayList filteredClassNames(String[] filenames, String packageName)
    {
        ArrayList classNamesArray = new ArrayList();

        if(filenames == null) return classNamesArray;

        for(int n = 0; n < filenames.length; n++)
        {
            if(!filenames[n].endsWith(".class")) continue;
            String packageString = filenames[n].replace('/', '.');
            if(!packageString.startsWith(packageName)) continue;
            classNamesArray.add(filenames[n]);
        }
        return classNamesArray;
    }
}
