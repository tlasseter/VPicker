package com.Sts.IO;

import com.Sts.Utilities.*;

import java.io.*;
import java.net.*;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

/** This is a concrete set of StsFile(s) in this directory which satisfy a filenameFilter or filters */
public class StsFileSet extends StsAbstractFileSet
{
    FilenameFilter[] filenameFilters;

    private StsFileSet(String directory) throws IOException
    {
        this.directory = directory;
    }

    static public StsFileSet constructor(String directory, FilenameFilter filenameFilter)
    {
        return constructor(directory, new FilenameFilter[] { filenameFilter });
    }

    static public StsFileSet constructor(String directory, FilenameFilter[] filenameFilters)
    {
        try
        {
            StsFileSet fileset = new StsFileSet(directory);
            fileset.filenameFilters = filenameFilters;
            fileset.addFiles();
            return fileset;
        }
        catch(Exception e)
        {
            StsException.outputException("StsFileSet.constructor() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public StsAbstractFile[] initializeAvailableFiles()
    {
        files = new StsAbstractFile[0];
        addFiles();
        return files;
    }

    public void addFiles()
    {
        for(FilenameFilter filenameFilter : filenameFilters)
            addFiles(filenameFilter);
    }

    public void addFiles(FilenameFilter filenameFilter)
    {
        File directoryFile = new File(directory);
        String[] filenames = directoryFile.list(filenameFilter);
        if(filenames == null) return;

        StsAbstractFile[] newFiles;
        int nFilenames = filenames.length;
        int nFiles = 0;
        if(filenameFilter instanceof StsFilenameEndingFilter)
        {
            files = new StsFile[nFilenames];
            for(int n = 0; n < nFilenames; n++)
            {
                StsFile file = StsFile.constructor(directory, filenames[n], (StsFilenameEndingFilter) filenameFilter);
                if(file.isAFileAndNotDirectory()) files[nFiles++] = file;
            }
        }
        else
        {
            newFiles = new StsFile[nFilenames];
            for(int n = 0; n < nFilenames; n++)
            {
                StsFile file = StsFile.constructor(directory, filenames[n]);
                if(file.isAFileAndNotDirectory()) newFiles[nFiles++] = file;
            }
            files = (StsAbstractFile[])StsMath.arrayAddArray(files, newFiles, nFiles);
        }
        files = (StsAbstractFile[])StsMath.trimArray(files, nFiles);
    }

    public String getDescription()
    {
        return new String("files in directory: " + directory);
    }

    public StsAbstractFile getFile(int index)
    {
        try
        {
            return files[index];
        }
        catch(Exception e)
        {
            StsException.outputException("StsFileSet.getStsFile() failed.",
                e, StsException.WARNING);
            return null;
        }
    }
}
