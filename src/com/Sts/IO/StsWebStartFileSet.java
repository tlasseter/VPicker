package com.Sts.IO;

import com.Sts.Utilities.*;

import javax.jnlp.*;
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

public class StsWebStartFileSet extends StsAbstractFileSet
{
    URL codebase;
    URL directoryURL; // files are in this server-side directory
    int[] tags;
    boolean cacheFiles = true;

    static final private boolean debug = false;

    // subdirectory must be terminated with a "/"
    // if there is no subdirectory, specifiy as ""
    // codebase terminates in a slash
    private StsWebStartFileSet(String subdirectory, String[] filenames, boolean cacheFiles) throws MalformedURLException, IOException
    {
        this.cacheFiles = cacheFiles;
        codebase = JNLPUtilities.getCodeBase();
        if(debug) System.out.println("StsWebStartFileSet.codebase: " + codebase.toString());

        directory = new String(codebase.toString() + subdirectory);
        directoryURL = new URL(directory);
        if(debug) System.out.println("StsWebStartFileSet.directoryURL: " + directory);

        // set the appropriate tag for each file we might be interested in.
        // Note that persistenceService.getTag(muffinURL) throws an exception if muffin does not exist
        // on the client-side. So classInitialize to NOT_CACHED and change tags of those found on the
        // client side.
        int nFilenames = filenames.length;
        files = new StsWebStartFile[nFilenames];
        for(int n = 0; n < nFilenames; n++)
            files[n] = StsWebStartFile.constructor(directory, filenames[n], cacheFiles);

        if(!cacheFiles) return;

        PersistenceService persistenceService = JNLPUtilities.getPersistenceService();

        try
        {
            String[] muffins = null;
            // muffins are file pathnames relative to the codebase
            // first time thru, this fails so ignore and check if muffins != null
            try { muffins = persistenceService.getNames(directoryURL); }
            catch(Exception e) { }

            if(muffins != null)
            {
                int nMuffins = muffins.length;
                if(debug) System.out.println("Number of muffins: " + nMuffins);

                // skip the first which is the subdirectory itself
                for(int c = 0; c < nMuffins; c++)
                {
                    String muffin = muffins[c];
                    if(debug) System.out.print("DEBUG: StsWebStartFileSet.muffins[" + c + "]:" + muffin);
                    URL muffinURL = new URL(directoryURL, muffin);
                    int tag = persistenceService.getTag(muffinURL);
                    if(debug) System.out.println("tag: " + StsWebStartFile.tagLabels[tag]);
                    for(int f = 0; f < nFilenames; f++)
                    {
                        StsWebStartFile file = (StsWebStartFile)files[f];
                        if(muffin.equals(file.getFilename()))
                        {
                            file.setTag(tag);
                            file.setFileContents(persistenceService.get(muffinURL));
                            break;
                        }
                    }
                }
            }
        }
        // mainDebug print
//        for(int f = 0; f < nFilenames; f++)
//            System.out.println("DEBUG: StsWebStartFileSet.filenames[" + f + "] name: " + filenames[f] + " tag:" + tagLabels[tags[f]]);
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFileSet.init() failed.",
                e, StsException.WARNING);
        }
    }

    static public StsWebStartFileSet constructor(String subdirectory, String[] filenames, boolean cacheFiles)
    {
        try
        {
            return new StsWebStartFileSet(subdirectory, filenames, cacheFiles);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFileSet.constructor() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public String getDescription()
    {
        return new String("WebStart files from: " + directory);
    }

    /** Doesn't make sense to initialize (again).  Here for superClass compatability. */
    public StsAbstractFile[] initializeAvailableFiles()
    {
        return files;   
    }
/*
    public Set getFilenameSet(String subdirectory, FilenameFilter filter)
    {
        try
        {
            URL url;

            URL codebase = JNLPUtilities.getCodeBase();
            if(subdirectory != null)
                url = new URL(codebase, subdirectory);
            else
                url = codebase;

            System.out.println("Getting server files from URL: " + url.toString());

            filenames = persistenceService.getNames(url);

            if(filenames != null)
            {
                System.out.println("Found files at URL: " + url.toString());
                for(int n = 0; n < filenames.length; n++)
                    System.out.println("file " + n + " : " + filenames[n]);
            }

            java.util.List list = Arrays.asList(filenames);
            return new HashSet(list);
*/
    /*
            String[] filteredFilenames = String.
            // get the attributes (tags) for each of these muffins.
            // update the server's copy of the data if any muffins
            // are dirty
            int [] tags = new int[muffins.length];
            URL [] muffinURLs = new URL[muffins.length];
            for(int i = 0; i < muffins.length; i++)
            {
                muffinURLs[i] = new URL(codebase.toString() + muffins[i]);
                tags[i] = persistenceService.getTag(muffinURLs[i]);
                // update the server if anything is tagged DIRTY
                if(tags[i] == PersistenceService.DIRTY) doUpdateServer(muffinURLs[i]);

            }

            // read in the contents of a muffin and then delete it
            FileContents fc = persistenceService.get(muffinURLs[0]);
            long maxsize = fc.getMaxLength();
            int length = (int)fc.getLength();
            byte [] buf = new byte[length];
            InputStream is = fc.getInputStream();
            int pos = 0;
            while((pos = is.read(buf, pos, length - pos)) > 0) {
                // just loop
            }
            is.close();

            persistenceService.delete(muffinURLs[0]);

            // re-create the muffin and repopulate its data
            persistenceService.create(muffinURLs[0], maxsize);
            fc = persistenceService.get(muffinURLs[0]);
            // don't append
            OutputStream os = fc.getOutputStream(false);
            os.write(buf);
            os.close();
            return true;
    */
/*
        }
        catch (Exception e)
        {
            StsException.outputException("StsWebStartFileSet.loadFileResources() failed.",
                e, StsException.WARNING);
            return null;
        }
    }
*/
    public StsAbstractFile getFile(int index)
    {
        try
        {
            StsWebStartFile file = (StsWebStartFile)files[index];
            if(!file.checkCache()) return null;
            return file;
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFileSet.getFile(index) failed." +
                "index: " + index,
                e, StsException.WARNING);
            return null;
        }
    }

    private void debugCheckMuffins(PersistenceService persistenceService)
    {
        try
        {
            if(debug) System.out.println("debugCheckMuffins.directory: " + directoryURL.toString());
//            PersistenceService persistenceService = JNLPUtilities.getPersistenceService();

            // muffins are file pathnames relative to the codebase
            String[] muffins = persistenceService.getNames(directoryURL);
            if(muffins == null)
            {
                System.out.println("debugCheckMuffins: Sorry no muffins");
                return;
            }
            int nMuffins = muffins.length;
            if(debug) System.out.println("debugCheckMuffins: Number of muffins: " + nMuffins);

            for(int c = 0; c < nMuffins; c++)
            {
                String muffin = muffins[c];
                if(debug) System.out.print("  muffins[" + c + "]:" + muffin);
                URL muffinURL = new URL(directoryURL, muffin);
                int tag = persistenceService.getTag(muffinURL);
                if(debug) System.out.println("  tag: " + StsWebStartFile.tagLabels[tag]);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartFileSet.debugCheckMuffins() failed.",
                e, StsException.WARNING);
        }
    }
}
