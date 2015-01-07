package com.Sts.IO;

import com.Sts.Utilities.*;

import javax.jnlp.*;
import java.net.*;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

public class StsWebStartJar extends StsJar
{
    private String codebase;
    private DownloadService downloadService;

    private StsWebStartJar(String jarFilename)
    {
        super("", jarFilename);
    }

    static public StsWebStartJar constructor(String jarFilename)
    {
        StsWebStartJar jar = new StsWebStartJar(jarFilename);
        if(!jar.initialize()) return null;
        return jar;
    }

    public boolean initialize()
//    public boolean classInitialize(boolean cacheIt)
    {
        try
        {
        /*
            URL codeBaseURL = JNLPUtilities.getBasicService().getCodeBase();
            codebase = codeBaseURL.toString();
            String jarPath = codebase + jarFilename;
            jarURL =  new URL(jarPath);
            System.out.println("URL for file " + jarFilename + ": " + jarPath);
        */
            downloadService = JNLPUtilities.getDownloadService();
            if(downloadService == null)
            {
                System.out.println("Couldn't find downloadService for javax.jnlp.DownloadService in StsSelectCurveTableMode.getResourceURL");
                return false;
            }
            URL codeBaseURL = JNLPUtilities.getBasicService().getCodeBase();
            codebase = codeBaseURL.toString();
            String jarPath = codebase + jarFilename;
            jarURL =  new URL(jarPath);
            boolean cached = downloadService.isResourceCached(jarURL, null);
            System.out.println("URL for file " + jarFilename + ": " + jarPath + " Cached: " + cached);
            if(!cached)
            {
                System.out.println("Not cached.");
                DownloadServiceListener downloadListener = downloadService.getDefaultProgressWindow();
                downloadService.loadResource(jarURL, null, downloadListener);
            }
            return createFileSet(jarPath);
        /*
            if((!cached) && (cacheIt))
            {
                System.out.println("Not cached and want to cache it.");
                DownloadServiceListener downloadListener = downloadService.getDefaultProgressWindow();
                downloadService.loadResource(jarURL, null, downloadListener);
                return createFileSet(jarPath);
            }
            else if(cached)
            {
                System.out.println("Cached");
                return createFileSet(jarPath);
            }
            else
            {
                System.out.println("Not cached and wont be cached.");
                return false;
            }
        */
//            return createFileSet(jarPath);
        }
        catch(Exception e)
        {
            String urlString;
            if(jarURL != null) urlString = jarURL.toString();
            else               urlString = "nullURL";
            StsException.outputException("StsWebStartJar.classInitialize() failed." +
                " Couldn't find file in cache and couldn't download it: " + urlString,
                e, StsException.WARNING);
            return false;
        }
    }

    public String getDescription()
    {
        return new String("WebStart downloaded jarfile: " + jarFilename);
    }

/*
    public URL getJarFileEntryURL(String entryName)
    {
        String urlEntryName;

        try
        {
            if(codebase == null) return null;
            urlEntryName = new String("jar:" + codebase + jarFilename + "!/" + entryName);

            System.out.println("StsWebStartJar.getJarFileEntryURL() urlEntryName: " + urlEntryName);
            return new URL(urlEntryName);
        }
        catch(Exception e)
        {
            StsException.outputException("StsWebStartJar.getJarFileEntryURL() failed." +
            "Couldn't find jar file entry: " + jarFilename + "/" + entryName,
            e, StsException.WARNING);
            return null;
        }
    }
*/
/*
    public String checkUnpack()
    {
        InputStream is = null;
        String directory = "dont-know";
        try
        {
            FileSaveService fss = JNLPUtilities.getFileSaveService();

            Set jarEntries = getJarFileEntries();
            int nJarEntries = jarEntries.size();
            if(nJarEntries == 0) return directory;
            Iterator iter = jarEntries.iterator();
            while(iter.hasNext())
            {
                String entryName = (String)iter.next();
                URL url = getJarFileEntryURL(entryName);
                URLConnection urlConnection = url.openConnection();
                is = urlConnection.getInputStream();
//                bis = new BufferedInputStream(is);
                FileContents fileContents = fss.saveFileDialog(null, null, is, entryName);
                String filename = fileContents.getName();
                System.out.println("WebStart file write to: " + filename);
            }
            return directory;
        }
        catch(Exception e)
        {
            StsException.outputException("StsJar.checkUnpack() failed.",
                e, StsException.WARNING);
            return directory;
        }
        finally
        {
            try
            {
                if(is != null) is.close();
            }
            catch(Exception e) { }
        }
    }
*/
}
