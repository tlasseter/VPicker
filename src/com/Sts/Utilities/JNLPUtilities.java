package com.Sts.Utilities;

import javax.jnlp.*;
import java.io.*;
import java.net.*;

// import com.sun.jnlp.*;
// import com.sun.javaws.jnl.*;

/**
 * Chapter 11 - Utilities
 * @author Mauro Marinilli
 * @version 1.0
 */
public class JNLPUtilities
{
    static private final boolean debug = false;

    /** JNLP services */
    static private ExtensionInstallerService eiService;
    static private BasicService basicService;
    static private DownloadService downloadService;

    /** Initialize basic services class */
    static private synchronized void initServices()
    {
        if (basicService == null)
        {
	        try
            {
//		        eiService = (ExtensionInstallerService)ServiceManager.lookup("javax.jnlp.ExtensionInstallerService");
		        basicService = (BasicService)ServiceManager.lookup("javax.jnlp.BasicService");
		        downloadService = (DownloadService)ServiceManager.lookup("javax.jnlp.DownloadService");
	        }
            catch(UnavailableServiceException use)
            {
		        System.out.println("Unable to locate service: " + use);
	        }
        }
    }

    static public BasicService getBasicService()
    {
	    if (basicService == null) initServices();
	    return basicService;
    }

    static public ExtensionInstallerService getInstallService()
    {
	    if (eiService == null) initServices();
	    return eiService;
    }

    static public DownloadService getDownloadService()
    {
	    if (downloadService == null) initServices();
	    return downloadService;
    }

    static public Object getService(String fullyQName)
    {
        Object service = null;
        try
        {
            service = ServiceManager.lookup(fullyQName);
        }
        catch (UnavailableServiceException use)
        {
            StsException.outputException("JNLPUtilities.getService(" + fullyQName + ") failed.",
                use, StsException.WARNING);
        }
        return service;
    }

    static public FileSaveService getFileSaveService()
    {
        return (FileSaveService)getService("javax.jnlp.FileSaveService");
    }

    static public URL getCodeBase()
    {
        BasicService basicService = getBasicService();
        return basicService.getCodeBase();
    }

    static public PersistenceService getPersistenceService()
    {
        return (PersistenceService)getService("javax.jnlp.PersistenceService");
    }

    static public final URL getJarURL(String jarname)
    {
        URL url = null;
        try
        {
            URL codebase = getCodeBase();
            String codebaseName = codebase.toString();
            if(debug) System.out.println("jar:" + codebaseName + jarname);

            url = new URL("jar:" + codebaseName + jarname + "!/");
            return url;
        }
        catch(Exception e)
        {
            String urlString;
            if(url != null) urlString = url.toString();
            else            urlString = "nullURL";
            StsException.outputException("JNLPUtilities.getWebStartJarURL() failed." +
                " Couldn't find file in cache and couldn't download it: " + urlString,
                e, StsException.WARNING);
            return null;
        }
    }

    static public final boolean isJarCached(URL url)
    {
        DownloadService downloadService;
        try
        {
            downloadService = (DownloadService)getService("javax.jnlp.DownloadService");
            if(downloadService == null)
            {
                StsException.systemError("JNLPUtilities.isJarCached() failed. Couldn't find downloadService");
                return false;
            }
            return downloadService.isResourceCached(url, null);
        }
        catch(Exception e)
        {
            String urlString;
            if(url != null) urlString = url.toString();
            else            urlString = "nullURL";
            StsException.outputException("JNLPUtilities.isJarCached() failed for file: " + urlString,
                e, StsException.WARNING);
            return false;
        }
    }

    static public final boolean isResourceCached(URL url)
    {
        DownloadService downloadService;
        try
        {
            downloadService = (DownloadService)getService("javax.jnlp.DownloadService");
            if(downloadService == null)
            {
                StsException.systemError("JNLPUtilities.isResourceCached() failed. Couldn't find downloadService");
                return false;
            }
            return downloadService.isResourceCached(url, null);
        }
        catch(Exception e)
        {
            String urlString;
            if(url != null) urlString = url.toString();
            else            urlString = "nullURL";
            StsException.outputException("JNLPUtilities.isResourceCached() failed for file: " + urlString,
                e, StsException.WARNING);
            return false;
        }
    }

    static public final URL downloadJarURL(String jarname)
    {
        URL url = null;

        DownloadService downloadService;
        try
        {
            downloadService = (DownloadService)getService("javax.jnlp.DownloadService");

            if(downloadService == null)
            {
                StsException.systemError("JNLPUtilities.getWebStartJarURL() failed. Couldn't find downloadService");
                return null;
            }
            URL codebase = getCodeBase();
            if(debug) System.out.println("jar: " + codebase.toString() + jarname);
            url = new URL(codebase, jarname);
            if(!downloadService.isResourceCached(url, null))
            {
                if(debug) System.out.println("Jar file is not cached, downloading...url: " + url.toString());
                DownloadServiceListener downloadListener = downloadService.getDefaultProgressWindow();
                downloadService.loadResource(url, null, downloadListener);
            }
            else
               if(debug) System.out.println("Jar file is cached, no download required...");
            return url;
        }
        catch(Exception e)
        {
            String urlString;
            if(url != null) urlString = url.toString();
            else            urlString = "nullURL";
            StsException.outputException("JNLPUtilities.downloadJarURL() failed." +
                " Couldn't find file in cache and couldn't download it: " + urlString,
                e, StsException.WARNING);
            return null;
        }
    }

    static public final URL downloadJNLP(String jnlp)
    {
        URL url = null;

        BasicService basicService;
        try
        {
            basicService = (BasicService)getService("javax.jnlp.BasicService");
            if(basicService == null)
            {
                StsException.systemError("JNLPUtilities.getWebStartJarURL() failed. Couldn't find basicService");
                return null;
            }
            URL codebase = getCodeBase();
            if(debug) System.out.println("jar: " + codebase.toString() + jnlp);
            url = new URL(codebase, jnlp);
            basicService.showDocument(url);
            return url;
        }
        catch(Exception e)
        {
            String urlString;
            if(url != null) urlString = url.toString();
            else            urlString = "nullURL";
            StsException.outputException("JNLPUtilities.downloadJarURL() failed." +
                " Couldn't find file in cache and couldn't download it: " + urlString,
                e, StsException.WARNING);
            return null;
        }
    }

    public String getInstalledPath()
    {
        ExtensionInstallerService eiService = getInstallService();
        return eiService.getInstallPath();
    }
/*
    static public final String[] getResourceFilenames(ClassLoader classLoader)
    {
        try
        {
            System.out.println("WebStart classloader: " + classLoader.getClass().getName());
            JNLPClassLoader jnlpClassLoader = (JNLPClassLoader)classLoader;
            LaunchDesc launchDesc = jnlpClassLoader.getLaunchDesc();
            ResourcesDesc resourcesDesc = launchDesc.getResources();
            JARDesc[] jarDescs = resourcesDesc.getEagerOrAllJarDescs(true);
            if(jarDescs == null || jarDescs.length == 0)
            {
                System.out.println("Failed to find any resources in JNLP");
                return new String[0];
            }

            String[] filenames = new String[jarDescs.length];
            int nFilenames = 0;
            for(int n = 0; n < jarDescs.length; n++)
            {
                if(jarDescs[n].isJavaFile())
                {
                    URL jarURL = jarDescs[n].getLocation();
                    String pathname = jarURL.getPath();
                    System.out.println("getResourceFilenames() pathname: " + pathname);
                    String filename = getFilenameFromPathname(pathname);
                    filenames[nFilenames++] = filename;
                }
            }
            return (String[])StsMath.trimArray(filenames, nFilenames);
        }
        catch(Exception e)
        {
            System.out.println("JNLPUtilities.getRsourceFilenbames() failed.");
            e.printStackTrace();
            return new String[0];
        }
    }
*/
    static public String getFilenameFromPathname(String filename) throws StsException
    {
        String fileSeparator = File.separator;
        int separatorIndex = filename.lastIndexOf(fileSeparator);
        if(separatorIndex == -1 && !fileSeparator.equals("/"))
            separatorIndex = filename.lastIndexOf("/");
        if(separatorIndex == -1)
            System.out.println("StsFile(pathname) failed. Didn't find separator.");
        int length = filename.length();
        return filename.substring(separatorIndex+1, length);
    }
/*
    public final URL getJarFileEntryURL(String jarFilename, String entryName)
    {
        String urlEntryName;

        try
        {
            classInitialize();
            String codebaseName = codebase.toString();
            urlEntryName = new String("jar:" + codebaseName + jarFilename + "!/" + entryName);
            return new URL(urlEntryName);
        }
        catch(Exception e)
        {
            StsException.outputException("StsProject.getJarFileEntryURL() failed." +
            "Couldn't find jar file entry: " + jarFilename + "/" + entryName,
            e, StsException.WARNING);
            return null;
        }
    }
*/
/*
    static public String getInstallDirectory()
    {
        ExtensionInstallerService eis;
        eis = (ExtensionInstallerService)JNLPUtilities.getService("javax.jnlp.ExtensionInstallerService");
        return eis.getInstallPath();
    }
*/
}
