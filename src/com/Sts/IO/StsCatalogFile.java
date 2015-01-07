
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.IO;

import com.Sts.Types.*;
import com.Sts.Utilities.*;

import java.io.*;

public class StsCatalogFile
{
    File file;
    String[] allFilenames = null;
    String[] validFilenames = null;

    public StsCatalogFile(String filename) throws StsException
    {
        file = new File(filename);
        if (!file.exists()) throw new StsException(StsException.FATAL,
                "StsCatalogFile.StsCatalogFile:  catalog file doesn't exist.");
        if (!file.canRead()) throw new StsException(StsException.FATAL,
                "StsCatalogFile.StsCatalogFile:  catalog file can't be read.");
    }

    public String getPath() { return file.getParent(); }

    public String[] getFilenames()
    {
        if (allFilenames == null)
        {
            if (!readFilenames()) return null;
        }
        return allFilenames;
    }

    public String[] getValidFilenames()
    {
        if (validFilenames == null)
        {
            if (!readFilenames()) return null;
        }
        return validFilenames;
    }

    private boolean readFilenames()
    {
     	FileReader fileReader = null;
        StsList allFiles = new StsList(1, 1);
        StsList validFiles = new StsList(1, 1);
        String path = getPath();
    	try
        {
      		fileReader = new FileReader(file.getPath());
            BufferedReader bufferedReader = new BufferedReader(fileReader);
            String line = bufferedReader.readLine();
      		while (line != null)
            {
                allFiles.add(line);
                File f = new File(path, line);
                if (f.exists() && f.canRead()) validFiles.add(line);
                line = bufferedReader.readLine();
      		}
    	}
    	catch (Exception e) { return false; }
    	finally
        {
      		try
        	{
        		fileReader.close();
                int nFiles = allFiles.getSize();
                int nValidFiles = validFiles.getSize();
                if (nFiles > 0)
                {
                    allFilenames = new String[nFiles];
                    for (int i=0; i<nFiles; i++)
                    {
                        String filename = (String)allFiles.getElement(i);
                        allFilenames[i] = new String(filename);
                    }
                }
                if (nValidFiles > 0)
                {
                    validFilenames = new String[nValidFiles];
                    for (int i=0; i<nValidFiles; i++)
                    {
                        String filename = (String)validFiles.getElement(i);
                        validFilenames[i] = new String(filename);
                    }
                }
      		}
      		catch (Exception e) { return false; }
    	}
        return true;
  	}
}
