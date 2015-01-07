
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
import java.util.*;

/** an ascii file holding fault cut lines */
abstract public class StsFaultCutFile
{
    private boolean read;
    private File file;
    private StsFaultCut[] faultCuts;
    private BufferedReader bufferedReader;
    private BufferedWriter bufferedWriter;

    /** constructor */
    public StsFaultCutFile(boolean read, String filename) throws StsException
    {
        this.read = read;
        file = new File(filename);
        if (read)
        {
            if (!file.exists()) throw new StsException(StsException.FATAL,
                    "StsFaultCutFile:  can't find fault cut file: " + filename);
            if (!file.canRead()) throw new StsException(StsException.FATAL,
                    "StsFaultCutFile:  can't read fault cut file: " + filename);
        }
        else  // write
        {
            if (file.exists() && !file.delete()) throw new StsException(StsException.FATAL,
                    "StsFaultCutFile:  unable to delete fault cut file: " + filename);
        }
    }


    /** accessors */
    public String getPath() { return file.getParent(); }
    public String getFilename() { return file.getName(); }

    public String readLine() throws java.io.IOException { return bufferedReader.readLine(); }
    public void writeLine(String line) throws java.io.IOException { bufferedWriter.write(line); }

    /** read and return fault cuts */
    public StsFaultCut[] getFaultCuts() throws StsException
    {
        if (faultCuts == null && !readFaultCuts()) return null;
        return faultCuts;
    }

    /** set and write fault cuts */
    public boolean setFaultCuts(StsFaultCut[] faultCuts)
    {
        if (faultCuts == null) return false;
        try
        {
            open();
            for (int i=0; i<faultCuts.length; i++)
            {
                if (!writeFaultCut(faultCuts[i])) return false;
            }
            close();
        }
        catch (IOException e) { return false; }
        return true;
    }

    // buffered writer methods
    public void open() throws IOException
    {
        bufferedWriter = new BufferedWriter(new FileWriter(file));
    }
    public void close() throws IOException
    {
        bufferedWriter.close();
    }

    // read fault cuts from the file
    private boolean readFaultCuts() throws StsException
    {
     	FileReader fileReader = null;
        ArrayList faultCutList = new ArrayList(5);
        String path = getPath();
        StsFaultCut fp;
    	try
        {
      		fileReader = new FileReader(file.getPath());
            bufferedReader = new BufferedReader(fileReader);
            while ((fp = readFaultCut()) != null)
            {
                faultCutList.add(fp);
            }
    	}
    	catch (StsException Stse)
        {
            throw Stse;
        }
    	catch (Exception e) { return false; }
    	finally
        {
      		try
        	{
        		fileReader.close();
                int n = faultCutList.size();
                faultCuts = new StsFaultCut[n];
                for (int i=0; i<n; i++)
                {
                    faultCuts[i] = (StsFaultCut)faultCutList.get(i);
                }
      		}
      		catch (Exception e) { return false; }
    	}
        return true;
  	}

    // read a fault cut from the file
    abstract public StsFaultCut readFaultCut() throws IOException, StsException;

    // write a fault cut
    abstract public boolean writeFaultCut(StsFaultCut faultCut) throws IOException;
}
