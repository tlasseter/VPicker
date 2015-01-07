
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

/** an ascii file holding fault polygons */
public class StsFaultPolygonFile
{
    // constants
    final static public int OPEN_WORKS = 0;

    final static private int OW_FIRST_POINT = 6;
    final static private int OW_INTERIOR_POINT = 7;
    final static private int OW_LAST_POINT = 8;

    private int type;
    private File file;
    private StsFaultPolygon[] faultPolygons;
    private BufferedReader bufferedReader;

    public StsFaultPolygonFile(String filename, int type) throws StsException
    {
        this.type = type;
        if (!typeIsValid()) throw new StsException(StsException.FATAL,
                "StsFaultPolygonFile:  unknown fault polygon file type.");

        file = new File(filename);
        if (!file.exists()) throw new StsException(StsException.FATAL,
                "StsFaultPolygonFile:  fault polygon file doesn't exist.");
        if (!file.canRead()) throw new StsException(StsException.FATAL,
                "StsFaultPolygonFile:  fault polygon file can't be read.");
    }

    // check for valid type
    private boolean typeIsValid()
    {
        switch (type)
        {
            case OPEN_WORKS: return true;
        }
        return false;
    }

    /** accessors */
    public String getPath() { return file.getParent(); }
    public String getFilename() { return file.getName(); }

    /** read and return fault polygons */
    public StsFaultPolygon[] getFaultPolygons()
    {
        if (faultPolygons == null && !readFaultPolygons()) return null;
        return faultPolygons;
    }

    // read fault polygons from the file
    private boolean readFaultPolygons()
    {
     	FileReader fileReader = null;
        ArrayList polygonList = new ArrayList(5);
        String path = getPath();
        StsFaultPolygon fp;
    	try
        {
      		fileReader = new FileReader(file.getPath());
            bufferedReader = new BufferedReader(fileReader);
            while ((fp = readFaultPolygon()) != null)
            {
                polygonList.add(fp);
            }
    	}
    	catch (Exception e) { return false; }
    	finally
        {
      		try
        	{
        		fileReader.close();
                int n = polygonList.size();
                faultPolygons = new StsFaultPolygon[n];
                for (int i=0; i<n; i++)
                {
                    faultPolygons[i] = (StsFaultPolygon)polygonList.get(i);
                }
      		}
      		catch (Exception e) { return false; }
    	}
        return true;
  	}

    // read a fault polygon from the file
    private StsFaultPolygon readFaultPolygon()
    {
        try
        {
            switch (type)
            {
                case OPEN_WORKS: return readOpenWorksFaultPolygon();
            }
        }
        catch (IOException e) { return null; }
        return null;
    }

    // read a fault polygon from an OpenWorks map file
    private StsFaultPolygon readOpenWorksFaultPolygon() throws IOException
    {
        ArrayList xList = new ArrayList(5);
        ArrayList yList = new ArrayList(5);

        String line;
        boolean gotFirst = false;
        boolean gotLast = false;
      	while (!gotLast && ((line = bufferedReader.readLine()) != null))
        {
            // parse the line
            StringTokenizer stok = new StringTokenizer(line);
            if (stok.countTokens() != 5) return null;
            xList.add(Double.valueOf(stok.nextToken()));
            yList.add(Double.valueOf(stok.nextToken()));
            stok.nextToken();  // Z not used
            int ptype = Integer.valueOf(stok.nextToken()).intValue();
            switch (ptype)
            {
                case OW_FIRST_POINT:
                    if (gotFirst) return null;
                    gotFirst = true;
                    break;
                case OW_INTERIOR_POINT:
                    if (!gotFirst) return null;
                    break;
                case OW_LAST_POINT:
                    if (!gotFirst) return null;
                    gotLast = true;
                    break;
                default: return null;
            }
  		}
        if (!gotLast) return null;

        int n = xList.size();
        if (n < 1) return null;
        double[] x = new double[n];
        double[] y = new double[n];
        for (int i=0; i<n; i++)
        {
            Double xObj = (Double)xList.get(i);
            x[i] = xObj.doubleValue();
            Double yObj = (Double)yList.get(i);
            y[i] = yObj.doubleValue();
        }
        return new StsFaultPolygon(x, y);
    }

    /** test program */
    public static void main(String[] args)
    {
        //final String filename = "e:/Sts/Data/Fault Polygons/FaultPolygons.txt";
        final String filename = "e:/Sts/Data/Fault Polygons/FaultPolygonsOut.txt";
        System.out.println("Polygon file: " + filename);
        try
        {
            StsFaultPolygonFile fpf = new StsFaultPolygonFile(filename, OPEN_WORKS);
            StsFaultPolygon[] fps = fpf.getFaultPolygons();
            if (fps == null) System.out.println("no fault polygons read");
            else System.out.println("read " + fps.length + " fault polygons");
        }
        catch (StsException e) { System.err.println(e.toString()); }
    }

}
