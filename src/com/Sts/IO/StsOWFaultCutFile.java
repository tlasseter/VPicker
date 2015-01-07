
//Title:        S2S: Seismic-to-simulation: Seismic-to-simulation
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

public class StsOWFaultCutFile extends StsFaultCutFile
{
    static final protected String NL = System.getProperty("line.separator");

    static final private int OW_FIRST_POINT = 6;
    static final private int OW_INTERIOR_POINT = 7;
    static final private int OW_LAST_POINT = 8;

    public StsOWFaultCutFile(boolean read, String filename) throws StsException
    {
        super(read, filename);
    }

    // read a fault cut from an OpenWorks map file
    public StsFaultCut readFaultCut() throws IOException, StsException
    {
        ArrayList xList = new ArrayList(5);
        ArrayList yList = new ArrayList(5);

        String line;
        boolean gotFirst = false;
        boolean gotLast = false;
      	while (!gotLast && ((line = readLine()) != null))
        {
            // parse the line
            StringTokenizer stok = new StringTokenizer(line);
            if (stok.countTokens() < 4)
            {
                throw new StsException(StsException.WARNING, "File " + getFilename() + " is corrupt or is not an OpenWorks format.");
            }
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
        return new StsFaultCut(x, y);
    }

    // write a fault cut to an OpenWorks map file
    public boolean writeFaultCut(StsFaultCut faultCut) throws IOException
    {
        if (faultCut == null) return false;
        double[][] absXYs = faultCut.getAbsoluteXYs();
        if (absXYs == null) return false;
        int nPnts = absXYs.length;
        if (nPnts < 2) return false;

        writeLine(OWLine(absXYs[0][0], absXYs[0][1], OW_FIRST_POINT));
        for (int i=1; i<nPnts-1; i++)
            writeLine(OWLine(absXYs[i][0], absXYs[i][1], OW_INTERIOR_POINT));
        writeLine(OWLine(absXYs[nPnts-1][0], absXYs[nPnts-1][1], OW_LAST_POINT));
        return true;
    }

    // one line in openworks format
    static private String OWLine(double x, double y, int owPointType)
    {
        return x + "\t" + y + "\t" + com.Sts.Utilities.StsParameters.nullValue
                + "\t" + owPointType + NL;
    }

    /** test program */
    public static void main(String[] args)
    {
        final String filename = "e:/Sts/texas_data/OW_fault_cuts.txt";
        System.out.println("fault cut file: " + filename);
        try
        {
            StsOWFaultCutFile f = new StsOWFaultCutFile(true, filename);
            StsFaultCut[] fcs = f.getFaultCuts();
            if (fcs == null) System.out.println("no fault cuts read");
            else System.out.println("read " + fcs.length + " fault cuts");
            for (int i=0; i<fcs.length; i++) fcs[i].print();
        }
        catch (StsException e) { System.err.println(e.toString()); }
        System.exit(0);
    }

}
