
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

public class CPS3FaultCutFile extends StsFaultCutFile
{
    public CPS3FaultCutFile(boolean read, String filename) throws StsException
    {
        super(read, filename);
    }

    // read a fault cut from an OpenWorks map file
    public StsFaultCut readFaultCut() throws IOException, StsException
    {
        ArrayList xList = new ArrayList(5);
        ArrayList yList = new ArrayList(5);

        String line;
      	while (((line = readLine()) != null))
        {
            // parse the line
            StringTokenizer stok = new StringTokenizer(line);
            if (stok.countTokens() == 0) break;
            if (stok.countTokens() < 2)
            {
                throw new StsException(StsException.WARNING, "File " + getFilename() + " is corrupt or is not an CPS3 format.");
            }
            xList.add(Double.valueOf(stok.nextToken()));
            yList.add(Double.valueOf(stok.nextToken()));
  		}

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

        for (int n=0; n<nPnts; n++) writeLine(CPS3Line(absXYs[n][0], absXYs[n][1], n+1));
        return true;
    }

    // one line in openworks format
    static private String CPS3Line(double x, double y, int id)
    {
        return x + "   " + y + "   " + "   " + (float)id;
    }

    /** test program */
    public static void main(String[] args)
    {
        final String filename = "e:/Sts/texas_data_hdf/OW_fault_cuts.txt";
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
