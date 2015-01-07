
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Types;

import com.Sts.Utilities.*;

import java.io.*;
import java.util.*;

/** this class encapsulates parameters for a 3D project */
public class StsProjectParms implements Serializable
{
    // constants
    transient static public final String filterString = "Project Files";
    transient static public final String fileFilter = "proj";

    /** input data - all required fields without defaults */
    protected double xOrigin = 519000.0;
    protected double xSize = 23500.0;
    protected double yOrigin = 5980000.0;
    protected double ySize = 16500.0;
    protected float xIncr = 1000.0f;
    protected float xGridIncr = 1000.0f;
    protected float yIncr = 1000.0f;
    protected float yGridIncr = 1000.0f;
    protected float zMin = 5400.0f;
    protected float zMax = 7200.0f;
    protected float zIncr = 50.0f;
    protected float zGridIncr = 25.0f;
    protected float zGrid = 7200.0f;

    /** 'transient' modifier used here to indicate 'optional' parms that have
        reasonable defaults;  but: these are still written out by packFields() */
    transient protected float xyAngle = 0.0f;
    transient protected boolean zDomain = true;
    transient protected float mapGenericNull = 1.0e30f;
	transient protected float logNull = -999.25f;

    transient protected File projectFileDir = new File(".");
    transient protected String name = "Q";

    /** constructTraceAnalyzer using parms */
    public StsProjectParms(double xOrigin, double xSize, float xIncr, float xGridIncr,
                         double yOrigin, double ySize, float yIncr, float yGridIncr,
                         float zMin, float zMax, float zIncr, float zGridIncr,
                         float zGrid)
    	throws StsException
  	{
        /** validate and save inputs; then derive some parms */
        setXParms(xOrigin, xSize, xIncr, xGridIncr);
        setYParms(yOrigin, ySize, yIncr, yGridIncr);
        setZParms(zMin, zMax, zIncr, zGridIncr, zGrid);
  	}

    /** default constructTraceAnalyzer */
    public StsProjectParms() throws StsException
    {
//    	this(-50.0, 100.0, 10.0f, 10.0f, -50.0, 100.0, 10.0f, 10.0f,
//          	 0.0f, 100.0f, 10.0f, 10.0f, 100.0f);
    }
/*
    public static boolean isProjectFile(File f)
    {
        boolean addProfileOk = false;
    	FileReader fileReader = null;
    	try
        {
      		fileReader = new FileReader(f);
      		char[] charBuffer = new char[8192];
      		StringBuffer strBuffer = new StringBuffer();

      		while (fileReader.read(charBuffer) != -1)
            {
        		strBuffer.append(charBuffer);
      		}
            String buf = strBuffer.toString();
            addProfileOk = true;
            String[] fields = new String[]
            {
                "xOrigin",
                "xSize",
                "xIncr",
                "xGridIncr",
                "yOrigin",
                "ySize",
                "yIncr",
                "yGridIncr",
                "zDomain",
                "zMin",
                "zMax",
                "zIncr",
                "zGridIncr",
                "zGrid",
                "xyAngle",
				"logNull"
            };
            for( int i=0; i<fields.length && addProfileOk; i++ )
                if( buf.indexOf(fields[i]) == -1 ) addProfileOk = false;

    	}
    	catch (Exception ex) {	}
    	finally
        {
      		try	{ fileReader.close(); }
      		catch (IOException ex) { }
    	}
        return addProfileOk;
    }
*/
    /** validate input X parms; derive size and no. of divisions */
    private void setXParms(double xOrigin, double xSize, float xIncr, float xGridIncr)
      	throws StsException
    {
		checkAxisParms(xOrigin, xOrigin+xSize, xIncr, xGridIncr);
    	this.xOrigin = xOrigin;
    	this.xSize = xSize;
    	this.xIncr = xIncr;
    	this.xGridIncr = xGridIncr;
  	}

    /** X  gets */
    public double getXOrigin() { return xOrigin; }
    public double getXSize() { return xSize; }
    public double getXMin() { return xOrigin; }
    public double getXMax() { return xOrigin+xSize; }
    public float getXIncrement() { return xIncr; }
    public float getXGridIncrement() { return xGridIncr; }

    /** validate input Y parms; derive size and no. of divisions */
    public void setYParms(double yOrigin, double ySize, float yIncr, float yGridIncr)
      	throws StsException
    {
        checkAxisParms(yOrigin, yOrigin+ySize, yIncr, yGridIncr);
        this.yOrigin = yOrigin;
        this.ySize = ySize;
        this.yIncr = yIncr;
        this.yGridIncr = yGridIncr;
    }

    /** Y  gets */
    public double getYOrigin() { return yOrigin; }
    public double getYSize() { return ySize; }
    public double getYMin() { return yOrigin; }
    public double getYMax() { return yOrigin+ySize; }
    public float getYIncrement() { return yIncr; }
    public float getYGridIncrement() { return yGridIncr; }

    /** validate input Z parms; derive size and no. of divisions */
    public void setZParms(float zMin, float zMax, float zIncr, float zGridIncr, float zGrid) throws StsException
    {
    	checkAxisParms(zMin, zMax, zIncr, zGridIncr);
    	if (zGrid<zMin || zGrid>zMax) // reference Z must be in Z range
        {
      		throw new StsException("zGrid must be within zMin and zMax");
    	}

        this.zMin = zMin;
        this.zMax = zMax;
        this.zIncr = zIncr;
        this.zGridIncr = zGridIncr;
        this.zGrid = zGrid;
    }

    /** Z  gets */
    public float getZMin() { return zMin; }
    public float getZMax() { return zMax; }
    public float getZIncrement() { return zIncr; }
    public float getZGridIncrement() { return zGridIncr; }
    public float getZGridReference() { return zGrid; }

    /** validate and set rotation angle */
    public void setXYRotationAngle(float angle) throws StsException
    {
    	if (angle < 0.0 || angle > 360.0)
        {
      		throw new StsException("StsProjectParms.setXYRotationAngle: "
                    			   + " Angle must be within 0 to 360");
    	}

    	xyAngle = angle;
    }

    /** get rotation angle */
    public float getXYRotationAngle() { return xyAngle; }

    /** set/get z domain flag */
    public void setZDomain(boolean zDomain) { this.zDomain = zDomain; }
    public boolean getZDomain() { return zDomain; }

    /** set/get null value for generic maps */
    public void setMapGenericNull(float nullValue) { mapGenericNull = nullValue; }
    public float getMapGenericNull() { return mapGenericNull; }

	public void setLogNull(float logNull) { this.logNull = logNull; }
	public float getLogNull() { return logNull; }

    /** get project file directory */
    public File getProjectFileDirectory() { return projectFileDir; }

    /** set/get project name */
    public void setName(String name) { this.name = name; }
    public String getName() { return name; }

    /** validate axis range and increment */
    private void checkAxisParms(double min, double max, float incr, float gridIncr) throws StsException
    {
		if (min >= max)  // max must exceed min
        {
      		throw new StsException(StsException.WARNING, "StsProjectParms.checkAxisParms: "
                  				   + "Max must be larger than Min");
    	}
    	if (incr <= 0.0)  // non-positive increment
        {
      		throw new StsException(StsException.WARNING, "StsProjectParms.checkAxisParms: "
            					   + "Increment must be larger than 0");
    	}
    	if (gridIncr <= 0.0)  // non-positive increment
        {
      		throw new StsException(StsException.WARNING, "StsProjectParms.checkAxisParms: "
                  				   + "Grid increment must be larger than 0");
    	}
    	if (incr > max-min)
        {
      		throw new StsException(StsException.WARNING, "StsProjectParms.checkAxisParms: "
                  				   + "Increment cannot exceed size of axis");
    	}
    	if (gridIncr > max-min)
        {
      		throw new StsException(StsException.WARNING, "StsProjectParms.checkAxisParms: "
                  				   + "Increment cannot exceed size of axis");
    	}
  	}

    /** build the parms object */
    static public StsProjectParms loadFile(File file)
    {
        StsProjectParms parms = null;
        String filePath = null;

      	try
        {
            if(file == null) return null;
            parms = new StsProjectParms();
            filePath = file.getPath();
        	StsParameterFile.initialReadObjectFields(filePath, parms, null, null);
            parms.projectFileDir = new File(file.getParent());
            parms.buildName(file.getName());
      	}
      	catch (StsException e)
        {
            StsException.outputException("StsProjectParms.loadFile() failed for file: " + filePath, e, StsException.WARNING);
        	return null;
      	}
        return parms;
    }

    static public StsProjectParms loadFile()
    {
        StsProjectParms parms = null;
      	try
        {
            parms = new StsProjectParms();
            String filename = StsParameterFile.askReadPackedFields(parms, fileFilter, filterString);
        	if (filename==null) return null;
            if (parms!=null)
            {
                // get project file's location
                File f = new File(filename);
                parms.projectFileDir = new File(f.getParent());
            }
            parms.buildName(filename);
      	}
      	catch (StsException ex)
        {
        	return null;
      	}
        return parms;
    }

    /** use the selected file, then try to open it and build a parms object */
    static public StsProjectParms loadFile(String filename)
    {
        StsProjectParms parms = null;
      	try
        {
            parms = new StsProjectParms();
        	StsParameterFile.initialReadObjectFields(filename, parms, null, null);
            if (parms!=null)
            {
                // get project file's location
                File f = new File(filename);
                parms.projectFileDir = new File(f.getParent());
            }
            parms.buildName(filename);
      	}
      	catch (StsException ex)
        {
        	return null;
      	}
        return parms;
    }

    private void buildName(String filename)
    {
        if (name==null)  // build project name
        {
            StringTokenizer stok = new StringTokenizer(filename, ".");
            int nTokens = stok.countTokens();
            for (int i=0; i<nTokens-1; i++) stok.nextToken();
            name = new String(stok.nextToken());
            if (name==null) name = new String("Default");
        }
    }

    /* test main */
  	public static void main(String[] args)
  	{
//            File f = new File("data" + File.separator + "logs" + File.separator + "junk");
            File f = new File(".." + File.separator + "Folder");
            System.out.println("relative name = " + f.toString() + "junk");
            System.out.println("File = " + f.toString());
            System.out.println("Path = " + f.getPath());
            System.out.println("Absolute path = " + f.getAbsolutePath());
            try { System.out.println("Canonical path = " + f.getCanonicalPath()); }
            catch (IOException e) { }
            System.out.println("Name = " + f.getName());
            System.out.println("Parent = " + f.getParent());
/*
   	try
        {
            StsProjectParms projParmsOut = new StsProjectParms();
            System.out.println("projParmsOut object:");
            System.out.println(projParmsOut.packFields());

            if (!StsFieldIO.writePackedFields(projParmsOut))
            {
            	System.out.println("writePackFields cancelled");
            }

            //StsFieldIO.writeSerializedFields(projParmsOut);
            StsProjectParms projParmsIn;
            //projParmsIn = (StsProjectParms)StsFieldIO.readSerializedFields();
            projParmsIn = new StsProjectParms(0f, 2f, 1f, 1f, 0f, 2f, 1f, 1f, 0f, 2f, 1f, 1f, 0f);

            System.out.println("projParmsIn object:");
            System.out.println(projParmsIn.packFields());

            if (StsFieldIO.askReadPackedFields(projParmsIn, null)!=null)
            {
            	System.out.println("projParmsIn object after readPackedFields():");
            	System.out.println(projParmsIn.packFields());
            }
            else
            {
            	System.out.println("readPackFields cancelled");
            }
        }
        catch (StsException e)
        {
            System.out.println("StsProjectParms main() -- StsException: " + e.toString());
            return;
        }
*/
  	}
}

