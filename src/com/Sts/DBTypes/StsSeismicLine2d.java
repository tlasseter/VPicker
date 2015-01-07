package com.Sts.DBTypes;

/**
 * <p>Title: Workflow development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2001</p>
 * <p>Company: 4D Systems LLC</p>
 * @author unascribed
 * @version 1.0
 */

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.IO.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.io.*;
import java.nio.*;
import java.util.*;

/**
 * This class is currently being used for poststack 2d lines which are added to an existing StsSeismicVolume.
 * Will no doubt be replaced/refactored with a design similar to prestack 2d and 3d lines and lineSets.
 * cdpX and cdpY have been changed from absolute doubles to relative floats to keep from converting all the time to draw,
 * but has not been tested yet.  TJL 1/18/07
 */
public class StsSeismicLine2d extends StsSeismicLine implements StsTextureSurfaceFace //, StsEfficientRepaintable
{
    /** line set this line belongs to */
    protected StsSeismicLineSet seismicLineSet;
    // these members are persistent, but not loaded from seis3d.txt.name file
    // protected StsColorscale colorscale;

    public StsFilterProperties filterProperties = null;
    public StsAGCProperties agcProperties = null;

    protected boolean displayAxis = true;
    protected boolean displayVAR = true;
    protected boolean displayWiggles = false;
    protected boolean contourColors = true;
    protected int wiggleToPixelRatio = 1;
    protected float zShift = 0.0f;

    transient boolean readoutEnabled = false;
    transient StsFile out = null;
    transient OutputStream os = null;
    transient BufferedOutputStream bos = null;
    transient ByteArrayOutputStream baos = null;
    transient DataOutputStream ds = null;
    transient byte exportType = StsSeismicExportPanel.BYTE;
    transient byte exportNullType = StsSeismicExportPanel.NULL_NONE;
    transient float exportScale = 1.0f;
    transient String exportTextHeader = null;

    transient StsSeismicLineSetClass lineSetClass;
    /** Display lists should be used (controlled by View:Display Options) */
    transient boolean useDisplayLists = true;
    /** Display lists currently being used for surface geometry */
    transient boolean usingDisplayLists = true;

    // the following are initialized after reading parameters file or reloading database
    transient protected StsSpectrumDialog spectrumDialog = null;
    transient protected boolean spectrumDisplayed = true;

    transient protected int colorListNum = 0;
    transient protected boolean colorListChanged = true;

    transient int sampleSize = 1;
    transient boolean attributeInMdepth = false;
    transient byte attributeDomain = StsParameters.TD_DEPTH;

    transient Vector itemListeners = null;
    transient StsActionListeners actionListeners = null;

    transient StsTimer timer = null;
    transient boolean runTimer = false;

    transient int displayListNum = 0;

    transient float[] floatData = null; // data read for this vsp
    transient float[] scaledFloatData = null; // scaled, agc'd etc trace data
    //transient byte[] byteData = null; // byte data converted from fdata

    transient StsTextureTiles textureTiles = null;
    transient private boolean textureChanged = true;
    transient float azimuth;
    transient int panelWidth = 100;
    transient float[][] panelRanges;
    transient boolean isDisplayWiggles = false;
    transient boolean isPixelMode = false;
    transient String exportDirectory = null;
    transient StsSeismicLine2dClass vspClass;

    static protected StsObjectPanel objectPanel = null;

    static public final byte nullByte = StsParameters.nullByte;
    static public final int nullUnsignedInt = 255;

    static public final String group = StsSeismicBoundingBox.group2d;

    static final boolean debug = false;
    static final boolean voxelDebug = false;
    static final boolean wiggleDebug = false;

//	static public final StsFloatFieldBean tpiBean = new StsFloatFieldBean(StsSeismicLine2d.class, "tracesPerInch", 1, 250, "Traces per Inch:");

    public String getGroupname()
    {
        return group2d;
    }

    public StsSeismicLine2d()
    {
    }

    public StsSeismicLine2d(boolean persistent)
    {
        super(persistent);
    }

    private StsSeismicLine2d(String directory, String filename, StsModel model, StsSeismicLineSet seismicLineSet) throws FileNotFoundException, StsException
    {
        super(false);
        String pathname = directory + filename;
        this.seismicLineSet = seismicLineSet;
        if(!(new File(pathname)).exists())
        {
            throw new FileNotFoundException();
        }
        StsParameterFile.initialReadObjectFields(pathname, this, StsSegyLine2d.class, StsMainObject.class);
        setName(getStemname());

        stsDirectory = directory;
        if(!initialize(model))
        {
            throw new FileNotFoundException(pathname);
        }
        isVisible = true;
    }

    static public StsSeismicLine2d constructor(StsFile file, StsModel model, StsSeismicLineSet lineSet)
    {
        return StsSeismicLine2d.constructor(file.getDirectory(), file.getFilename(), model, lineSet);
    }

    static public StsSeismicLine2d constructor(String directory, String filename, StsModel model, StsSeismicLineSet lineSet)
    {

        StsSeismicLine2d line = null;
        try
        {
            line = new StsSeismicLine2d(directory, filename, model, lineSet);
            if(line != null) model.setCurrentObject(line);
//			volume.refreshObjectPanel();
            return line;
        }
        catch(FileNotFoundException e)
        {
            return null;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSeismicLine2d.constructor() failed.", e, StsException.WARNING);
            return null;
        }
    }

    public boolean initialize(StsFile file, StsModel model)
    {
        try
        {
            String pathname = file.getDirectory() + file.getFilename();
            StsParameterFile.initialReadObjectFields(pathname, this, StsSegyLine2d.class, StsBoundingBox.class);
            setName(getStemname());

            stsDirectory = file.getDirectory();
            if(!initialize(model))
            {
                return false;
            }
            isVisible = true;


            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSeismicLine2d.loadFile() failed.",
                    e, StsException.WARNING);
            return false;
        }
    }

    public boolean initialize(StsModel model)
    {
        try
        {
            super.initialize(currentModel);
            lineSetClass = seismicLineSet.lineSetClass;
            if(seismicLineSet != null) seismicLineSet.addActionListener(this);
        }
        catch(Exception e)
        {
            StsException.outputException("StsSeismicLine2d.classInitialize(model) failed.", e, StsException.WARNING);
//			StsMessage.printMessage("Failed to find file. Error: " + e.getMessage());
            return false;
        }
        return true;
    }

    private boolean initialize()
    {
        super.initialize(currentModel);

        byte zDomainByte = StsParameters.getZDomainFromString(zDomain);
        zDomain = StsParameters.TD_ALL_STRINGS[zDomainByte];
        return currentModel.getProject().checkSetZDomain(zDomainByte, zDomainByte);
    }
/*
	public void displayWiggleProperties()
	{
		getWiggleDisplayProperties().displayWiggleProperties("Edit " + getName() + " Wiggle Properties");
	}

	public void displayAGCProperties()
	{
		displayWiggleProperties();
	}
*/

    public String getGroup()
    {
        return group;
    }
    public String createFloatDataFilename()
    {
        return createFilename(floatFormat);
    }

    public String createByteDataFilename()
    {
        return createFilename(byteFormat);
    }

    public float getZShift()
    {
        return zShift;
    }

    public void setZShift(float zshift)
    {
        zShift = zshift;
        textureChanged(true);
        dbFieldChanged("zShift", zShift);
    }

    /**
     * Given a line into the screen at the mousePoint clipped to the project box,
     * find its XY intersection with points along the 2d line.
     * Compute the intersection point (XYZ) and interpolation factor along line (F)
     * and returned the sequence index of the nearest gather.
     */

    public int getMouseSeismicLine2dIntersect(double[][] mouseLine, double[] intersectionXYZF)
    {
        double mouseLineLenSq;
        double[] lineVector;
        double[] pointVector0, pointVector1;
        double x0, y0, x1, y1;

        if(nCols < 1)
        {
            intersectionXYZF[3] = StsParameters.largeDouble;
            return -1;
        }
        try
        {
            mouseLineLenSq = StsMath.distanceSq(mouseLine[0], mouseLine[1], 2);
            lineVector = StsMath.subtract(mouseLine[1], mouseLine[0]);
            pointVector1 = StsMath.vector2(mouseLine[0][0], mouseLine[0][1], cdpX[0], cdpY[0]);
            y1 = computePointCoordinatesY(pointVector1, lineVector, mouseLineLenSq);
            int nPointNearest = 0;
            double fLineNearest = StsParameters.largeDouble;
            for(int n = 1; n < nCols; n++)
            {
                pointVector0 = pointVector1;
                pointVector1 = StsMath.vector2(mouseLine[0][0], mouseLine[0][1], cdpX[n], cdpY[n]);
                y0 = y1;
                y1 = computePointCoordinatesY(pointVector1, lineVector, mouseLineLenSq);

                if(y1 * y0 < 0.0)
                {
                    x0 = computePointCoordinatesX(pointVector0, lineVector, mouseLineLenSq);
                    x1 = computePointCoordinatesX(pointVector1, lineVector, mouseLineLenSq);

                    double fPoint = -y0 / (y1 - y0);
                    double fLine = x0 + fPoint * (x1 - x0);
                    double z = mouseLine[0][2] + fLine * (mouseLine[1][2] - mouseLine[0][2]);
                    if(fLine < fLineNearest && StsMath.betweenInclusive(z, zMin, zMax))
                    {
                        fLineNearest = fLine;

                        if(fPoint <= 0.5)
                            nPointNearest = n - 1;
                        else
                            nPointNearest = n;
                    }
                }
            }
            if(!StsMath.interpolate(mouseLine[0], mouseLine[1], fLineNearest, 3, intersectionXYZF))
                return -1;
            intersectionXYZF[3] = fLineNearest;
            return nPointNearest;
        }
        catch(Exception e)
        {
            return -1;
        }
    }

    private double computePointCoordinatesX(double[] pointVector, double[] lineVector, double lineLenSq)
    {
        return StsMath.dot(lineVector, pointVector) / lineLenSq;
    }

    private double computePointCoordinatesY(double[] pointVector, double[] lineVector, double lineLenSq)
    {
        return StsMath.cross2(lineVector, pointVector) / lineLenSq;
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        byte projectZDomain = StsObject.getCurrentModel().getProject().getZDomain();
        byte volumeZDomain = getZDomain();

        if(projectZDomain == StsParameters.TD_TIME)
        {
            if(volumeZDomain == StsParameters.TD_TIME)
                drawTextureTileTimeSurface(tile, gl);
            else
            {
                StsSeismicVelocityModel velocityVolume = currentModel.getProject().velocityModel;
                if(velocityVolume == null) return;
                drawTextureTileDepthSurfaceInTime(velocityVolume, tile, gl);
            }
        }
        else if(projectZDomain == StsParameters.TD_DEPTH)
        {
            // volumeZDomain is TD_TIME

            if(volumeZDomain == StsParameters.TD_DEPTH) // seismic already in depth, don't need to convert so draw as if in time
                drawTextureTileTimeSurface(tile, gl);
            else
            {
                StsSeismicVelocityModel velocityVolume = currentModel.getProject().velocityModel;
                if(velocityVolume == null) return;
                drawTextureTileTimeSurfaceInDepth(velocityVolume, tile, gl);
            }
        }
    }

    public void drawTextureTileDepthSurface(StsPreStackVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        System.out.println("2D seismic depth display not implemented yet");
        /*
           float cursorXInc, cursorYInc;
           if (dirNo == StsCursor3d.ZDIR)
           {
            return;
           }

           gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
           double rowTexCoor = tile.minRowTexCoor;
           double dRowTexCoor = tile.dRowTexCoor;
           double dColTexCoor = tile.dColTexCoor;
           double[] xyz = tile.xyzPlane[0];
           double x1 = xyz[0];
           double y1 = xyz[1];
           double t1 = xyz[2];
           int volumeRow = this.getNearestBoundedRowCoor( (float)y1);
           int volumeCol = getNearestBoundedColCoor( (float)x1);

  //		StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
  //		if (velocityVolume == null)return;
  //		float depthMin = velocityModel.depthDatum;
           float depthMin = 0.0f;
           int volumeRowInc = 0;
           int volumeColInc = 0;
           if (dirNo == XDIR)
           {
            cursorXInc = 0;
            cursorYInc = getYInc();
            volumeRowInc = 1;
           }
           else // dirNo == StsCursor3d.YDIR
           {
            cursorXInc = getXInc();
            cursorYInc = 0;
            volumeColInc = 1;
           }
           double tInc = getZInc();
           for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
           {
            double x0 = x1;
            double y0 = y1;
            x1 += cursorXInc;
            y1 += cursorYInc;

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            double t = t1 + tile.colMin * tInc;

            for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
            {
             float v0 = velocityVolume.getValue(volumeRow, volumeCol, col);
             float z0 = (float) (v0 * t + depthMin);
             gl.glTexCoord2d(colTexCoor, rowTexCoor);
             gl.glVertex3d(x0, y0, z0);
             float v1 = velocityVolume.getValue(volumeRow + volumeRowInc, volumeCol + volumeColInc, col);
             float z1 = (float) (v1 * t + depthMin);
             gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
             gl.glVertex3d(x1, y1, z1);
            }
            gl.glEnd();
            volumeRow += volumeRowInc;
            volumeCol += volumeColInc;
           }
           */
    }

    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if(textureTiles == null) return;
        textureTiles.deleteTextures(gl);
        textureTiles.deleteDisplayLists(gl);
        textureChanged = true;
    }
/*
    public void addTextureToDeleteList()
    {
        if (textureTiles != null)
        {
            StsTextureList.addTextureToDeleteList(this);
        }
	    textureChanged = true;
    }
*/
    private void deleteDisplayLists(GL gl)
    {
        if(textureTiles != null)
            textureTiles.deleteDisplayLists(gl);
    }

    public boolean textureChanged()
    {
	    textureChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void geometryChanged()
    {
        if(textureTiles == null) return;
//        StsTextureList.addTextureToDeleteList(this);
        textureChanged = true;
    }

    private void drawTextureTileTimeSurface(StsTextureTile tile, GL gl)
    {
        StsGridCrossingPoint gridCrossingPoint;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
        double dRowTexCoor = tile.dRowTexCoor;
        float tileZMin = zMin - zShift + tile.colMin * zInc;
        float tileZMax = zMin - zShift + tile.colMax * zInc;
        gl.glBegin(GL.GL_QUAD_STRIP);
        for(int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float x = cdpX[row];
            float y = cdpY[row];
            gl.glTexCoord2d(tile.minColTexCoor, rowTexCoor);
            gl.glVertex3f(x, y, tileZMin);
            gl.glTexCoord2d(tile.maxColTexCoor, rowTexCoor);
            gl.glVertex3f(x, y, tileZMax);
        }
        gl.glEnd();
    }

    public void drawTextureTileDepthSurfaceInTime(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        //            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
        //		double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        float x1 = cdpX[tile.rowMin];
        float y1 = cdpY[tile.rowMin];
        float tileZMin = zMin - zShift + tile.colMin * zInc;
        float[] times1 = velocityModel.getTimesFromIntervalVelocity(x1, y1, tileZMin, zInc, tile.nCols);
        for(int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float x0 = x1;
            float y0 = y1;
            float[] times0 = times1;
            x1 = cdpX[row];
            y1 = cdpY[row];
            times1 = velocityModel.getTimesFromIntervalVelocity(x1, y1, tileZMin, zInc, tile.nCols);

            gl.glBegin(GL.GL_QUAD_STRIP);
            double colTexCoor = tile.minColTexCoor;
            float z = tileZMin;
            for(int n = 0, col = tile.colMin; col <= tile.colMax; col++, n++, z += zInc, colTexCoor += dColTexCoor)
            {
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex3f(x0, y0, times0[n]);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex3f(x1, y1, times1[n]);
            }
            gl.glEnd();
        }
    }

    public void drawTextureTileTimeSurfaceInDepth(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        if (velocityModel == null)return;

        if(velocityModel.getInputVelocityVolume() != null)
            drawTextureTileDepthSurfaceFromVolume(velocityModel, tile, gl);
        else
            drawTextureTileDepthSurfaceFromIntervalVelocities(velocityModel, tile, gl);
    }

    public void drawTextureTileDepthSurfaceFromVolume(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
        float x1 = cdpX[tile.rowMin];
        float y1 = cdpY[tile.rowMin];
        StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
		int r1 = velocityVolume.getNearestBoundedRowCoor(y1);
		int c1 = velocityVolume.getNearestBoundedColCoor(x1);
        float velocityTInc = velocityVolume.zInc;
        float[] velocityTrace1 = velocityVolume.getTraceValues(r1, c1);
        float depthDatum = velocityModel.depthDatum;
        float timeDatum = velocityModel.timeDatum;
		for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
			float x0 = x1;
			float y0 = y1;
			float[] velocityTrace0 = velocityTrace1;
            x1 = cdpX[row];
            y1 = cdpY[row];
			r1 = velocityVolume.getNearestBoundedRowCoor(y1);
			c1 = velocityVolume.getNearestBoundedColCoor(x1);
			velocityTrace1 = velocityVolume.getTraceValues(r1, c1);

			gl.glBegin(GL.GL_QUAD_STRIP);

			double colTexCoor = tile.minColTexCoor;
			float t = zMin + tile.colMin*zInc;

            for (int col = tile.colMin; col <= tile.colMax; col++, t += zInc, colTexCoor += dColTexCoor)
            {
                float v0 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace0);
                float z0 = (v0 * (t - timeDatum) + depthDatum);
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex3f(x0, y0, z0);
                float v1 = StsMath.interpolateValue(t, timeDatum, velocityTInc, velocityTrace1);
                float z1 = (v1 * (t - timeDatum) + depthDatum);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex3f(x1, y1, z1);
            }
		   gl.glEnd();
		}
	}

    public void drawTextureTileDepthSurfaceFromIntervalVelocities(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
        float x1 = cdpX[tile.rowMin];
        float y1 = cdpY[tile.rowMin];
        float displayTMin = zMin + tile.colMin*zInc;
        int nSlices = tile.colMax - tile.colMin + 1;
        float[] depths1 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, zInc, nSlices, x1, y1);
        for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
			float x0 = x1;
			float y0 = y1;
            x1 = cdpX[row];
            y1 = cdpY[row];
            float[] depths0 = depths1;
            depths1 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, zInc, nSlices, x1, y1);

			gl.glBegin(GL.GL_QUAD_STRIP);

			double colTexCoor = tile.minColTexCoor;

            for (int n = 0, col = tile.colMin; col <= tile.colMax; n++, col++, colTexCoor += dColTexCoor)
            {
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex3f(x0, y0, depths0[n]);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex3f(x1, y1, depths1[n]);
            }
		    gl.glEnd();
		}
	}
/*
    public void drawTextureTileTimeSurfaceInDepth(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
        //            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
        //		double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        float x1 = cdpX[tile.rowMin];
        float y1 = cdpY[tile.rowMin];
        StsSeismicVolume velocityVolume = velocityModel.getVelocityVolume();
        int r1 = velocityVolume.getNearestBoundedRowCoor(y1);
        int c1 = velocityVolume.getNearestBoundedColCoor(x1);
        float[] velocityTrace1 = velocityVolume.getTraceValues(r1, c1);
        float depthDatum = velocityModel.depthDatum;
        float timeDatum = velocityModel.timeDatum;
        float tileZMin = zMin - zShift + tile.colMin * zInc;
        float tileZMax = zMin - zShift + tile.colMax * zInc;
        for(int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float x0 = x1;
            float y0 = y1;
            float[] velocityTrace0 = velocityTrace1;

            x1 = cdpX[row];
            y1 = cdpY[row];
            r1 = velocityVolume.getNearestBoundedRowCoor(y1);
            c1 = velocityVolume.getNearestBoundedColCoor(x1);
            velocityTrace1 = velocityVolume.getTraceValues(r1, c1);

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            float t = tileZMin;

            for(int col = tile.colMin; col <= tile.colMax; col++, t += zInc, colTexCoor += dColTexCoor)
            {
                float v0 = velocityTrace0[col];
                float z0 = (t - timeDatum)*v0 + depthDatum;
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex3f(x0, y0, z0);
                float v1 = velocityTrace1[col];
                float z1 = (t - timeDatum)*v1 + depthDatum;
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                gl.glVertex3f(x1, y1, z1);
            }
            gl.glEnd(

            );
        }
    }
*/
    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
    public boolean getUseShader() { return lineSetClass.getContourColors(); }

    protected boolean initializeTextureTiles(StsGLPanel3d glPanel3d)
    {
        if (textureTiles == null || textureChanged)
        {
            textureTiles = StsTextureTiles.constructor(currentModel, this, nCols, nSlices, isPixelMode);
            if (textureTiles == null) return false;
            textureChanged = true;
//            geometryChanged = true;
        }
        else if (!textureTiles.isSameSize(this))
        {
            textureTiles.constructTiles(this);
            textureChanged = true;
//            geometryChanged = true;
        }

        if (textureTiles.shaderChanged())
        {
            textureChanged = true;
        }
        return true;
    }

    public void textureChanged(boolean changed) {textureChanged = changed;}

    public boolean displayTexture(StsGLPanel3d glPanel3d)
    {
        GL gl = glPanel3d.getGL();
        try
        {
            if(isPixelMode != lineSetClass.getIsPixelMode())
            {
                deleteTexturesAndDisplayLists(gl);
                textureChanged = true;
                isPixelMode = !isPixelMode;
            }
            /*
            if(textureChanged)
            {
                deleteTexturesAndDisplayLists(gl);
                textureTiles = null;
            }
            */
            if (!currentModel.getProject().supportsZDomain(getZDomain()))return false;

            if(!initializeTextureTiles(glPanel3d)) return false;

            byte projectZDomain = currentModel.getProject().getZDomain();
            if (projectZDomain != zDomainDisplayed)
            {
                zDomainDisplayed = projectZDomain;
                deleteDisplayLists(gl);
                textureTiles.constructSurface(this, gl, useDisplayLists, true);
            }

            gl.glDisable(GL.GL_LIGHTING);
            gl.glEnable(GL.GL_TEXTURE_2D);
            gl.glEnable(GL.GL_BLEND);
            gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_FLAT);

            ByteBuffer planeData = null;

            if(!seismicLineSet.setGLColorList(gl, false, textureTiles.shader)) return false;

            if(useDisplayLists && !usingDisplayLists)
                usingDisplayLists = true;
            else if(!useDisplayLists && usingDisplayLists)
            {
                deleteDisplayLists(gl);
                usingDisplayLists = false;
            }

            if(textureChanged)
            {
                planeData = createByteBufferFromFloatData();
                if(planeData == null) return false;
                textureTiles.displayTiles(this, gl, isPixelMode, planeData, nullByte);
                textureChanged = false;
            }
            else
                textureTiles.displayTiles(this, gl, isPixelMode, (ByteBuffer)null, nullByte);

            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "displayTexture", e);
            return false;
        }
        finally
        {
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glDisable(GL.GL_BLEND);
            gl.glEnable(GL.GL_LIGHTING);
            if(textureTiles.shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
        }
    }

    public byte[] readByteData()
    {
        try
        {
            StsMappedByteBuffer byteBuffer = StsMappedByteBuffer.constructor(stsDirectory, rowCubeFilename, "r");
            int nValues = nCols * nSlices;
            if(!byteBuffer.map(0, nValues)) return null;
            byte[] bytes = new byte[nValues];
            byteBuffer.get(bytes);
            byteBuffer.close();
            return bytes;
        }
        catch(Exception e)
        {
            StsException.systemError("StsSeismicLine2d.readFloatData() failed to find file " + stsDirectory + rowCubeFilename);
            return null;
        }
    }

    public ByteBuffer readByteDataBuffer()
    {
        StsMappedByteBuffer stsByteBuffer = createMappedByteRowBuffer("r");
        int nSamples = nCols * nSlices;
        if(!stsByteBuffer.map(0, nSamples)) return null;
        return stsByteBuffer.byteBuffer;
    }

    public ByteBuffer createByteBufferFromFloatData()
    {
        try
        {
            StsMappedFloatBuffer floatBuffer = StsMappedFloatBuffer.openRead(stsDirectory, rowFloatFilename);
            int nValues = nCols * nSlices;
            if(!floatBuffer.map(0, nValues)) return null;
            ByteBuffer byteBuffer = floatBuffer.getScaledByteBuffer(seismicLineSet.dataMin, seismicLineSet.dataMax, nValues);
            floatBuffer.close();
            return byteBuffer;
        }
        catch(Exception e)
        {
            StsException.systemError("StsSeismicLine2d.readFloatData() failed to find file " + stsDirectory + rowFloatFilename);
            return null;
        }
    }

    public float[] readFloatData()
    {
        try
        {
            StsMappedFloatBuffer floatBuffer = StsMappedFloatBuffer.openRead(stsDirectory, rowFloatFilename);
            int nValues = nCols * nSlices;
            if(!floatBuffer.map(0, nValues)) return null;
            float[] values = new float[nValues];
            floatBuffer.get(values);
            floatBuffer.close();
            return values;
        }
        catch(Exception e)
        {
            StsException.systemError("StsSeismicLine2d.readFloatData() failed to find file " + stsDirectory + rowFloatFilename);
            return null;
        }
    }

     public void display(StsGLPanel3d glPanel3d)
     {
        displayTexture(glPanel3d);
     }

    public Class getDisplayableClass() { return StsSeismicLine2d.class; }
/*
    public void display(StsGLPanel3d glPanel3d)
    {
        GL gl = glPanel3d.getGL();

        if(cdpX == null || cdpY == null) return;

        try
        {
            gl.glDisable(GL.GL_LIGHTING);

            gl.glColor3f(0.8f, 0.f, 0.f);
            gl.glLineWidth(2.f);
            gl.glBegin(GL.GL_LINE_STRIP);
            for(int n = 0; n < nCols; n++)
            {
                gl.glVertex3f(cdpX[n], cdpY[n], zMin - zShift);
            }
            for(int n = 0; n < nCols; n++)
            {
                gl.glVertex3f(cdpX[nCols - n - 1], cdpY[nCols - n - 1], zMax - zShift);
            }
            gl.glVertex3f(cdpX[0], cdpY[0], zMin - zShift);
            gl.glEnd();

            gl.glColor3f(1.0f, 0.0f, 0.0f);
            gl.glPointSize(1.f);
            gl.glBegin(GL.GL_POINTS);
            for(int n = 0; n < nCols; n++)
            {
                gl.glVertex3f(cdpX[n], cdpY[n], zMin - zShift);
            }
            for(int n = 0; n < nCols; n++)
            {
                gl.glVertex3f(cdpX[n], cdpY[n], zMax - zShift);
            }
            gl.glEnd();
            if(!displayTexture(glPanel3d))
            {
                gl.glColor3f(0.8f, 0.8f, 0.8f);
                gl.glBegin(GL.GL_QUAD_STRIP);
                for(int n = 0; n < nCols; n++)
                {
                    gl.glVertex3f(cdpX[n], cdpY[n], zMin - zShift);
                    gl.glVertex3f(cdpX[n], cdpY[n], zMax - zShift);
                }
                gl.glEnd();
            }
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "display", e);
        }
        finally
        {
            gl.glEnable(GL.GL_LIGHTING);
        }
    }
*/
/*
    public StsWiggleDisplayProperties getWiggleDisplayProperties()
    {
        return lineSetClass.getWiggleDisplayProperties();
    }
*/
}
