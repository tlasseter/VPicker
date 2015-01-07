package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.util.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author not attributable
 * @version 1.0
 */

/** A seismicCurtain is a series of traces which follow a line or edge.
 *  Points are ordered beginning with down the first trace.
 */

public class StsSeismicCurtain implements StsTextureSurfaceFace, ActionListener
{
	/** seismic volume displayed on curtain */
	StsSeismicVolume seismicVolume;
	/** gridPoints at grid row & col crossings for line */
    public StsGridPoint[] cellGridCrossingPoints;

    transient StsModel model;
	transient StsSeismicVolumeClass seismicVolumeClass;
	transient boolean isPixelMode;
	transient boolean lighting = true;
	transient public float tMin, tMax, tInc;
	transient boolean textureChanged = true;
    transient boolean geometryChanged = true;
    transient boolean deleteTexture = false;
	transient StsTextureTiles textureTiles = null;
	/** rows correspond to traces */
	int nRows = 0;
	/** cols are samples down the trace */
	transient int nCols = 0;
    /** z coordinates (top down) for 2d texture */
    public double[] zCoordinates;
	/** Display lists should be used (controlled by View:Display Options) */
	boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	boolean usingDisplayLists = true;

	public byte zDomain = StsParameters.TD_NONE;

    static final byte nullByte = StsParameters.nullByte;
	static final boolean debug = false;

	public StsSeismicCurtain(StsModel model, StsPoint[] rotatedPoints, StsSeismicVolume seismicVolume)
	{
        this.model = model;
		// Each grid crossing is on a row or column of the seismic volume
        ArrayList gridCrossingPointsList = StsLine.getCellCrossingPoints(seismicVolume, rotatedPoints);
        cellGridCrossingPoints = StsGridCrossings.getStsGridPoints(gridCrossingPointsList);
        computeHorizontalArcLengths();
//        gridCrossingPointsList = line.getGridCrossingPoints(line2d);
//        gridCrossingPoints = StsGridCrossings.getStsGridPoints(gridCrossingPointsList);
        this.seismicVolume = seismicVolume;
        initializeGeometry();

		seismicVolumeClass = (StsSeismicVolumeClass) model.getStsClass(seismicVolume.getClass());
		isPixelMode = seismicVolumeClass.getIsPixelMode();

        seismicVolume.addActionListener(this);
	}

	public StsSeismicCurtain()
	{
	}

    public void initializeGeometry()
    {
       tMin = seismicVolume.getZMin();
       tMax = seismicVolume.getZMax();
       tInc = seismicVolume.getZInc();
       nRows = cellGridCrossingPoints.length;
       nCols = seismicVolume.getNSlices();
    }

    public boolean initialize(StsModel model)
	{
		this.model = model;
		return true;
	}


	public void actionPerformed(ActionEvent e)
	{
		textureChanged = true;
		textureChanged();
	}

    public boolean delete()
    {
        textureChanged();
        //line2d.removeActionListener(this);
		//return super.delete();
		return true;
    }

	/** Different from other texture drawing in that texture rowCoor begins at 0 instead of minRowTexCoor
	 *  because we are not drawing a half-cell on the edge, but a full cell between first and second traces.
	 */
    // TODO may want to add ability to draw DEPTH volume curtain in TIME in general
	public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
		if(zDomain == seismicVolume.getZDomain())
            drawTextureTileTimeSurface3d(tile, gl);
		else
		{
			StsSeismicVelocityModel velocityVolume = model.getProject().velocityModel;
			if (velocityVolume == null)return;
			drawTextureTileDepthSurface3d(velocityVolume, tile, gl);
		}
	}

	private void drawTextureTileTimeSurface3d(StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		double rowTexCoor = 0;
		double dRowTexCoor = tile.dRowTexCoor;
		float tileZMin = tMin + tile.colMin*tInc;
		float tileZMax = tMin + tile.colMax*tInc;
		gl.glBegin(GL.GL_QUAD_STRIP);
		for (int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
            float[] xyz = cellGridCrossingPoints[row].getXYZorT();
			float x = xyz[0];
			float y = xyz[1];

			gl.glTexCoord2d(tile.minColTexCoor, rowTexCoor);
			gl.glVertex3f(x, y, tileZMin);
			gl.glTexCoord2d(tile.maxColTexCoor, rowTexCoor);
			gl.glVertex3f(x, y, tileZMax);
		}
		gl.glEnd();
	}

    public void drawTextureTileDepthSurface3d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        if (velocityModel == null)return;

        if(velocityModel.getInputVelocityVolume() != null)
            drawTextureTileDepthSurfaceFromVolume3d(velocityModel, tile, gl);
        else
            drawTextureTileDepthSurfaceFromIntervalVelocities3d(velocityModel, tile, gl);
    }

    public void drawTextureTileDepthSurfaceFromVolume3d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
		float[] xyz = cellGridCrossingPoints[tile.rowMin].getXYZorT();
		float x1 = xyz[0];
		float y1 = xyz[1];
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

			xyz = cellGridCrossingPoints[row].getXYZorT();
			x1 = xyz[0];
			y1 = xyz[1];
			r1 = velocityVolume.getNearestBoundedRowCoor(y1);
			c1 = velocityVolume.getNearestBoundedColCoor(x1);
			velocityTrace1 = velocityVolume.getTraceValues(r1, c1);

			gl.glBegin(GL.GL_QUAD_STRIP);

			double colTexCoor = tile.minColTexCoor;
			float t = tMin + tile.colMin*tInc;

            for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
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

    public void drawTextureTileDepthSurfaceFromIntervalVelocities3d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
	{
		StsGridCrossingPoint gridCrossingPoint;

		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
		double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
		double dColTexCoor = tile.dColTexCoor;
        StsPoint point = cellGridCrossingPoints[tile.rowMin].getPoint();
		float x1 = point.getX();
		float y1 = point.getY();
        float displayTMin = tMin + tile.colMin*tInc;
        int nSlices = tile.colMax - tile.colMin + 1;
        float[] depths1 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, tInc, nSlices, x1, y1);
        for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
			float x0 = x1;
			float y0 = y1;

			point = cellGridCrossingPoints[row].getPoint();
			x1 = point.getX();
		    y1 = point.getY();
            float[] depths0 = depths1;
            depths1 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, tInc, nSlices, x1, y1);

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

    private void drawTextureTileTimeSurface2d(StsTextureTile tile, GL gl)
    {
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
        double dRowTexCoor = tile.dRowTexCoor;
        float tileZMin = tMin + tile.colMin*tInc;
        float tileZMax = tMin + tile.colMax*tInc;
        gl.glBegin(GL.GL_QUAD_STRIP);
        for (int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float[] xyz = cellGridCrossingPoints[row].getXYZorT();
            float x = xyz[0];
            float y = xyz[1];

            gl.glTexCoord2d(tile.minColTexCoor, rowTexCoor);
            gl.glVertex3f(x, y, tileZMin);
            gl.glTexCoord2d(tile.maxColTexCoor, rowTexCoor);
            gl.glVertex3f(x, y, tileZMax);
        }
        gl.glEnd();
    }

    public void drawTextureTileDepthSurface2d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        if (velocityModel == null)return;

        if(velocityModel.getInputVelocityVolume() != null)
            drawTextureTileDepthSurfaceFromVolume2d(velocityModel, tile, gl);
        else
            drawTextureTileDepthSurfaceFromIntervalVelocities2d(velocityModel, tile, gl);
    }

    public void drawTextureTileDepthSurfaceFromVolume2d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        StsGridCrossingPoint gridCrossingPoint;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        float[] xyz = cellGridCrossingPoints[tile.rowMin].getXYZorT();
        float x1 = xyz[0];
        float y1 = xyz[1];
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

            xyz = cellGridCrossingPoints[row].getXYZorT();
            x1 = xyz[0];
            y1 = xyz[1];
            r1 = velocityVolume.getNearestBoundedRowCoor(y1);
            c1 = velocityVolume.getNearestBoundedColCoor(x1);
            velocityTrace1 = velocityVolume.getTraceValues(r1, c1);

            gl.glBegin(GL.GL_QUAD_STRIP);

            double colTexCoor = tile.minColTexCoor;
            float t = tMin + tile.colMin*tInc;

            for (int col = tile.colMin; col <= tile.colMax; col++, t += tInc, colTexCoor += dColTexCoor)
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

    public void drawTextureTileDepthSurfaceFromIntervalVelocities2d(StsSeismicVelocityModel velocityModel, StsTextureTile tile, GL gl)
    {
        StsGridCrossingPoint gridCrossingPoint;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
//            gl.glColor4f(0.5f, 0.5f, 0.5f, 1.0f);
        double rowTexCoor = 0;
//		double rowTexCoor = tile.minRowTexCoor;
        double dRowTexCoor = tile.dRowTexCoor;
        double dColTexCoor = tile.dColTexCoor;
        StsPoint point = cellGridCrossingPoints[tile.rowMin].getPoint();
        float x1 = point.getX();
        float y1 = point.getY();
        float displayTMin = tMin + tile.colMin*tInc;
        for (int row = tile.rowMin + 1; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
        {
            float x0 = x1;
            float y0 = y1;

            point = cellGridCrossingPoints[row].getPoint();
            x1 = point.getX();
            y1 = point.getY();
            int nSlices = tile.colMax - tile.colMin + 1;
            float[] depths0 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, tInc, nSlices, x0, y0);
            float[] depths1 = velocityModel.getDepthTraceFromIntervalVelocities(displayTMin, tInc, nSlices, x1, y1);

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

    public void displayTexture3d(StsGLPanel3d glPanel3d)
	{
		displayTexture3d(glPanel3d, false);
	}
	
	public void displayTexture3d(StsGLPanel3d glPanel3d, boolean transparent)
	{
		GL gl = glPanel3d.getGL();

        if(!seismicVolumeClass.getIsVisibleOnCurtain()) return;

        if(seismicVolumeChanged())
        {
            if(seismicVolume == null) return;
            deleteTexturesAndDisplayLists(gl);
            initializeGeometry();
        }
         if (isPixelMode != seismicVolumeClass.getIsPixelMode())
         {
			 deleteTexturesAndDisplayLists(gl);
             isPixelMode = !isPixelMode;
			 textureChanged = true;
         }
		 byte projectZDomain = model.getProject().getZDomain();
		 if(projectZDomain != zDomain)
		 {
			 deleteDisplayLists(gl);
			 zDomain = projectZDomain;
			 usingDisplayLists = false;
		 }

		if (textureChanged || textureTiles == null)
		{
			initializeTextureTiles(glPanel3d, gl);
		}

		useDisplayLists = model.useDisplayLists;
		if(!useDisplayLists && usingDisplayLists)
		{
			deleteDisplayLists(gl);
			usingDisplayLists = false;
		}

        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glEnable(GL.GL_TEXTURE_2D);
            //gl.glEnable(GL.GL_BLEND);
            gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
            //gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_FLAT);

            if(transparent)
            {
                gl.glEnable(GL.GL_POLYGON_STIPPLE);
                gl.glPolygonStipple(StsGraphicParameters.getNextStipple(), 0);
            }

            if(!seismicVolume.seismicColorList.setGLColorList(gl, false, textureTiles.shader)) return;

            byte[] data = null;
            if (textureChanged)
            {
                data = seismicVolume.getSeismicCurtainData(cellGridCrossingPoints);
            }

/*
		textureTiles.displayTiles(this, gl, isPixelMode, data, nullByte);
		gl.glDisable(GL.GL_TEXTURE_2D);
		//gl.glDisable(GL.GL_BLEND);
		gl.glEnable(GL.GL_LIGHTING);
        if(transparent) gl.glDisable(GL.GL_POLYGON_STIPPLE);
*/
            if(textureChanged || useDisplayLists && !usingDisplayLists)
            {
                if(textureChanged)       textureChanged = false;
                else if(useDisplayLists) usingDisplayLists = true;
                if(textureTiles == null) StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
                textureTiles.constructSurface(this, gl, useDisplayLists, true);
           }
           textureTiles.displayTiles(this, gl, isPixelMode, data, nullByte);
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "displayTexture", e);
        }
        finally
        {
            gl.glDisable(GL.GL_TEXTURE_2D);
            //gl.glDisable(GL.GL_BLEND);
            gl.glEnable(GL.GL_LIGHTING);
            if(transparent) gl.glDisable(GL.GL_POLYGON_STIPPLE);
        }
	}

    private boolean seismicVolumeChanged()
    {
        StsSeismicVolume currentSeismicVolume = seismicVolumeClass.getCurrentSeismicVolume();
        if(seismicVolume == currentSeismicVolume) return false;
        if(currentSeismicVolume == null) return true;
        seismicVolume.removeActionListener(this);
        seismicVolume = currentSeismicVolume;
        seismicVolume.addActionListener(this);
        return true;
    }

    public int getDefaultShader() { return StsJOGLShader.NONE; }
    public boolean getUseShader() { return false; }
/*
  public void displayTexture(StsGLPanel3d glPanel3d)
  {
	  GL gl = glPanel3d.getGL();

	   if (isPixelMode != lineSet.getIsPixelMode())
	   {
		   deleteTextureTileSurface(gl);
		   isPixelMode = !isPixelMode;
	   }

	  if (textureChanged || textureTiles == null)
	  {
		  checkTextureTiles(glPanel3d);
	  }

	  gl.glDisable(GL.GL_LIGHTING);
	  gl.glEnable(GL.GL_TEXTURE_2D);
	  gl.glEnable(GL.GL_BLEND);
	  gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
	  gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
	  gl.glShadeModel(GL.GL_FLAT);

	  int colorListNum = line2d.getColorDisplayListNum(gl, true);
	  if (mainDebug)
	  {
		  System.out.println("Using colorListNum " + colorListNum);
	  }
	  gl.glCallList(colorListNum);

	  byte[] data = null;
	  if (textureChanged)
	  {
		  data = (byte[]) line2d.getSeismicCurtainData(cellGridCrossingPoints);
		  textureChanged = false;
	  }
	  textureTiles.displayTiles(this, gl, isPixelMode, data);
	  gl.glDisable(GL.GL_TEXTURE_2D);
	  gl.glDisable(GL.GL_BLEND);
	  gl.glEnable(GL.GL_LIGHTING);
  }
*/


	/** This puts texture display on delete list.  Operation is performed
	 *  at beginning of next draw operation.
	 */
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
	/** Called to actually delete the displayables on the delete list. */
	public void deleteTexturesAndDisplayLists(GL gl)
	{
		if(textureTiles == null) return;
		textureTiles.deleteTextures(gl);
		textureTiles.deleteDisplayLists(gl);
		textureChanged = true;
	}

    /** This puts texture display on delete list.  Operation is performed
	 *  at beginning of next draw operation.
	 */

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
        geometryChanged = true;
    }

    private void deleteDisplayLists(GL gl)
	{
		if(textureTiles != null)
			textureTiles.deleteDisplayLists(gl);
	}

/*
    public void deleteTextures()
    {
		if (textureTiles != null)
	   {
		   model.win3d.glPanel3d.addDisplayDeleteableObject(textureTiles);
//			textureTiles.setDeleteTextures();
	   }
       deleteTexture = true;
    }
*/
	protected void initializeTextureTiles(StsGLPanel3d glPanel3d, GL gl)
	{
		// if(!glPanel3d.initialized) return;
		if (textureTiles != null) deleteTexturesAndDisplayLists(gl);
		textureTiles = StsTextureTiles.constructor(model, this, nRows, nCols, isPixelMode);
	}

    public void display(StsGLPanel3d glPanel3d)
    {
        StsArrayList displayableClasses = model.getCursorDisplayableClasses();
        int nDisplayableClasses = displayableClasses.size();
         for(int n = 0; n < nDisplayableClasses; n++)
         {
             StsClassCursorDisplayable displayableClass = (StsClassCursorDisplayable)displayableClasses.get(n);
             displayableClass.drawOn3dCurtain(glPanel3d, this);
         }
    }

    public StsGridPoint[] getCellGridPoints()
    {
        return cellGridCrossingPoints;
    }

    public void computeHorizontalArcLengths()
    {
        int nPoints = cellGridCrossingPoints.length;
        double arcLength = 0.0;
        StsPoint point1 = cellGridCrossingPoints[0].point;

        float[] xy1 = cellGridCrossingPoints[0].point.getXYZ();
        cellGridCrossingPoints[0].point = new StsPoint(xy1, 0.0f);
        for(int n = 1; n < nPoints; n++)
        {
            float[] xy0 = xy1;
            xy1 = cellGridCrossingPoints[n].point.getXYZ();
            double xyDistance = (double)StsMath.distance(xy0, xy1, 2);
            arcLength += xyDistance;
            cellGridCrossingPoints[n].point = new StsPoint(xy1, (float)arcLength);
        }
    }

    public StsSeismicVolume getSeismicVolume() { return seismicVolume; }
	public StsSeismicVolumeClass getSeismicVolumeClass() { return seismicVolumeClass; }
	public float getZMin() { return tMin; }
	public float getZMax() { return tMax; }
	public float getZInc() { return tInc; }
	public int getNCols() { return nRows; }
	public int getNSlices() { return nCols; }
	public byte [] getData2D(StsSeismicVolume seismicV)
	{
		return  seismicV.getSeismicCurtainData(cellGridCrossingPoints);
	}

    public Class getDisplayableClass() { return StsSeismicVolume.class; }


    public String getName()
    {
        return "seismic curtain for: " + seismicVolume.getName();
    }
}