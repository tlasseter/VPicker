package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;

/** A seismicCurtain is a series of traces which follow a line or edge.
 *  Points are ordered beginning with down the first trace.
 */

public class StsTimeLogCurtain implements StsTextureSurfaceFace, ActionListener
{
	/** Dts data displayed on this curtain */
	StsTimeLogCurve dataset;
    StsModel model;
    boolean isPixelMode;
	boolean lighting = true;
	boolean textureChanged = true;
    boolean deleteTexture = false;
	StsTextureTiles textureTiles = null;
    /** Indicates we draw a single line at the current time rather than a curtain of all times. */
    public boolean displayAsLine = false;
    /** rows correspond to traces */
	int nRows = 0;
	/** cols are samples down the trace */
	int nCols = 0;
    /** display width of panel in pixels away from well */
    float displayWidth;
    /** direction indicating which side of the well to draw: 1 for right, -1 for left */
    float direction;
    /** horizontal length of texture display in world coordinates */
    float horizLength;
    /** horizontal increment between traces in world coordinates */
    float horizInc;
    /** world coordinates along well */
    double[][] wellXYZs;
    /** normalized slopes along well */
    double[][] wellSlopes;
    /** coordinates along panel edge */
    double[][] edgeXYZs;
    /** delta between trace points along each line from wellXYZ to edgeXYZ */
    double[][] dXYZs;
    // not implemented yet
	int shader = StsJOGLShader.NONE;
	/** Display lists should be used (controlled by View:Display Options) */
	boolean useDisplayLists;
	/** Display lists currently being used for surface geometry */
	boolean usingDisplayLists = true;

	byte zDomain = StsParameters.TD_NONE;

    StsTimeLogCurveClass datasetClass;

    static final byte nullByte = StsParameters.nullByte;
    static final boolean debug = false;


	public StsTimeLogCurtain()
	{
	}
    private StsTimeLogCurtain(StsModel model, StsTimeLogCurve dataset)
	{
        this.model = model;
        this.dataset = dataset;
		isPixelMode = dataset.getIsPixelMode();
        initialize();
    }

    static public StsTimeLogCurtain constructor(StsModel model, StsTimeLogCurve dataset)
    {
        try
        {
            return new StsTimeLogCurtain(model, dataset);
        }
        catch(Exception e)
        {
            return null;
        }
    }

    private void initialize()
    {
        datasetClass = (StsTimeLogCurveClass)dataset.getStsClass();
        nRows = dataset.getNSurveys();
        StsPoint[] points = dataset.getWellPoints();
        StsPoint[] slopes = dataset.getWellSlopes();
        int nWellPoints = points.length;
        nCols = nWellPoints;
        wellXYZs = new double[nCols][];
        wellSlopes = new double[nCols][];
        edgeXYZs = new double[nCols][];
        dXYZs = new double[nCols][];
        float[] xyz;
        for(int n = 0; n < nCols; n++)
        {
            xyz = points[n].getXYZorT();
            wellXYZs[n] = StsMath.copyDouble(xyz, 3);
            xyz = slopes[n].getXYZorT();
            wellSlopes[n] = StsMath.copyDouble(xyz, 3);
        }
    }

    private void computePanelEdge(StsGLPanel3d glPanel3d)
    {
        for(int n = 0; n < nCols; n++)
        {
            double[] screenPoint = model.getGlPanel3d().getScreenCoordinates(wellXYZs[n]);
            double[] slopePoint = StsMath.multByConstantAddPointStatic(wellSlopes[n], 1000.0, wellXYZs[n]);
            double[] screenSlopePoint = glPanel3d.getScreenCoordinates(slopePoint);

            // screen normal is negative reciprocal of screen slope
            double dsx = -(screenSlopePoint[1] - screenPoint[1]);
		    double dsy =  (screenSlopePoint[0] - screenPoint[0]);

            double s = Math.sqrt(dsx * dsx + dsy * dsy );

            if(s == 0.0)
            {
                dsx = 1.0;
                dsy = 0.0;
            }
            else
            {
                dsx /= s;
                dsy /= s;
            }
            screenPoint[0] += displayWidth*dsx*direction;
            screenPoint[1] += displayWidth*dsy*direction;
            edgeXYZs[n] = glPanel3d.getWorldCoordinates(screenPoint);
            if(!displayAsLine)
            {
                if(direction == 1.0f)
                    dXYZs[n] = StsMath.subtractDivide(edgeXYZs[n], wellXYZs[n], (double)(nRows-1));
                else
                    dXYZs[n] = StsMath.subtractDivide(wellXYZs[n], edgeXYZs[n], (double)(nRows-1));
            }
        }
    }

    public void actionPerformed(ActionEvent e)
	{
		textureChanged = true;
//		addTextureToDeleteList(model.win3d.getGlPanel3d());
	}

    public boolean delete()
    {
//        addTextureToDeleteList(model.win3d.getGlPanel3d());
		return true;
    }

    /** for debugging purposes: call this instead of routine below to draw edges of texture */
    public void drawTextureTileSurfaceLines(StsTextureTile tile, GL gl, boolean is3d)
	{
        // computePanelEdge(glPanel3d);
        StsGLDraw.drawLine(gl, StsColor.RED, false, wellXYZs);
        StsGLDraw.drawLine(gl, StsColor.BLUE, false, edgeXYZs);
    }
    
    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
	{
        if(displayAsLine)
            drawTextureTileLine(tile, gl, is3d);
        else
            drawTextureTileCurtain(tile, gl, is3d);
    }

    private void drawTextureTileCurtain(StsTextureTile tile, GL gl, boolean is3d)
    {
        // computePanelEdge(glPanel3d);
        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		double minRowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
        double minColTexCoor = tile.minColTexCoor;
		double dColTexCoor = tile.dColTexCoor;
        gl.glBegin(GL.GL_QUAD_STRIP);

        double[][] edge;
        if(direction == 1.0f)
            edge = StsMath.copy(wellXYZs);
        else
            edge = StsMath.copy(edgeXYZs);

        double rowTexCoor = minRowTexCoor;
        for (int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
            double colTexCoor = minColTexCoor;
            for(int col = tile.colMin; col <= tile.colMax; col++, colTexCoor += dColTexCoor)
            {
                double[] xyz = edge[col];
                gl.glTexCoord2d(colTexCoor, rowTexCoor);
                gl.glVertex3dv(xyz, 0);
                gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
                StsMath.increment(xyz, dXYZs[col]);
                gl.glVertex3dv(xyz, 0);
            }
        }
		gl.glEnd();
	}

    private void drawTextureTileLine(StsTextureTile tile, GL gl, boolean is3d)
	{
        int nCurtainLine = dataset.checkSetNSelectedSurvey();
        int rowMin = tile.rowMin;
        if(nCurtainLine < rowMin) return;
        int rowMax = tile.rowMax;
        if(nCurtainLine > rowMax) return;

        gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		double minRowTexCoor = tile.minRowTexCoor;
		double dRowTexCoor = tile.dRowTexCoor;
        double minColTexCoor = tile.minColTexCoor;
		double dColTexCoor = tile.dColTexCoor;
        gl.glBegin(GL.GL_QUAD_STRIP);

        int rowOffset = nCurtainLine - rowMin;
        double rowTexCoor = minRowTexCoor + rowOffset*dRowTexCoor;
        double colTexCoor = minColTexCoor;
        for(int col = tile.colMin; col <= tile.colMax; col++, colTexCoor += dColTexCoor)
        {
            gl.glTexCoord2d(colTexCoor, rowTexCoor);
            gl.glVertex3dv(wellXYZs[col], 0);
        //    gl.glTexCoord2d(colTexCoor, rowTexCoor + dRowTexCoor);
            gl.glVertex3dv(edgeXYZs[col], 0);
        }
		gl.glEnd();
	}

    public void displayTextureLine(StsGLPanel3d glPanel3d, float displayWidth, float origin)
    {
        displayAsLine = true;
        displayTexture(glPanel3d, displayWidth, origin);
    }

    public void displayTextureCurtain(StsGLPanel3d glPanel3d, float displayWidth, float origin)
    {
        displayAsLine = false;
        displayTexture(glPanel3d, displayWidth, origin);
    }

    public void displayTexture(StsGLPanel3d glPanel3d, float displayWidth, float origin)
	{
        this.displayWidth = displayWidth;
        if(origin == 0.0f)
            direction = 1.0f;
        else
            direction = -1.0f;
        displayTexture(glPanel3d, false);
        if(!displayAsLine) drawCurrentTime(glPanel3d);
    }

	private void displayTexture(StsGLPanel3d glPanel3d, boolean transparent)
	{
		GL gl = glPanel3d.getGL();

         if (isPixelMode != dataset.getIsPixelMode())
         {
			 deleteTextureAndSurface(gl);
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

		if(!dataset.setGLColorList(gl, false, shader)) return;

    /*
		if(textureChanged || useDisplayLists && !usingDisplayLists)
		{
			if(textureChanged)       textureChanged = false;
			else if(useDisplayLists) usingDisplayLists = true;
			if(textureTiles == null) StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
			textureTiles.constructSurface(this, gl, useDisplayLists, true);
	   }
   */
        computePanelEdge(glPanel3d);

        if(textureChanged)
        {
            textureTiles.displayTiles(this, gl, isPixelMode, dataset.getData(), nullByte);
            textureChanged = false;
        }
        else
            textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null, nullByte);
        gl.glDisable(GL.GL_TEXTURE_2D);
		//gl.glDisable(GL.GL_BLEND);
		gl.glEnable(GL.GL_LIGHTING);
        if(transparent) gl.glDisable(GL.GL_POLYGON_STIPPLE);
	}

    private void drawCurrentTime(StsGLPanel3d glPanel3d)
    {
        GL gl = glPanel3d.getGL();
        gl.glLineWidth(1.0f);
		StsColor.BLACK.setGLColor(gl);
        float f = dataset.getDeltaTFraction();

        if(direction != 1.0f)
            f = 1.0f - f;
        double[][] points = StsMath.interpolatePoints(wellXYZs, edgeXYZs, f);
        if(points == null) return;
        glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
        StsGLDraw.drawLine(gl, StsColor.BLACK, false, points);
        glPanel3d.resetViewShift(gl);
    }


    /** This puts texture display on delete list.  Operation is performed
	 *  at beginning of next draw operation.
	 */
/*
	public void addTextureToDeleteList(StsGLPanel glPanel)
	{
		if (textureTiles != null)
		{
			glPanel.addTextureToDeleteList(this);
		}
		textureChanged = true;
	}
*/
	/** Called to actually delete the displayables on the delete list. */
	public void deleteTextureAndSurface(GL gl)
	{
		if(textureTiles == null) return;
		textureTiles.deleteTextures(gl);
		textureTiles.deleteDisplayLists(gl);
		textureChanged = true;
	}

	private void deleteDisplayLists(GL gl)
	{
		if(textureTiles != null)
			textureTiles.deleteDisplayLists(gl);
	}

	protected void initializeTextureTiles(StsGLPanel3d glPanel3d, GL gl)
	{
		// if(!glPanel3d.initialized) return;
		if (textureTiles != null) deleteTextureAndSurface(gl);
		textureTiles = StsTextureTiles.constructor(model, this, nRows, nCols, isPixelMode);
	}

    public void setTextureChanged()
    {
        textureChanged = true;
    }

    public String getName()
    {
        return "seismic curtain for: " + dataset.getName();
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
    }

    public Class getDisplayableClass() { return StsTimeLogCurve.class; }


    /** Called to actually delete the displayables on the delete list. */
    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if (textureTiles == null) return;
        if (debug) StsException.systemDebug(this, "deleteTextureTileSurface");
        textureTiles.deleteTextures(gl);
        textureTiles.deleteDisplayLists(gl);
        textureChanged = true;
    }

    public boolean getUseShader() { return datasetClass.getContourColors(); }
    public int getDefaultShader() { return StsJOGLShader.NONE; }

    protected boolean getIsPixelMode() { return datasetClass.getIsPixelMode(); }

    protected void setGLColorList(GL gl)
    {
        dataset.setGLColorList(gl, false, StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS);
    }

    protected void clearShader()
    {
    }

    protected byte[] getData()
    {
        return dataset.getData();
    }
}
