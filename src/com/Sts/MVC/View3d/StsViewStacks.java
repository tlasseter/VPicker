package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class defining two-dimensional cursor view. The cursor is a reference to the
 * three planes that are displayed in the 3D view. An object based on this class wuold result in a 2D
 * view of any one of the three cursor planes.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.Seismic.*;
import com.Sts.Utilities.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import java.awt.*;
import java.nio.*;

abstract public class StsViewStacks extends StsViewPreStack implements StsTextureSurfaceFace, StsSerializable
{
    public float[][] dataRanges = new float[2][2];
    public float[][] originalAxisRanges = new float[2][2];
	public boolean isPixelMode;
    /** temporary gather used in computing stacked traces for cvs */
    transient protected StsSuperGather cvsGather;
    /** CVS properties have changed: recompute */
    transient public boolean cvsPropertiesChanged = true;
	/** compass angle displayed in legend at top of 2D view */
	transient public float displayAngle;
	/** legend display type at top of 2D view */
	transient public byte legendType = S2S_LEGEND;

//    float minVelocity = 5000.0f;
//    float maxVelocity = 15000.0f;
    transient byte displayType;

    double[][] traceData = null; // [nCroppedSlices][2]
    byte[] tracePointTypes;
    double maxError = 0.01f;

    /** Tiles on which texture is generated */
	transient public StsTextureTiles textureTiles;
    transient public int nTextureRows, nTextureCols;
    transient int nBackgroundRows, nBackgroundCols;

    transient private boolean textureChanged = true;
	/** Display lists should be used (controlled by View:Display Options) */
	transient boolean useDisplayLists = true;
	/** Display lists currently being used for surface geometry */
	transient boolean usingDisplayLists = false;
//    transient public float[][][] cvsTraces = null;
    transient int nTracesInView = 0;
    /** trace data converted to texture data */
    transient ByteBuffer panelData;
    /** transparent trace inserted between panels */
    transient byte[] transparentTrace;
    /** scratch byte traces */
    transient byte[] traceBytes;
    /** display list used for wiggle traces and variable area */
    transient int displayListNum = 0;
    /** Indicates display has changed and display list needs to be deleted */
    transient boolean displayListChanged = true;
    /** Number of panels in display: each panel is constant velocity (CVS) or a perturbed velocity profile (VVS) */
    transient int nPanels;
    /** Number of traces in each panel */
    transient int nPanelTraces;
    /** Spacing in deltaVelocity per trace */
    transient float horizontalTraceSpacing;
    /** deltaVelocity between panels (includes gap) */
    transient float velocityStep;
    /** velocity of first panel */
    transient float velocityMin;
   /** velocity of last panel */
    transient float velocityMax;
    
    transient boolean ignoreSuperGather;

	static final byte CORE_LABS_LEGEND = 1;
	static final byte S2S_LEGEND = 2;

	static final String[] compassDirections = new String[]
		{"E", "NE", "N", "NW", "W", "SW", "S", "SE"};

    transient boolean runTimer = false;
    transient StsTimer timer;

    abstract boolean rangeChanged();
    abstract protected void stackPanelsChanged();
    abstract protected void recomputeAxisRanges();
    abstract protected boolean computeStacksData();
    abstract protected float[][][] getData();
    abstract protected boolean  isDataNull();
    /**
	 * Default constructor
	 */
	public StsViewStacks()
	{
//		limitPan = true; // no-arg constructor used only by database when loading in which case limitPan will be read in and overwrites this.  TJL 1/20/07
	}

	/**
	 * StsPreStackView2d constructor
	 * @param glPanel3d the graphics context
	 */
	public StsViewStacks(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
		limitPan = true;
        setLineSet(StsPreStackLineSetClass.currentProjectPreStackLineSet);
		initialize();
		// matchLockedWindows();
	}

    public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        isViewGLInitialized = true;
    }

    public void initializeTransients(StsGLPanel3d glPanel3d)
	{
		super.initializeTransients(glPanel3d);
        setLineSet(StsPreStackLineSetClass.currentProjectPreStackLineSet);
        initialize();
        initializeAxisLabels();
        recomputeAxisRanges();
        limitPan = true;
	}

    protected void initializeAxisLabels()
    {
        axisLabels = new String[] {"Velocity", "Time"};
	}

	public void initialize()
	{
        lineSetClass = (StsPreStackLineSetClass)model.getCreateStsClass(lineSet);
        isPixelMode = lineSetClass.getIsPixelMode();
		superGather = lineSet.getSuperGather(glPanel3d.window);
        if(superGather != null)
        superGather.initializeSuperGather(superGather.superGatherRow, superGather.superGatherCol);
        initializeAxisLabels();
		setAxisRanges();
        setMouseModeZoom();
        // StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
        // if(toolbar != null) toolbar.zoom();
	}

	public void computeProjectionMatrix()
	{
		if(axisRanges == null)return;
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
		glu.gluOrtho2D(axisRanges[0][0], axisRanges[0][1], axisRanges[1][0], axisRanges[1][1]);
		gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
    }

	/**
	 * Is the cursor viewable
	 * @return true if viewable
	 */
	public boolean isViewable()
	{
		return true;
	}

    public boolean initializeDefaultAction()
    {
        setDefaultAction(null);
		return true;
    }

    protected void clearView()
    {
        //StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        //glPanel3d.setClearColor(wiggleProperties.getWiggleBackgroundColor());
        //gl.glDrawBuffer(GL.GL_BACK);
        // glPanel3d.applyClearColor();
        gl.glDisable(GL.GL_LIGHTING);
        gl.glDisable(GL.GL_DEPTH_TEST);
        // gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
        clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
    }

    protected void displayData(GL gl, GLU glu)
    {
        if(superGather == null) return;
        byte displayType = lineSet.getCVSProperties().getDisplayType();
        if(this.displayType != displayType)
        {
            if(this.displayType == StsCVSProperties.DISPLAY_VAR_DENSITY)
                StsTextureList.addTextureToDeleteList(this);
            else
                deleteDisplayLists(gl);
        }
        this.displayType = displayType;
        if(displayType == StsCVSProperties.DISPLAY_VAR_DENSITY)
            displayPanelTexture(gl, glu);
        else
            displayPanel(gl, glu, displayType);
    }

    protected synchronized boolean checkComputeData()
    {
        if(isCursor3dDragging()) return false;

        if(superGather.computingStacks) return true;
        if(!textureChanged) return false;
        if(!isDataNull()) return false;
        if(debug) debugThread("computeCVSProcess started");
        return computeStacksData();
    }

    public void adjustCursor(int dir, float dirCoor)
	{
		if(lineSet == null)return;
        if(superGather != null)
        superGather.resetReprocessSemblanceBlockFlag();
        textureChanged();
    }

    protected boolean  isDataChanged() { return this.cvsPropertiesChanged; }

    protected void debugThread(String message)
    {
        System.out.println(Thread.currentThread().getName() + " window " + glPanel3d.window.getTitle() + " view class " + StsToolkit.getSimpleClassname(this) + " " + message +
                    " computingSemblance " + superGather.computingSemblance + " currentRow " + superGather.superGatherRow + " currentCol " + superGather.superGatherCol);
    }

    public boolean viewObjectChanged(Object source, Object object)
	{
        // TODO determine more precisely when we need to set cvsPropertiesChanged to true
//        cvsPropertiesChanged = true;
        if(object instanceof StsPreStackLineSet)
		{
			gatherTracesChanged();
            setLineSet((StsPreStackLineSet)object);
            return true;
		}
		else if(object instanceof StsPreStackVelocityModel)
		{
			return true;
		}
        else if(object instanceof StsAGCProperties || object instanceof StsFilterProperties)
        {
            gatherTracesChanged();
            return true;
        }
		else if(object instanceof StsWiggleDisplayProperties)
		{
            StsWiggleDisplayProperties wiggleProperties = (StsWiggleDisplayProperties)object;
            if(wiggleProperties.isRangeChanged())
            {
                setAxisRanges();
                computePixelScaling();
                gatherDisplayChanged();
                return true;
            }
            if(wiggleProperties.stretchMuteChanged)
            {
                gatherDisplayChanged();
                return true;
            }
            else
            {
                gatherDisplayChanged();
                return true;
            }
        }
        else if(object instanceof StsCVSProperties)
        {
            StsCVSProperties cvsProperties = (StsCVSProperties)object;
            if(cvsProperties.isRangeChanged() && rangeChanged())
			{
                // setAxisRanges();
                recomputeAxisRanges();
                computePixelScaling();
                gatherTracesChanged();
            }
            else if(cvsProperties.ignoreSuperGatherChanged())
                gatherTracesChanged();
            else if(cvsProperties.isChanged())
            {
                cvsPropertiesChanged = true;
                gatherDisplayChanged();
            }
            return true;
		}
        else if(object instanceof StsSuperGatherProperties)
        {
            setAxisRanges();
            recomputeAxisRanges();
            gatherTracesChanged();
            return true;
		}
        else if(object instanceof StsDatumProperties)
        {
            if(superGather.checkSetDatumShift())
            {
                gatherTracesChanged();
                return false;
            }
            return true;
        }
        else if(object instanceof StsSemblanceRangeProperties)
        {
            StsSemblanceRangeProperties rangeProperties = (StsSemblanceRangeProperties)object;
            if(rangeProperties.isRangeChanged() && rangeChanged())
            {
                stackPanelsChanged();
                setAxisRanges();
                recomputeAxisRanges();
                gatherTracesChanged();
            }
        }
        else if (object instanceof StsColorList)
        {
            if(object == lineSet.seismicColorList)
                textureChanged();                       //Note: never called. ActionPerformed event handled by StsPreStackLineSet only calls viewObjectRepaint() SWC 8/6/09
            return true;
        }
        return false;
	}

    public boolean viewObjectRepaint(Object source, Object object)
	{
		if(object instanceof StsPreStackLineSet)
		{
//			currentRow = superGather.superGatherRow;
//			currentCol = superGather.superGatherCol;
            glPanel3d.repaint();
			return true;
		}
		else if(object instanceof StsVelocityProfile)
		{
			glPanel3d.repaint();
			return true;
		}
        else if(object instanceof StsPanelProperties)  // Catch all property changes.
        {
            glPanel3d.repaint();
            return true;
		}
        else if (object == lineSet.seismicColorList )
        {
            textureChanged();       //textureChanged() must be called here because viewObjectChanged() never called. (see note in viewObjectChanged) SWC 8/6/09
            glPanel3d.repaint();
            return true;
        }
        return false;
	}

    protected StsColor getGridColor()
	{
		StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
		return new StsColor(wiggleProperties.getLineColor());
	}

	/*
	 private boolean checkSetCurrentRowCol()
	 {
	  StsGLPanel3d parentGlPanel3d = model.getGLPanel3d(glPanel3d.window);
	  StsCursor3d cursor3d = parentGlPanel3d.cursor3d;
	  return cursor3d.rowOrColChanged(currentRowCol);
	 }
	 */
        /*
    private int getPixelsPerInch()
    {
        int pixelsPerInch = 96;
        try
        {
            pixelsPerInch = Toolkit.getDefaultToolkit().getScreenResolution();
        }
        catch(Exception e)
        { }
        return pixelsPerInch;
    }
*/
/*
    private void rescaleInchesPerSecond()
    {
        // with our current size, how many samples can we fit ?
        float inches = (float)glPanel3d.glc.getHeight() / (float)getPixelsPerInch();
        float nSeconds = inches / inchesPerSecond;

        axisRanges[1][1] = totalAxisRanges[1][1];
        axisRanges[1][0] = totalAxisRanges[1][1] + (nSeconds * 1000.0f);
    }
 */
/*
	private void rescaleTracesPerInch()
	{
		StsCVSProperties cvsProperties = lineSet.cvsProperties;
		// with our current size, how many traces can we fit ?
		float inches = (float)glPanel3d.glc.getWidth() / (float)getPixelsPerInch();
		float nt = inches * cvsProperties.getTracesPerInch();

		axisRanges[0][0] = totalAxisRanges[0][0];
		axisRanges[0][1] = totalAxisRanges[0][0] + nt;
        if(axisRanges[0][1] > totalAxisRanges[0][1])
            axisRanges[0][1] = totalAxisRanges[0][1];
	}
*/
	public void setAxisRanges()
	{
	    //ignoreSuperGather = lineSet.cvsProperties.ignoreSuperGather;
	    recomputeAxisRanges();
//       rescaleInchesPerSecond();
//       matchLockedWindows();
//       glPanel3d.viewChanged = true;
	}

    public float[][] getAxisRanges()
    {
        float[][] ranges = new float[2][2];
        if(axisRanges[0][0] < originalAxisRanges[0][0])
            ranges[0][0] = originalAxisRanges[0][0];
        else
            ranges[0][0] = axisRanges[0][0];

        if(axisRanges[0][1] > originalAxisRanges[0][1])
            ranges[0][1] = originalAxisRanges[0][1];
        else
            ranges[0][1] = axisRanges[0][1];

        ranges[1] = new float[] {axisRanges[1][0], axisRanges[1][1]};
        return ranges;
	}

	/**
	 * Output the mouse tracking readout to the information panel on the main screen
	 * @param glPanel3d
     * @param mouse mouse object
	 */
	public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse)
	{
		float[] xyz = new float[3];

		StsCursorPoint cursorPoint = glPanel3d.getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
		if(cursorPoint == null) return null;

		StsPoint point = glPanel3d.getPointInPlaneAtMouse(mouse);
		xyz[0] = glPanel3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.XDIR);
		xyz[1] = glPanel3d.getCursor3d().getCurrentDirCoordinate(StsCursor3d.YDIR);
		xyz[2] = point.v[1];

		setCursorXOR(glPanel3d, mouse, xyz);

		int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
		if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
			StsMessageFiles.infoMessage(axisLabels[0] + "= " + point.v[0] +
										" " + model.getProject().getZDomainString() + "= " + point.v[1]);
		return cursorPoint;
	}

    private void displayPanel(GL gl, GLU glu, byte displayType)
    {
        if(displayListChanged)
            deleteDisplayLists(gl);
        useDisplayLists = StsObject.getCurrentModel().useDisplayLists;
		if(!useDisplayLists && usingDisplayLists)
			deleteDisplayLists(gl);

        if (displayListNum == 0 && useDisplayLists) // build display list
        {
            displayListNum = gl.glGenLists(1);
            if (displayListNum == 0)
            {
                StsMessageFiles.logMessage("System Error in StsGrid.displayGrid: " +
                                        "Failed to allocate a display list");
                return;
            }
            gl.glNewList(displayListNum, GL.GL_COMPILE_AND_EXECUTE);
            float[][][] displayData = getData();
            displayPanel(displayData, gl, glu, displayType);
            gl.glEndList();
        }
        else if (useDisplayLists) // use existing display list
        {
            gl.glCallList(displayListNum);
        }
        else // immediate mode draw
        {
            if (displayListNum > 0)
            {
                gl.glDeleteLists(displayListNum, 1);
                displayListNum = 0;
            }
            float[][][] displayData = getData();
            displayPanel(displayData, gl, glu, displayType);
        }
    }

    private void displayPanel(float[][][] displayData, GL gl, GLU glu, byte displayType)
    {
       // StsCVSProperties cvsProperties = lineSet.cvsProperties;
        StsWiggleDisplayProperties panelWiggleProperties = lineSet.getPanelWiggleDisplayProperties();

        int overlapPercent = panelWiggleProperties.getWiggleOverlapPercent();
        boolean wiggleSmoothCurve = panelWiggleProperties.getWiggleSmoothCurve();
        boolean reversePolarity = panelWiggleProperties.getWiggleReversePolarity();

        float horizScale = computeWiggleHorizontalScaling(overlapPercent);  // horizontal size per trace in velocity units
        if (reversePolarity) horizScale = -horizScale;
        float[] scaleAndOffset = StsMath.floatToNormalizedFloatScaleAndOffset(lineSet.dataMin, lineSet.dataMax); 
        float valueOffset = scaleAndOffset[1];
        double rmsAmp = StsMath.rmsIgnoreZero(displayData);
        float valueScale = (float)(1.0/rmsAmp);//scaleAndOffset[0]; //we don't scale to absolute min/max of lineset - could be spikes! RMS much more stable SWC 11/23/09
        valueScale *= horizScale;
        valueOffset *= horizScale;

        int nSlices = superGather.nVolumeSlices;
        float halfPanelWidth = horizontalTraceSpacing*(nPanelTraces-1)/2;
        float zMax = lineSet.zMax;
        float zInc = lineSet.zInc;
        float zMin = lineSet.zMin;
        float velocity = velocityMin;
        double[] muteRange = new double[] { 0.0, zMax };
        int[] displayDataRange = new int[] {0, nSlices-1};//StsTraceUtilities.getGatherDataRange(axisRanges, traceData); //traceData always null!! SWC 11/20/09
        int displayMin = displayDataRange[0];
        int displayMax = displayDataRange[1];
        int nDisplaySlices = displayMax - displayMin + 1;
        float displayZMin = zMin + displayMin*zInc;
        float displayZMax = zMin + displayMax*zInc;
        for(int i = 0; i < nPanels; i++, velocity += velocityStep)
        {
            float x0 = velocity - halfPanelWidth;
            for(int j = 0; j < nPanelTraces; j++, x0 += horizontalTraceSpacing)
            {
                if(displayData[i][j] == null)  // No gather or no traces in gather at this location
                    drawFlatLine(gl, x0, displayZMin, displayZMax);
                else
                    StsTraceUtilities.displayInterpolatedPoints(gl, displayData[i][j], x0, displayMin, displayZMin, zInc, nDisplaySlices, valueScale, valueOffset, panelWiggleProperties, muteRange, 1);
            }
        }
//        cvsPropertiesChanged = false;
    }

    private void drawFlatLine(GL gl, float x, float zMin, float zMax)
    {
        StsColor.BLACK.setGLColor(gl);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glBegin(GL.GL_LINES);
        gl.glVertex2f(x, zMin);
        gl.glVertex2f(x, zMax);
        gl.glEnd();
        gl.glEnable(GL.GL_LIGHTING);
    }

    protected void displayPanelTexture(GL gl, GLU glu)
    {
        superGather = lineSet.getSuperGather(glPanel3d.window);

        if(isPixelMode != lineSetClass.getIsPixelMode())
        {
            isPixelMode = !isPixelMode;
            StsTextureList.addTextureToDeleteList(this);
        }

        // size of textureTiles may have changed above, so check and rebuild if necessary
        if(!checkTextureTiles(gl)) return;
        if(textureTiles.shaderChanged()) textureChanged = true;
        
        try
        {
            gl.glDisable(GL.GL_LIGHTING);
            gl.glEnable(GL.GL_TEXTURE_2D);
            gl.glEnable(GL.GL_BLEND);
            //        gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
            gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
            gl.glShadeModel(GL.GL_SMOOTH);
//			gl.glDisable(GL.GL_DEPTH_TEST);

            if(lineSet == null || lineSet.seismicColorList == null)
                return;
//            if(!lineSet.seismicColorList.setGLColorList(gl, true, shader))
            if(!lineSet.seismicColorList.setGLColorList(gl, true, textureTiles.shader, Color.WHITE))
                return;

            useDisplayLists = StsObject.getCurrentModel().useDisplayLists;
            if(!useDisplayLists && usingDisplayLists)
            {
                deleteDisplayLists(gl);
            }
            else if(useDisplayLists && !usingDisplayLists)
            {
                if(textureTiles == null)
                    StsException.systemError("StsSurface.displaySurfaceFill() failed. textureTiles should not be null.");
                textureTiles.constructSurface(this, gl, useDisplayLists, false);
                usingDisplayLists = true;
            }

            if(debug) System.out.println("draw CVS texture");
            if(textureTiles.shader != StsJOGLShader.NONE)
                StsJOGLShader.enableARBShader(gl, textureTiles.shader);
            if(textureChanged && !isDataNull())
            {
                if(debug) debugThread("rebuilding new texture");
                textureTiles.displayTiles2d(this, gl, false, isPixelMode, getPanelData(), nullByte);
                textureChanged = false;
//				gather.semblanceBytes = null;
            }
            else
            {
                if(debug)
                    debugThread("displaying current texture");
                textureTiles.displayTiles2d(this, gl, false, isPixelMode, (byte[])null, nullByte);
            }
            if(runTimer) timer.stopPrint("display CVS texture");
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
        finally
        {
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glDisable(GL.GL_BLEND);
            if(textureTiles.shader != StsJOGLShader.NONE)
                StsJOGLShader.disableARBShader(gl);
            gl.glDisable(GL.GL_TEXTURE_2D);
            gl.glEnable(GL.GL_LIGHTING);
        }
    }

    public int getDefaultShader() { return StsJOGLShader.NONE; }
    public boolean getUseShader() { return false; }

    protected boolean checkTextureTiles(GL gl)
    {
        if(lineSet == null) return false;
        if(textureTiles != null && !tilesChanged()) return true;

        float[][] textureRanges = new float[2][2];
        textureRanges[0][0] = totalAxisRanges[1][1];
        textureRanges[0][1] = totalAxisRanges[1][0];
        textureRanges[1][0] = totalAxisRanges[0][0];
        textureRanges[1][1] = totalAxisRanges[0][1];
        if(textureTiles != null) deleteTexturesAndDisplayLists(gl);
        textureTiles = StsTextureTiles.constructor(model, this, nTextureRows, nTextureCols, isPixelMode, textureRanges);
        textureChanged = true;
//			if (textureTiles == null)return;
//			nTextureRows = textureTiles.nTotalRows;
//			nTextureCols = textureTiles.nTotalCols;
        return true;
    }

    public void deleteDisplayLists(GL gl)
    {
        if (textureTiles != null)
        {
            textureTiles.deleteDisplayLists(gl);
        }
        if (displayListNum > 0)
        {
            gl.glDeleteLists(displayListNum, 1);
            displayListNum = 0;
        }
        displayListChanged = false;
        usingDisplayLists = false;
    }

    private ByteBuffer getPanelData()
    {
        StsCVSProperties cvsProperties = lineSet.cvsProperties;
        StsSemblanceRangeProperties rangeProperties = lineSet.semblanceRangeProperties;
        float[][][] cvsTraces = getData();
        int nPanels = cvsProperties.getNumberPanels();
        int nSlices = superGather.nVolumeSlices;
        checkScratchTraces(nSlices);
        panelData.rewind();
        double gainFactor = lineSet.getCVSProperties().getWiggleOverlapPercent()/100.0;
        for(int p = 0; p < nPanels; p++)
        {
            panelData.put(transparentTrace);
            panelData.put(transparentTrace);
            for(int t = 0; t < nPanelTraces; t++)
            {
                StsSeismicBoundingBox.scaleNormalizedTraceFloatsToBytes(StsMath.arrayMult(cvsTraces[p][t], gainFactor), traceBytes);
                panelData.put(traceBytes);
            }
 //           StsSeismicBoundingBox.scaleNormalizedTraceFloatsToBytes(cvsTraces[p][nPanelTraces-1], traceBytes0);
 //           StsMath.interpolateUnsignedBytes254(traceBytes0, traceBytes1, 1.5f, nCroppedSlices, traceBytes);
            panelData.put(transparentTrace);
        }
        panelData.put(transparentTrace);
        panelData.rewind();
        return panelData;
    }

    private void checkScratchTraces(int nSlices)
    {
        if(transparentTrace != null && transparentTrace.length == nSlices) return;
        transparentTrace = StsToolkit.constructByteArray(nSlices, (byte)-1);
        traceBytes = new byte[nSlices];
    }
    
    private boolean tilesChanged()
    {
        if(textureTiles == null) return true;
        if(textureTiles.nTotalRows != nTextureRows) return true;
        if(textureTiles.nTotalCols != nTextureCols) return true;
        float[][] textureRanges = textureTiles.axisRanges;
        if(textureRanges[0][0] != dataRanges[1][1]) return true;
        if(textureRanges[0][1] != dataRanges[1][0]) return true;
        if(textureRanges[1][0] != dataRanges[0][0]) return true;
        if(textureRanges[1][1] != dataRanges[0][1]) return true;
        return false;
    }

    protected void checkTextureSizeAndRange()
	{
        if(!textureChanged) return;
        if(debug) debugThread("setting texture range");
        nTextureRows = (nPanels * (nPanelTraces + 3)) + 1;
        nTextureCols = superGather.nVolumeSlices;
        int size = nTextureRows*nTextureCols;
        if(panelData != null && panelData.capacity() == size) return;
        panelData = BufferUtil.newByteBuffer(size);
    }

    /** computes horizontalScaling of normalized amplitude in units of velocity.
     * An overlap of zero puts the max amplitude half way between two traces; 100 percent has amplitude touching baseline of next trace
     */
    private float computeWiggleHorizontalScaling(float overlapPercentage)
    {
        //return ((overlapPercentage + 100.0f)/200.0f)*horizontalTraceSpacing;
        return (overlapPercentage/100)*horizontalTraceSpacing;
    }


    public double[][] getWigglePoints(float[] floats, int panelIdx, int traceIdx)
    {
        int nSlices = superGather.nVolumeSlices;
        int minSlice = 0;
        double[][] points = new double[nSlices+1][3];
        float tMin = lineSet.getZCoor(minSlice);
        float tInc = lineSet.getZInc();
        int nPoints = nSlices;
        int index = 0;
        boolean started = false;
        double v1 = floats[index++];

//        if (Math.abs(v1) < StsParameters.roundOff) v1 = 0;
        int n = 0;
        float t = tMin;
        int nLastNonZero = 0;
 //       points[0][0] = t;
 //       points[0][1] = v1;
 //       n=1;

        for (int s = 1; s < nSlices; s++, t += tInc)
        {
            double v = floats[index++];
            if (Math.abs(v) < StsParameters.roundOff) v = 0;

            if (!started)
            {
                double v0 = v1;
                v1 = v;
                if (v1 != 0)
                {
                    started = true;
                    points[n][0] = t - tInc;
                    points[n][1] = v0;
                    n++;
                }
            }
            if (started)
            {
                points[n][0] = t;
                points[n][1] = v;
                if (v != 0) nLastNonZero = n;
                n++;
            }
        }
        nPoints = n;
        int nTrimmedPoints = nLastNonZero + 1;
        if (nTrimmedPoints < 2)
        {
            points = new double[2][3];
            points[0][0] = tMin;
            points[1][0] = tMin + tInc * (nSlices - 1);
            nPoints = 2;
        }
        else if (nTrimmedPoints < nPoints)
        {
            points = (double[][]) StsMath.trimArray(points, nTrimmedPoints + 1); // include a zero point at the end
            nPoints = nTrimmedPoints + 1;
        }

        // compute slopes (0 curvature as end conditions)
        for (n = 1; n < nPoints - 1; n++)
            points[n][2] = (points[n + 1][1] - points[n - 1][1]) / 2;
        points[0][2] = 1.5 * (points[1][1] - points[0][1]) - points[1][2] / 2;
        points[nPoints - 1][2] = 1.5 * (points[nPoints - 1][1] - points[nPoints - 2][1]) - points[nPoints - 2][2] / 2;

        return points;
    }

    /** panel or traces have been fundamentally changed, so we need to recompute panel */
    protected void gatherTracesChanged()
    {
        superGather.gatherTracesChanged();
        gatherDisplayChanged();
    }

    /** panel or trace display properties have changed, so we don't need to recompute panel, but do need to redisplay it */
    protected void gatherDisplayChanged()
    {
        textureChanged();
        glPanel3d.viewChanged = true;
        displayListChanged = true;
    }

    public void drawForeground(GL gl)
	{
		gl.glDisable(GL.GL_BLEND);
		gl.glDrawBuffer(GL.GL_FRONT);
		gl.glEnable(GL.GL_COLOR_LOGIC_OP);
		gl.glLogicOp(GL.GL_XOR);
		gl.glDepthMask(false);
		gl.glDepthFunc(GL.GL_ALWAYS);
		drawForegroundCursor();
		gl.glFlush();
		if(cursorButtonState == StsMouse.RELEASED)
		{
			gl.glLogicOp(GL.GL_COPY);
			gl.glDisable(GL.GL_COLOR_LOGIC_OP);
			gl.glDepthMask(true);
			gl.glDepthFunc(GL.GL_LESS);
			gl.glDrawBuffer(GL.GL_BACK);
			cursorButtonState = StsMouse.CLEARED;
		}
	}

	/** Draw the cursor in the front buffer.  If previously drawn,
	 *  draw previous one again to erase it; then drawn new one.
	 */
	private void drawForegroundCursor()
	{
		if(previousXYZ != null)
			drawCursorPoint(previousXYZ, this.glPanel3d);
		if(currentXYZ != null)
			drawCursorPoint(currentXYZ, glPanel3d);
		previousXYZ = currentXYZ;
		currentXYZ = null;
	}

	/** Draw  horizontal lines thru the cursor point.
	 */
	private void drawCursorPoint(float[] xyz, StsGLPanel3d glPanel3d)
	{
        System.out.println("Draw cursor on CVS X=" + axisRanges[0][0] + " Y=" + xyz[2]);
		StsColor.GRAY.setGLColor(gl);
		gl.glLineWidth((float)3.f);
		gl.glBegin(GL.GL_LINES);
		gl.glVertex2f((float)axisRanges[0][0], xyz[2]);
		gl.glVertex2f((float)axisRanges[0][1], xyz[2]);
		gl.glEnd();

	}

	/** For Corelabs legend type, display the line direction for a vertical section at the top;
	 *  and the view direction for a horizontal section.  Line direction is the
	 *  compass direction to the left and right.  View direction is the compass direction
	 *  in the direction from bottom to top of the view.  Compass direction is one
	 *  of the 8 compass points.
	 *  <p>For S2S legend type, display the view direction and angle for both horizontal
	 *  and vertical sections.
	 *
	 */
	public void drawTopLegend(GL gl)
	{
		int x, y;

		if(legendType == CORE_LABS_LEGEND)
		{
			String[] legends = getClbView2dLegend();
			y = glPanel3d.getHeight() - insets.top / 2 + halfWidth;
			if(legends.length == 2)
			{
				x = insets.left;
				StsGLDraw.fontOutput(gl, x, y, legends[0], horizontalFont);
				x = glPanel3d.getWidth() - insets.right - StsGLDraw.getFontStringLength(horizontalFont, legends[1]);
				StsGLDraw.fontOutput(gl, x, y, legends[1], horizontalFont);
			}
			else if(legends.length == 1)
			{
				x = glPanel3d.getWidth() - insets.right - StsGLDraw.getFontStringLength(horizontalFont, legends[0]) - 150;
				StsGLDraw.fontOutput(gl, x, y, legends[0], horizontalFont);
			}
		}
		else // S2S_LEGEND
		{
			String legend = getStsView2dLegend();
			int axisMaxX = glPanel3d.getWidth() - insets.right;
			x = axisMaxX - StsGLDraw.getFontStringLength(horizontalFont, legend);
			y = glPanel3d.getHeight() - insets.top / 2 + halfWidth;
			StsGLDraw.fontOutput(gl, x, y, legend, horizontalFont);
		}
	}

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        drawTextureTileSurface(tile, gl);
    }

    private void drawTextureTileSurface(StsTextureTile tile, GL gl)
	{
		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		double rowTexCoor = 0;
		double dRowTexCoor = tile.dRowTexCoor;
        float zInc = lineSet.zInc;
        float zMin = lineSet.zMin;
        float tileZMin = zMin + tile.colMin*zInc;
		float tileZMax = zMin + tile.colMax*zInc;
        gl.glBegin(GL.GL_QUAD_STRIP);
		for (int row = tile.rowMin; row <= tile.rowMax; row++, rowTexCoor += dRowTexCoor)
		{
			float x = getTileX(row);
			gl.glTexCoord2d(tile.minColTexCoor, rowTexCoor);
			gl.glVertex2f(x, tileZMin);
			gl.glTexCoord2d(tile.maxColTexCoor, rowTexCoor);
			gl.glVertex2f(x, tileZMax);
		}
		gl.glEnd();
	}

    /** VS view consists of a series of panels with equally-spaced traces
     *  with a gap on left, right, and between panels with the same trace spacing.
     *  For texture, each panel consists of the given number of traces;
     *  each gap consists of two null traces.  The position of these gap traces is
     *  the same as the panel on the left and right (i.e., two traces at the boundary between
     *  panel and gap are at the same X location in window.
     *
     * @param row
     * @return
     */
    private float getTileX(int row)
    {
        int nPanel = (row + 1)/(nPanelTraces + 3);
        float velocityMin = totalAxisRanges[0][0];
        int nPanelTrace = row - nPanel*(nPanelTraces + 3) - 2;
        if(nPanelTrace == -1)
            return velocityMin + horizontalTraceSpacing + nPanel*velocityStep;
        if(nPanelTrace == -2)
            return velocityMin + nPanel*velocityStep;
        if(nPanelTrace == -3)
            return velocityMin - horizontalTraceSpacing + nPanel*velocityStep;
        else
            return velocityMin + horizontalTraceSpacing + nPanel*velocityStep + nPanelTrace*horizontalTraceSpacing;
    }

    public boolean textureChanged()
    {
        textureChanged = true;
        if(glPanel3d != null)
            glPanel3d.viewChanged = true;
        displayListChanged = true;
        return true;
    }

    public boolean dataChanged()
    {
        return textureChanged();
    }

    public void geometryChanged()
    {
    }
    
    public StsTextureTiles getTextureTiles() { return textureTiles; }
    public void setTextureTiles(StsTextureTiles textureTiles) {this.textureTiles = textureTiles; }

    public void addTextureToDeleteList()
    {
        if(textureTiles != null)
            StsTextureList.addTextureToDeleteList(this);
        textureChanged = true;
        glPanel3d.viewChanged = true;
    }

    public void deleteTexturesAndDisplayLists(GL gl)
    {
        if (textureTiles == null) return;
        textureTiles.deleteTextures(gl);
        textureChanged = true;
    }

    /** Defines in which quadrant the +Inline direction is which determines the layout of the window.
	 * possible orientations and angles: 0: East (315-45), 1: North (45-135),
	 *  2: West (135-225), 3: South (225-315).  Angles are measured counterclockwise from
	 *  the +X axis.  The displayAngle is compass angle measured clockwise from North.
	 */
	private void setViewOrientation()
	{
		xAxisReversed = false; // range on project x axis is reversed
		yAxisReversed = false; // range on project y axis is reversed

		float angle = model.getProject().getAngle();
		int orientation = (int)((angle + 45.0f) / 90.0f);
		orientation = orientation % 4;

		switch(orientation)
		{
			case 0: // East
				axesFlipped = false;
				xAxisReversed = false;
				yAxisReversed = false;
				break;
			case 1: // North
				axesFlipped = true;
				xAxisReversed = false;
				yAxisReversed = true;
				break;
			case 2: // West
				axesFlipped = false;
				xAxisReversed = true;
				yAxisReversed = true;
				break;
			case 3: // South
				axesFlipped = true;
				xAxisReversed = true;
				yAxisReversed = false;
				break;
		}

//		axisLabels[0] = glPanel3d.cursor3d.getHorizontalAxisLabel(currentDirNo);
//		axisLabels[1] = glPanel3d.cursor3d.getVerticalAxisLabel(currentDirNo);

		displayAngle = model.getProject().getAngle();
		/*
		  switch(currentDirNo)
		  {
		   case XDIR: // vertical section defined by x value in y direction

		 if(yAxisReversed)
			flipAxisDirection(0);
		 else
		  displayAngle = (displayAngle + 180f)%360f;
		 break;
		   case YDIR: // vertical section defined by y value in x direction
		 if(xAxisReversed)
		 {
		  flipAxisDirection(0);
		  displayAngle = (displayAngle + 270f)%360f;
		 }
		 else
		  displayAngle = (displayAngle + 90f)%360f;

		 break;
		   case ZDIR: // horizontal section
		 if(yAxisReversed) flipAxisDirection(1);
		 if(xAxisReversed) flipAxisDirection(0);
		 if(axesFlipped)
		 {
		  flipAxes(totalAxisRanges);
		  flipAxes(axisRanges);
//                    flipAxes(axisDisplayRanges);
		  flipLabels(axisLabels);
		 }
		 displayAngle = (displayAngle+90f - orientation*90f)%360f;
		  }
		 */
	}

	private String getStsView2dLegend()
	{
		int compassAngle = Math.round((450f - displayAngle) % 360);
		String compassDirection = get8PointCompassDirection(displayAngle);
		return new String(compassDirection + " (" + compassAngle + ")");
	}

	private String[] getClbView2dLegend()
	{
		String[] directionInfo;
		/*
		  if(currentDirNo == ZDIR)
		  {
		   directionInfo = new String[1];
		   directionInfo[0] = get8PointCompassDirection(displayAngle);
		  }
		  else
		  {
		   directionInfo = new String[2];
		   float leftAngle = (displayAngle + 90f)%360f;
		   directionInfo[0] = get8PointCompassDirection(leftAngle);
		   float riteAngle = (displayAngle + 270f)%360f;
		   directionInfo[1] = get8PointCompassDirection(riteAngle);
		  }
		  return directionInfo;
		 */
		return null;
	}

	private String get8PointCompassDirection(float angle)
	{
		if(angle < 0f)
		{
			StsException.systemError("Angle passed to get9PointCompassDirection() is < 0");
			return compassDirections[0];
		}
		int octant = ((Math.round(angle + 22.5f)) % 360) / 45;
		return compassDirections[octant];
	}

	private void flipAxisDirection(int index)
	{
		flipAxisDirection(totalAxisRanges[index]);
		flipAxisDirection(axisRanges[index]);
//        flipAxisDirection(axisDisplayRanges[index]);
	}

	private void flipAxisDirection(float[] axisRanges)
	{
		float temp = axisRanges[0];
		axisRanges[0] = axisRanges[1];
		axisRanges[1] = temp;
	}

	private void flipAxes(float[][] axisRanges)
	{
		float[] temp = axisRanges[0];
		axisRanges[0] = axisRanges[1];
		axisRanges[1] = temp;
	}

	private void flipLabels(String[] axisLabels)
	{
		String temp = axisLabels[0];
		axisLabels[0] = axisLabels[1];
		axisLabels[1] = temp;
	}

	/**
	 * Set the crossplot associated with this crossplot view to a new object
	 * @param object the new crossplot object
	 * @return true if successfully changed, false if object sent in was not a crossplot object
	 */
	public boolean setViewObject(Object object)
	{
		if(object instanceof StsPreStackLineSet3d)
		{
			return volumeChanged((StsPreStackLineSet3d)object);
		}
		else
			return false;

	}

	/**
	 * Set the volume associated with this gather view to a new volume object
	 * @param seismicVolume the new volume object
	 * @return true if successfully changed
	 */
	public boolean volumeChanged(StsPreStackLineSet3d seismicVolume)
	{
		volumeChanged = (this.lineSet != seismicVolume);
		if(!volumeChanged)return false;
        setLineSet(seismicVolume);
		initialize();
//        setVolumesAndAxisRanges();
		return volumeChanged;
	}

	public void showPopupMenu(StsMouse mouse)
	{
		lineSet.cvsPropertiesDialog(glPanel3d.window, mouse);
	}

	/*
	 * custom serialization requires versioning to prevent old persisted files from barfing.
	 * if you add/change fields, you need to bump the serialVersionUID and fix the
	 * reader to handle both old & new
	 */
	static final long serialVersionUID = 1l;

    public void rescaleVertical()
    {
        rescaleInchesPerSecond();
        glPanel3d.viewChanged = true;
    }

    public void stretch(StsMouse mouse, int index, float factor)
    {
        super.stretch(mouse, index, factor, VERTICAL);
    }

    public void setDefaultView()
    {
        super.setDefaultView(VERTICAL);
        rescaleInchesPerSecond();
    }
    
    protected void displayProfiles()
    {

        StsSemblanceDisplayProperties semblanceDisplayProperties = lineSet.getSemblanceDisplayProperties();
        currentVelocityProfile = getVelocityProfile(superGather.superGatherRow, superGather.superGatherCol);
        // display picked velocity profile
        if (currentVelocityProfile != null)
        {
            // StsVelocityProfile.debugPrintProfile(this, "displayProfiles", currentVelocityProfile);
            StsPoint[] displayPoints = getDisplayedProfilePoints(currentVelocityProfile);
            //StsVelocityProfile.displayOnSemblance(glPanel3d, StsColor.RED, true, false, semblanceDisplayProperties.showLabels, displayPoints, pickedIndex); //showing labels is redundant since already on semblance display
            StsVelocityProfile.displayOnSemblance(glPanel3d, StsColor.RED, true, true, false, displayPoints, pickedIndex);
        }
    }
}