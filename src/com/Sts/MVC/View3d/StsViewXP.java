package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class to draw and manage the 2D view of cross plotted data</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.1
 */

import com.Sts.Actions.Crossplot.*;
import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;
import com.sun.opengl.util.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.nio.*;

public class StsViewXP extends StsView2d implements StsSerializable
{
    /** Scaling factors when zooming */
	public float[][] rangeScaling;
	/** Active cursor plane XDIR, YDIR or ZDIR */
	public int currentPlaneIndex = -1;
    boolean isPixelMode;

    /** Current crossplot being viewed */
	transient public StsCrossplot crossplot = null;
	/** Model crossplot class which maintains a list of all crossplot objects */
	transient public StsCrossplotClass crossplotClass = null;
	/** Has the crossplot changed and therefore needs to be re-drawn */
	transient public boolean crossplotChanged = true;
	/** The two seismic volumes associated with the current crossplot */
	transient public StsSeismicVolume[] seismicVolumes;

	transient int nRows = 255, nCols = 255;
	transient int nBackgroundRows = 256, nBackgroundCols = 256;
	transient float maxTexCoor = 255.0f / 256.0f;

    transient private ByteBuffer densityData;
	transient private int texture = 0;
	transient boolean deleteTexture = true;
	transient boolean runTimer = false;
	transient StsTimer timer;

	/** Size of color point box. Total size is 6 with 1 pixel border around it. */
	static final int minPointSize = 4;

	static final boolean debug = false;

	static public final String entireVolumeLabel = "Entire PostStack3d";
	static public final String viewNameXP = "XP View";
    static public final String shortViewNameXP = "XP";

    /** Default constructor */
	public StsViewXP()
	{
	}

    /**
     * Crossplot constructor
     * @param glPanel3d the graphics context to use for the plot
     */
    public StsViewXP(StsGLPanel3d glPanel3d)
    {
        super(glPanel3d);
		// System.out.println("StsViewXP: " + toString() + " constructor() glPanel3d: " + glPanel3d.toString());
		initialize();
	}

    static public String getStaticViewName()
    {
        return viewNameXP;
    }

    static public String getStaticShortViewName()
    {
        return shortViewNameXP;
    }
    /**
	 * Initialize the plot by picking up the current crossplot and associated
	 * seismic volumes and then setting up the scaling and axis ranges.
	 */
	public void initialize()
	{
		doInitialize();
	}

	public void doInitialize()
	{
		if(isInitialized)return;

		//super.classInitialize();
		crossplotClass = (StsCrossplotClass)model.getCreateStsClass(StsCrossplot.class);
		crossplot = (StsCrossplot)crossplotClass.getCurrentObject();

		isPixelMode = crossplotClass.getIsPixelMode();
		setAxisRanges();
		// set mouse mode to zoom
        setMouseModeZoom();
        setDefaultView();
//		StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//		if(toolbar != null) toolbar.zoom();
//        labelFormat = model.getProject().getLabelFormat();
	}

	public void initializeTransients(StsGLPanel3d glPanel3d)
	{
		super.initializeTransients(glPanel3d);
		crossplotClass = (StsCrossplotClass)model.getCreateStsClass(StsCrossplot.class);
		crossplot = (StsCrossplot)crossplotClass.getCurrentObject();
		seismicVolumes = crossplot.getVolumes();
		isPixelMode = crossplotClass.getIsPixelMode();
		deleteTexture = true;

//		StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//		if(toolbar != null) toolbar.zoom();
	}

	public void setAxisRanges()
	{
		seismicVolumes = crossplot.getVolumes();
		if(seismicVolumes == null)return;
		totalAxisRanges = new float[2][2];
		axisLabels = new String[2];
		rangeScaling = new float[2][2];
		for(int n = 0; n < 2; n++)
		{
			float dataMin = seismicVolumes[n].dataMin;
			float dataMax = seismicVolumes[n].dataMax;
			totalAxisRanges[n][0] = dataMin;
			totalAxisRanges[n][1] = dataMax;
			rangeScaling[n][0] = dataMin;
			rangeScaling[n][1] = 254.0f / (dataMax - dataMin);
//            rangeScaling[n][1] = (float)nTextureRows/(dataMax - dataMin);
			axisLabels[n] = seismicVolumes[n].getName();
		}
	}

	public void computeProjectionMatrix()
	{
		float[][] axisRanges = getAxisRanges();
		if(axisRanges == null)return;
		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glLoadIdentity();
		glu.gluOrtho2D(axisRanges[0][0], axisRanges[0][1], axisRanges[1][0], axisRanges[1][1]);
//        glu.gluOrtho2D(0, glPanel3d.winRectGL.width, 0, glPanel3d.winRectGL.height);
		gl.glGetDoublev(GL.GL_PROJECTION_MATRIX, glPanel3d.projectionMatrix, 0);
        gl.glMatrixMode(GL.GL_MODELVIEW);
//        gl.glLoadIdentity();
	}

	/**
	 * Get the crossplot that is being isVisible in the crossplot view
	 * @return current crossplot
	 */
	public StsCrossplot getCrossplot()
	{
		return crossplot;
	}

	/**
	 * Get the current X axis label
	 * @return axis label string
	 */
	public String getXAxisLabel()
	{
		return axisLabels[0];
	}

	/**
	 * Get the current Y axis label
	 * @return axis label string
	 */
	public String getYAxisLabel()
	{
		return axisLabels[1];
	}

	/**
	 * Get the two axis ranges in view (factoring in zoom)
	 * @return axis range [2][2](xMin, xMax, yMin, yMax)
	 */
	public float[][] getAxisRanges()
	{
		return axisRanges;
	}

	/**
	 * Get the two total axis ranges, excluding zoom
	 * @return total axis range [2][2](xMin, xMax, yMin, yMax)
	 */
	public float[][] getTotalAxisRanges()
	{
		return totalAxisRanges;
	}

	/**
	 * Set or Reset the view to the default view
	 */
	public void setDefaultView()
	{
		setDefaultAxisRanges();
		computePixelScaling();
	}

	/**
	 * Reset the axis ranges back to full extent
	 */
	public void setDefaultAxisRanges()
	{
		axisRanges = StsMath.copyFloatArray(crossplot.getTotalAxisRanges());
	}

	/**
	 * Is the crossplot viewable
	 * @return true if viewable
	 */
	public boolean isViewable()
	{
		return crossplot != null;
	}

	/**
	 * Set the view as changed. This will force a redraw.
	 */
	public void viewChanged()
	{
        if(glPanel3d != null)
            glPanel3d.viewChanged = true;
	}

	public boolean viewObjectChanged(Object source, Object object)
	{
		if(source instanceof StsCursor3d)
		{
			clearTextureDisplay();
			return true;
		}
		return false;
	}
	public boolean viewObjectRepaint(Object source, Object object)
    {
        if(object instanceof StsCrossplot)
            glPanel3d.repaint();
        return false;
    }

	public void reshape(GLAutoDrawable drawable, int x, int y, int w, int h)
	{}

	/**
	 * Display the current crossplot view
	 * @param component the drawable component
	 */
	public void display(GLAutoDrawable component)
	{
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }
        //System.out.println("Display"+this);
		//try {Thread.sleep(10000);} catch (Exception e) {}
		if(crossplot == null)
		{
			// gl.glDrawBuffer(GL.GL_BACK);
			// glPanel3d.applyClearColor();
			// gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            return;
		}

		float[][] axisRanges = getAxisRanges();
		if(axisRanges == null)return;

        if(crossplotChanged)
		{
			setAxisRanges();
			deleteTexture = true;
//            crossplot.dataDisplayed = false;
			computeProjectionMatrix();
			computeModelViewMatrix();
			crossplotChanged = false;
		}

		if(glPanel3d.viewPortChanged)
		{
			glPanel3d.resetViewPort();
			glPanel3d.viewPortChanged = false;
		}

		if(glPanel3d.viewChanged)
		{
			computeProjectionMatrix();
			computeModelViewMatrix();
			glPanel3d.viewChanged = false;
		}

		if(isPixelMode != crossplot.getIsPixelMode())
		{
			deleteTexture(gl);
			isPixelMode = !isPixelMode;
			deleteTexture = true;
		}

		try
		{
            // gl.glDrawBuffer(GL.GL_BACK);
            // glPanel3d.applyClearColor();
            clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            // gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glDisable(GL.GL_DEPTH_TEST);
            if(crossplotClass == null) // jbw
                crossplotClass = (StsCrossplotClass)model.getCreateStsClass(StsCrossplot.class);
            boolean displayCrossplotAxis = false;
            if(crossplotClass != null)
                displayCrossplotAxis = crossplotClass.getDisplayCrossplotAxis();
            setInsets(displayCrossplotAxis);
            if(displayCrossplotAxis)insetViewPort();

            displayCrossplot();

            StsCursor3d cursor3d = glPanel3d.getCursor3d();
            int currentDirNo = cursor3d.getCurrentDirNo();
            float currentDirCoordinate = cursor3d.getCurrentDirCoordinate();
            //            currentPlaneIndex = cursor3d.getCurrentPlaneIndex();
            displayCrossplotPoints(StsXPolygonAction.getSeismicCrossplotPoints(currentDirNo, currentDirCoordinate));
            displayCrossplotPoints(StsXPolygonAction.getSeismicCrossplotPolygonPoints(currentDirNo, currentDirCoordinate));

            if(displayCrossplotAxis)
            {
                resetViewPort();
                try
                {
					gl.glDisable(GL.GL_LIGHTING);
                    gl.glMatrixMode(GL.GL_PROJECTION);
                    gl.glPushMatrix();
                    gl.glLoadIdentity();
                    //System.out.println("StsVewXP.display() StsViewXP: " + toString() + " glPanel3d: " + glPanel3d.toString());
                    glu.gluOrtho2D(0, glPanel3d.getWidth(), 0, glPanel3d.getHeight());
                    gl.glMatrixMode(GL.GL_MODELVIEW);

                    String axisLabel, titleLabel;
                    axisLabel = getXAxisLabel();
                    // Title needs to be Seismic Format not Crossplot
                    String label = getLabelFormat().format(seismicVolumes[0].getNumFromCoor(currentDirNo, currentDirCoordinate));
                    if(crossplotClass.getDisplayEntireVolume())
                        titleLabel = entireVolumeLabel;
                    else
                        titleLabel = glPanel3d.window.getCursor3d().getTitleLabel(currentDirNo, label);
                    if(crossplotClass.getDisplayOnSubVolumes())
                        titleLabel = titleLabel + " - SubVolume Limited";
                    drawHorizontalAxis(titleLabel, axisLabel, axisRanges[0], crossplotClass.getDisplayGridLinesOnCrossplot());
                    axisLabel = getYAxisLabel();
                    drawVerticalAxis(axisLabel, axisRanges[1], crossplotClass.getDisplayGridLinesOnCrossplot());
                    gl.glMatrixMode(GL.GL_PROJECTION);
                    gl.glPopMatrix();
                    gl.glMatrixMode(GL.GL_MODELVIEW);
                }
                catch(Exception e)
                {
                    StsException.outputException("StsViewXP.display() failed.", e, StsException.WARNING);
                }
                finally
                {

                }
            }
            insetViewPort(); // sets up viewport for picking

            //               gl.glDisable(GL.GL_LIGHTING);

        }
		catch(Exception e)
		{
			StsException.outputException("StsViewXP.display() failed.",
										 e, StsException.FATAL);
		}
		finally
		{
			gl.glEnable(GL.GL_DEPTH_TEST);
			gl.glEnable(GL.GL_LIGHTING);
		}
	}

	// for debugging to look at data
	private float[] getData()
	{
		int width = 255;
		int height = 255;
		FloatBuffer data = BufferUtil.newFloatBuffer(4 * width * height);
		gl.glReadBuffer(GL.GL_BACK);
		gl.glReadPixels(0, 0, width, height, GL.GL_RGBA, GL.GL_FLOAT, data);
		return data.array();
	}

	/**
	 * Key Release event handling
	 *
     *
     * @param mouse mouse object
     * @param e key event
     * @return true if successful
	 */
	public void keyReleased(StsMouse mouse, KeyEvent e)
	{
		int releasedKeyCode = e.getKeyCode();
        super.keyReleased(mouse, e);

        if(mouse.isButtonDown(StsMouse.RIGHT))
		{
			glPanel3d.restoreCursor();
			return;
		}

		if(releasedKeyCode == KeyEvent.VK_UP)
		{
			glPanel3d.checkAddView(StsView3d.class);
			//           glPanel3d.window.getViewSelectToolbar().setButtonEnabled(glPanel3d.getCurrentView(), true);
		}
		else if(releasedKeyCode == KeyEvent.VK_DOWN)
		{
			glPanel3d.checkAddView(StsViewCursor.class);
			//           glPanel3d.window.getViewSelectToolbar().setButtonEnabled(glPanel3d.getCurrentView(), true);
		}
		else if(releasedKeyCode == KeyEvent.VK_LEFT)
			glPanel3d.window.displayPreviousObject(StsCrossplot.class);
		else if(releasedKeyCode == KeyEvent.VK_RIGHT)
			glPanel3d.window.displayNextObject(StsCrossplot.class);
//        else if(releasedKeyCode==KeyEvent.VK_C) compressColorscale();
//        else if(releasedKeyCode==KeyEvent.VK_U) uncompressColorscale();
	}

	// displays the seismic sample point density on a crossPlot
	private void displayCrossplot()
	{
        int dirNo = glPanel3d.getCursor3d().getCurrentDirNo();
        float dirCoordinate = glPanel3d.getCursor3d().getCurrentDirCoordinate();
//        int dirNo = (model.getGLPanel3d(glPanel3d.window)).cursor3d.getCurrentDirNo();
//        float dirCoordinate = (model.getGLPanel3d(glPanel3d.window)).cursor3d.getCurrentDirCoordinate();
//        int planeIndex = (model.getGLPanel3d(glPanel3d.window)).cursor3d.getPlaneIndex(dirNo);



		gl.glDisable(GL.GL_LIGHTING);
		gl.glEnable(GL.GL_TEXTURE_2D);
		gl.glEnable(GL.GL_BLEND);
//        gl.glTexEnvf( GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE );
		gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
		gl.glShadeModel(GL.GL_FLAT);
		gl.glDisable(GL.GL_DEPTH_TEST);

        byte[] byteData;
        if(crossplot.getDisplayDensity())
            byteData = crossplot.getDensityData(dirNo, dirCoordinate, glPanel3d.getCursor3d());
        else
            byteData = crossplot.getAttributeData(dirNo, dirCoordinate, glPanel3d.getCursor3d());
		if(byteData == null)return;

        gl.glCallList(crossplot.getColorListNum(gl, crossplot.getColorscale()));
		densityData = BufferUtil.newByteBuffer(byteData.length);
		densityData.rewind();
		densityData.put(byteData);
		densityData.rewind();

		bindTexture(gl);
		//deleteTexture = true;
		if(deleteTexture)
		{

			if(debug)System.out.println("display subImage changed");
			gl.glTexSubImage2D(GL.GL_TEXTURE_2D, 0, 0, 0, nCols, nRows, GL.GL_COLOR_INDEX, GL.GL_UNSIGNED_BYTE, densityData);
			densityData = null;
			deleteTexture = false;
		}

		if(runTimer)timer.start();
		doDisplay(gl);
		if(runTimer)timer.stopPrint("display crossplot");

		gl.glDisable(GL.GL_TEXTURE_2D);
		gl.glDisable(GL.GL_BLEND);
		gl.glDisable(GL.GL_DEPTH_TEST);
		gl.glEnable(GL.GL_LIGHTING);

		crossplot.drawPolygons(gl, glu, false);
//        displayCrossplotPoints(gl, glu);
	}

	private void displayCrossplotPoints(StsCrossplotPoint[] crossplotPoints)
	{
		if(crossplotPoints == null)return;

		float rowMinF = (axisRanges[1][0] - rangeScaling[1][0]) * rangeScaling[1][1];
		float rowMaxF = (axisRanges[1][1] - rangeScaling[1][0]) * rangeScaling[1][1];
		float colMinF = (axisRanges[0][0] - rangeScaling[0][0]) * rangeScaling[0][1];
		float colMaxF = (axisRanges[0][1] - rangeScaling[0][0]) * rangeScaling[0][1];

		int viewPortHeight = glPanel3d.getHeight();
		float pixelsPerRow = Math.abs(viewPortHeight / (rowMaxF - rowMinF));

		float halfBoxHeight = 0.5f;
		if(pixelsPerRow < minPointSize)
			halfBoxHeight = minPointSize / pixelsPerRow / 2;
		float halfOutlineBoxHeight = halfBoxHeight + 1.0f / pixelsPerRow;

		int viewPortWidth = glPanel3d.getWidth();
		float pixelsPerCol = Math.abs(viewPortWidth / (colMaxF - colMinF));

		float halfBoxWidth = 0.5f;
		if(pixelsPerCol < minPointSize)
			halfBoxWidth = minPointSize / pixelsPerCol / 2;
		float halfOutlineBoxWidth = halfBoxWidth + 1.0f / pixelsPerCol;
		/*
				int height = glPanel3d.getViewPortHeight();
				float pixelsPerRow = (float)height/(rowMaxF - rowMinF + 1);
				int rowMult = 1;
				if(pixelsPerRow < minPointSize)
					rowMult = StsMath.ceiling(minPointSize/pixelsPerRow);

				int width = glPanel3d.getViewPortWidth();
				float pixelsPerCol = (float)width/(colMaxF - colMinF + 1);
				int colMult = 1;
				if(pixelsPerCol < minPointSize)
					colMult = StsMath.ceiling(minPointSize/pixelsPerCol);
		 */
//        GLU glu = glPanel3d.getGLU();

		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPushMatrix();
		gl.glLoadIdentity();
		glu.gluOrtho2D(colMinF, colMaxF, rowMinF, rowMaxF);
		gl.glMatrixMode(GL.GL_MODELVIEW);

		int nPoints = crossplotPoints.length;
		gl.glDisable(GL.GL_LIGHTING);
		gl.glBegin(GL.GL_QUADS);
		for(int n = 0; n < nPoints; n++)
		{
//            StsColor color =  crossplotPoints[n].color;
//            if(color != null) StsColor.setGLColor(gl, color);
			float rowCenter = crossplotPoints[n].crossplotRow;
			float colCenter = crossplotPoints[n].crossplotCol;

			gl.glBegin(GL.GL_QUADS);
			StsColor.BLACK.setGLColor(gl);
			gl.glVertex2f(colCenter - halfOutlineBoxWidth, rowCenter - halfOutlineBoxHeight);
			gl.glVertex2f(colCenter + halfOutlineBoxWidth, rowCenter - halfOutlineBoxHeight);
			gl.glVertex2f(colCenter + halfOutlineBoxWidth, rowCenter + halfOutlineBoxHeight);
			gl.glVertex2f(colCenter - halfOutlineBoxWidth, rowCenter + halfOutlineBoxHeight);

			crossplotPoints[n].stsColor.setGLColor(gl);
			gl.glVertex2f(colCenter - halfBoxWidth, rowCenter - halfBoxHeight);
			gl.glVertex2f(colCenter + halfBoxWidth, rowCenter - halfBoxHeight);
			gl.glVertex2f(colCenter + halfBoxWidth, rowCenter + halfBoxHeight);
			gl.glVertex2f(colCenter - halfBoxWidth, rowCenter + halfBoxHeight);
			gl.glEnd();
			/*
				gl.glVertex2i(col, row);
				gl.glVertex2i(col+colMult, row);
				gl.glVertex2i(col+colMult, row+rowMult);
				gl.glVertex2i(col, row+rowMult);
			 */
		}
		gl.glEnd();
		gl.glEnable(GL.GL_LIGHTING);

		gl.glMatrixMode(GL.GL_PROJECTION);
		gl.glPopMatrix();
		gl.glMatrixMode(GL.GL_MODELVIEW);
	}

	void bindTexture(GL gl)
	{
		if(runTimer)timer.start();
//        if(deleteTexture) deleteTexture(gl);
		if(texture == 0)texture = getTexture(gl);
		gl.glBindTexture(GL.GL_TEXTURE_2D, texture);
		if(runTimer)timer.stopPrint("bindTexture.");
	}

	int getTexture(GL gl)
	{
		gl.glColor4f(0.0f, 0.0f, 0.0f, 0.0f);
		// ByteBuffer background = BufferUtil.newByteBuffer(nBackgroundRows*nBackgroundCols);
		int[] textures = new int[1];
		gl.glGenTextures(1, textures, 0);
		texture = textures[0];
		gl.glPixelStorei(GL.GL_UNPACK_ALIGNMENT, 1);
		gl.glBindTexture(GL.GL_TEXTURE_2D, texture);
		gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_S, GL.GL_CLAMP);
		gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_WRAP_T, GL.GL_CLAMP);
		boolean isPixelMode = crossplot.getIsPixelMode();
		if(!isPixelMode)
		{
			gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_LINEAR);
			gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_LINEAR);
		}
		else
		{
			gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MAG_FILTER, GL.GL_NEAREST);
			gl.glTexParameteri(GL.GL_TEXTURE_2D, GL.GL_TEXTURE_MIN_FILTER, GL.GL_NEAREST);
		}
		gl.glTexImage2D(GL.GL_TEXTURE_2D, 0, GL.GL_RGBA8, nBackgroundCols, nBackgroundRows, 0, GL.GL_COLOR_INDEX, GL.GL_UNSIGNED_BYTE, null); // background
		//background = null;
		return texture;
	}

	void deleteTexture(GL gl)
	{
		if(texture == 0)return;
//        deleteTexture = false;
		gl.glDeleteTextures(1, new int[]
							{texture}, 0);
		texture = 0;
	}

	public void clearTextureDisplay()
	{
		deleteTexture = true;
//        clearDataDisplay();
	}

	public void clearDataDisplay()
	{
		if(crossplot == null)return;
		crossplot.clearData();
	}

	/**
	 * If the crossplot has changed, then clear the colors and display and force redraw
	 * @param changedCrossplot the new crossplot
	 */
	public void checkClearTexture(StsCrossplot changedCrossplot)
	{
		if(crossplot != changedCrossplot)return;
		deleteTexture = true;
	}

	private void doDisplay(GL gl)
	{
		if(runTimer)timer.start();
		gl.glColor4f(1.0f, 1.0f, 1.0f, 1.0f);
		gl.glBegin(GL.GL_QUADS);
		gl.glTexCoord2f(0.0f, 0.0f);
		gl.glVertex2f(totalAxisRanges[0][0], totalAxisRanges[1][0]); // 0,0
		gl.glTexCoord2f(0.0f, maxTexCoor);
		gl.glVertex2d(totalAxisRanges[0][0], totalAxisRanges[1][1]); // 0,1
		gl.glTexCoord2f(maxTexCoor, maxTexCoor);
		gl.glVertex2d(totalAxisRanges[0][1], totalAxisRanges[1][1]); // 1,1
		gl.glTexCoord2f(maxTexCoor, 0.0f);
		gl.glVertex2d(totalAxisRanges[0][1], totalAxisRanges[1][0]); // 1,0
		gl.glEnd();
		if(runTimer)timer.stopPrint("display crossplot");
	}

	/**
	 * Display the popup menu (currently none)
     * @param mouse
     */
	public void showPopupMenu(StsMouse mouse) {}

	/**
	 * Set the crossplot associated with this crossplot view to a new object
	 * @param object the new crossplot object
	 * @return true if successfully changed, false if object sent in was not a crossplot object
	 */
	public boolean setViewObject(Object object)
	{
		if(object instanceof StsCrossplot)
		{
			return crossplotChanged((StsCrossplot)object);
		}
		else
			return false;

	}

	/**
	 * Set the crossplot associated with this crossplot view to a new crossplot object
	 * @param crossplot the new crossplot object
	 * @return true if successfully changed
	 */
	public boolean crossplotChanged(StsCrossplot crossplot)
	{
		if(this.crossplot == crossplot)return false;
		crossplotChanged = true;
		this.crossplot = crossplot;
		setAxisRanges(); // jbw
		//setVolumesAndAxisRanges();
		return crossplotChanged;
	}

	/**
	 * Output the mouse tracking readout to the information panel on the main screen
	 * @param glPanel3d
     * @param mouse mouse object
	 */
	public StsCursorPoint logReadout(StsGLPanel3d glPanel3d, StsMouse mouse)
	{
		int leftButtonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);

		if(leftButtonState == StsMouse.PRESSED || leftButtonState == StsMouse.DRAGGED)
		{
			StsCursorPoint cursorPoint = this.glPanel3d.getCursor3d().getNearestPointInCursorPlane(glPanel3d, mouse);
			StsPoint point = this.glPanel3d.getPointInPlaneAtMouse(mouse);
			String nameX = seismicVolumes[0].getName();
			String nameY = seismicVolumes[1].getName();
			StsMessageFiles.infoMessage(" horizontal Axis (" + nameX + "): " + point.v[0] +
										" vertical Axis (" + nameY + "): " + point.v[1]);
			return cursorPoint;
		}
		return null;
	}
    public void resetToOrigin() { }
	public void removeDisplayableClass(StsObject object)
	{
		if(object.getClass() != StsCrossplot.class)return;
//        Component buttonXP = glPanel3d.window.getToolbarComponentNamed(StsViewSelectToolbar.NAME, StsViewSelectToolbar.CROSSPLOT);
//        buttonXP.setVisible(false);
//        glPanel3d.checkAddView(StsView3d.class);
	}

	/*
	 * custom serialization requires versioning to prevent old persisted files from barfing.
	 * if you add/change fields, you need to bump the serialVersionUID and fix the
	 * reader to handle both old & new
	 */
	static final long serialVersionUID = 1l;

    public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }
    public byte getVerticalAxisType() { return AXIS_TYPE_NONE; }
}
