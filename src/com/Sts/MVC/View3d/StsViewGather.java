package com.Sts.MVC.View3d;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: Class defining two-dimensional cursor view. The cursor is a reference to the
 * three planes that are isVisible in the 3D view. An object based on this class wuold result in a 2D
 * view of any one of the three cursor planes.</p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;

abstract public class StsViewGather extends StsViewPreStack implements StsSerializable
{
    boolean isPixelMode;
    /** Has the prestack volume changed and therefore needs to be re-drawn */
	transient public boolean changed = true;
	/** compass angle isVisible in legend at top of 2D view */
	transient public float displayAngle;
	/** legend display type at top of 2D view */
	transient public byte legendType = S2S_LEGEND;
	transient int texture = 0;
	transient boolean deleteTexture = true;
	transient boolean is2dSeismic = false;

	/** String indicating whether offsetAxis is actual offset, or an index from min to max offset, or an index of traces sorted by absolute value */
//	String offsetAxisType = null;
	/** vertical plane at a constant x going y+ and z+ directions */
//	static final int XDIR = StsCursor3d.XDIR;
	/** vertical plane at a constant y going x+ and z+ directions */
//	static final int YDIR = StsCursor3d.YDIR;
	/** horizontal plane at a constant z going x+ and y+ directions  */
//	static final int ZDIR = StsCursor3d.ZDIR;

	static final byte CORE_LABS_LEGEND = 1;
	static final byte S2S_LEGEND = 2;

	static public final String NAME = "Gather_View";

	static final String[] compassDirections = new String[]
		{"E", "NE", "N", "NW", "W", "SW", "S", "SE"};

    static final boolean debug = false;

	static private StsVelocityProfile scratchVelocityProfile = null;

	public StsViewGather()
	{
	}

	/**
	 * StsPreStackView2d constructor
	 * @param glPanel3d the graphics context
	 */
	public StsViewGather(StsGLPanel3d glPanel3d)
	{
		super(glPanel3d);
    }
	
	public StsViewGather(StsGLPanel3d glPanel3d, StsPreStackLineSet lineSet) {
		super(glPanel3d);
		setLineSet(lineSet);
		initialize();
	}
	
	public void init(GLAutoDrawable drawable)
    {
        if (isViewGLInitialized) return;
        super.init(drawable);
        viewChanged();
        //isViewGLInitialized = true;
        StsVelocityCursor.getStsVelocityCursor().addView(this);
        glPanel3d.glc.addMouseMotionListener(StsGatherCursor.createGatherCursorMouseMotionListener(this));
    }
	
    public void initialize()
    {
		isPixelMode = lineSetClass.getIsPixelMode();
		superGather = lineSet.getSuperGather(glPanel3d.window);
//		currentRow = superGather.superGatherRow;
//		currentCol = superGather.superGatherCol;
//		gather.classInitiavlize(currentRow, currentCol);
//		setOffsetAxisType();
        setAxisRanges();
//		StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//		if(toolbar != null) toolbar.zoom();
//		matchLockedWindows();
	}

    public void initializeTransients(StsGLPanel3d glPanel3d)
	{
		super.initializeTransients(glPanel3d);
//		resurrect(this);
//		glPanel3d.mouseMode = StsCursor.ZOOM; 
        setLineSet(StsPreStackLineSetClass.currentProjectPreStackLineSet);
        // set mouse mode to zoom
//		StsMouseActionToolbar toolbar = glPanel3d.window.getMouseActionToolbar();
//		if(toolbar != null)toolbar.zoom();
	}

	public void computeProjectionMatrix()
	{
		if(axisRanges == null) return;
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

	public void reshape(GLAutoDrawable drawable, int x, int y, int w, int h)
	{}

	/**
	 * Display the current cursor plane view
	 * @param component the drawable component
	 * @see StsSeismicVolume#getSeismicVolumeClass()
	 * @see StsSeismicVolumeClass#getCurrentSeismicVolume()
	 * @see #drawHorizontalAxis
	 * @see #drawVerticalAxis
	 */
	public void display(GLAutoDrawable component)
	{
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }
        if(debug)System.out.println("StsViewGather.display() entered in " + toString());
		if(lineSet == null)return;
		if(axisRanges == null) return;

//        if(isCursor3dDragging()) return;

//		setAxisRanges();

        boolean titleOnly = false;
        if(!superGather.filesOk)
            titleOnly = true;

        if (superGather.isComputingGather()) return;

        boolean displayGatherAxis = true;
		if(glPanel3d.viewPortChanged)
		{
            computePixelScaling();
            //resetLimits();
            glPanel3d.viewChanged = true;
            glPanel3d.resetViewPort();
            setInsets(displayGatherAxis);
            if(displayGatherAxis)
                insetViewPort();
            glPanel3d.viewPortChanged = false;
		}

        if(changed)
		{
			computeProjectionMatrix();
//			computeModelViewMatrix();
			changed = false;
		}
//                glPanel3d.viewPortChanged = true;

		if(glPanel3d.viewChanged)
		{
			computeProjectionMatrix();
//			computeModelViewMatrix();
			glPanel3d.viewChanged = false;
		}
		/*
		  if(isPixelMode != lineSet.getIsPixelMode())
		  {
		   deleteTexture(gl);
		   isPixelMode = !isPixelMode;
		   deleteTexture = true;
		  }
		 */
		try
		{
            displayGatherAxis = true;
            setInsets(displayGatherAxis);
			if(displayGatherAxis)
                insetViewPort();

			if(cursorButtonState != StsMouse.CLEARED)
			{
				// if this is window where cursor is being dragged and we have focus, draw foreground cursor.
				// If not the window where cursor is being dragged, but we are displaying cursor here,
				// draw the windows;
				if(isCursorWindow && glPanel3d.hasFocus() || !isCursorWindow)
				{
                    //System.out.println("draw cursor on Gather");
					drawForeground(gl);
					if(cursorButtonState != StsMouse.CLEARED)return;
				}
			}

            //StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
            //glPanel3d.setClearColor(wiggleProperties.getWiggleBackgroundColor());
            clearToBackground(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
            // gl.glDrawBuffer(GL.GL_BACK);
			// glPanel3d.applyClearColor();
			// gl.glClear(GL.GL_DEPTH_BUFFER_BIT | GL.GL_COLOR_BUFFER_BIT);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glDisable(GL.GL_DEPTH_TEST);

			if((superGather != null) && (!titleOnly))
               displayGather(glPanel3d, gl, glu);

			StsCursor3d cursor3d = glPanel3d.getCursor3d();
			int currentDirNo = cursor3d.getCurrentDirNo();
			float currentDirCoordinate = cursor3d.getCurrentDirCoordinate();
//            currentPlaneIndex = cursor3d.getCurrentPlaneIndex();


			if(displayGatherAxis)
			{
				resetViewPort();

				gl.glMatrixMode(GL.GL_PROJECTION);
				gl.glPushMatrix();
				gl.glLoadIdentity();
				glu.gluOrtho2D(0, glPanel3d.getWidth(), 0, glPanel3d.getHeight());
				gl.glMatrixMode(GL.GL_MODELVIEW);

				String axisLabel, titleLabel;
				titleLabel = superGather.getGatherDescription();
				drawHorizontalAxis(titleLabel, axisLabels[0], axisRanges[0], false);
				drawVerticalAxis(axisLabels[1], axisRanges[1], false);

				gl.glMatrixMode(GL.GL_PROJECTION);
				gl.glPopMatrix();
				gl.glMatrixMode(GL.GL_MODELVIEW);

				insetViewPort(); // sets up viewport for picking
			}

			gl.glEnable(GL.GL_DEPTH_TEST);
			gl.glEnable(GL.GL_LIGHTING);
//            getData();
//            glPanel3d.swapBuffers();
		}
		catch(Exception e)
		{
			StsException.outputException("StsViewGather.display() failed.", e, StsException.FATAL);
		}
		finally
		{
			if(debug)System.out.println("StsViewGather.display() exited in " + toString());
		}
	}

    public boolean checkSetOrder(byte order)
	{
		if(superGather.order == order) return true;
        superGather.order = order;
        return false;
    }

	public boolean viewObjectChanged(Object source, Object object)
	{
	    if (object instanceof StsSuperGatherProperties)
	    {
	        this.initialize();
	        return true;
	    }
	    else if(object instanceof StsVelocityProfile)
        {
		    superGather.velocityProfileChanged((StsVelocityProfile)object);
        }
		else if(object instanceof StsAGCProperties)
		{
            superGather.gatherTracesChanged();
            return true;
        }
        else if(object instanceof StsFilterProperties)
        {
            superGather.gatherTracesChanged();
            return true;
        }
        else if(object instanceof StsPreStackLineSet)
        {
            superGather = ((StsPreStackLineSet)object).getSuperGather();
            return true;
        }
        else if(object instanceof StsDatumProperties)
        {
            if(superGather.checkSetDatumShift())
            {
                superGather.gatherTracesChanged();
                return true;
            }
            else
                return false;
        }
        else if(object instanceof StsWiggleDisplayProperties)
        {
            StsWiggleDisplayProperties wiggleProperties = (StsWiggleDisplayProperties)object;
            if(wiggleProperties.isRangeChanged())
            {
                resetRangeWithScale();
                return true;
            }
            if(wiggleProperties.axisTypeChanged)
            {
                setHorizontalAxisRange();
                return true;
            }
            if(wiggleProperties.stretchMuteChanged)
            {
                setStretchMute();
                return true;
            }
            superGather.gatherDataChanged();
            this.viewChanged();
        }
        return false;
    }

	public boolean viewObjectRepaint(Object source, Object object)
	{
        if(object instanceof StsSuperGather)
		{
			glPanel3d.repaint();
			return true;
		}
        else if(object instanceof StsSensor)
		{
			glPanel3d.repaint();
			return true;
		}        
        else if(object instanceof StsPreStackLineSet)
		{
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
        return false;
	}

	protected StsColor getGridColor()
	{
		StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
		return new StsColor(wiggleProperties.getLineColor());
	}

    private void rescaleTracesPerInch()
	{

        float inches = (float)getInsetWidth() / (float)getPixelsPerInch();
        float nTraces = inches * lineSet.getWiggleDisplayProperties().getTracesPerInch();

//        axisRanges[0][0] = totalAxisRanges[0][0];
        axisRanges[0][1] = axisRanges[0][0] + nTraces;
        if(axisRanges[0][1] > totalAxisRanges[0][1])
            axisRanges[0][1] = totalAxisRanges[0][1];
//        glPanel3d.viewChanged = true;
	}

    public void setAxisRanges()
	{
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
		int offsetAxisType = wiggleProperties.getOffsetAxisType();
        System.out.println("StsViewGather offsetAxisType " + offsetAxisType);
        if(debug)System.out.println(toString() + " setting offsetAxisType to " + offsetAxisType);

        if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
		{
			totalAxisRanges = new float[][] { { lineSet.traceOffsetMin, lineSet.traceOffsetMax }, { lineSet.zMax, lineSet.zMin } };
			axisLabels = new String[] { "Trace Offset", "Time" };
		}
		else if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_INDEX)
		{
			totalAxisRanges = new float[][] { { 0.0f, lineSet.maxNTracesPerGather*superGather.nGathersInSuperGather }, { lineSet.zMax, lineSet.zMin } };
			axisLabels = new String[] { "Trace Offset Index", "Time" };
		}
		else
		{
			totalAxisRanges = new float[][] { { 0.0f, lineSet.maxNTracesPerGather*superGather.nGathersInSuperGather }, { lineSet.zMax, lineSet.zMin } };
			axisLabels = new String[] { "Trace Offset Abs Index", "Time" };
		}
		axisRanges = StsMath.copyFloatArray(totalAxisRanges);


        // Need to set scaling properly each time
//        rescaleTracesPerInch();
//        rescaleInchesPerSecond();
//		TJL: commented out temporarily so gatherWindow displayTypes changes work ok
		//setInsets(true);
	}

    private void setStretchMute()
    {
        if(lineSet.velocityModel != null)
        {
            StsVelocityProfile velocityProfile = lineSet.velocityModel.getComputeVelocityProfile(superGather.superGatherRow, superGather.superGatherCol);
            if(velocityProfile != null)
            {
                if(!velocityProfile.isInterpolated())
                {
                    velocityProfile.changeType = StsVelocityProfile.CHANGE_MUTE;
                    superGather.velocityProfileChanged(velocityProfile);
                }
            }
        }
    }

    private void setHorizontalAxisRange()
    {
		StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
		int offsetAxisType = wiggleProperties.getOffsetAxisType();
		if(debug)System.out.println(toString() + " setting offsetAxisType to " + offsetAxisType);

        if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
		{
			totalAxisRanges[0] = new float[] { lineSet.traceOffsetMin, lineSet.traceOffsetMax };
            axisRanges[0] = StsMath.copy(totalAxisRanges[0]);
            axisLabels = new String[] { "Trace Offset", "Time" };
		}
		else if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_INDEX)
		{
			totalAxisRanges[0] = new float[] { 0.0f, lineSet.maxNTracesPerGather*superGather.nGathersInSuperGather };
            axisRanges[0] = StsMath.copy(totalAxisRanges[0]);
            axisLabels = new String[] { "Trace Offset Index", "Time" };
            rescaleTracesPerInch();
        }
		else
		{
			totalAxisRanges[0] = new float[] { 0.0f, lineSet.maxNTracesPerGather*superGather.nGathersInSuperGather };
            axisRanges[0] = StsMath.copy(totalAxisRanges[0]);
            axisLabels = new String[] { "Trace Offset Abs Index", "Time" };
            rescaleTracesPerInch();
        }
        glPanel3d.viewChanged = true;
        computeHorizontalPixelScaling();
//		TJL: commented out temporarily so gatherWindow displayTypes changes work ok
		//setInsets(true);
    }

    public String getHorizontalAxisLabel()
	{
		StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
		int offsetAxisType = wiggleProperties.getOffsetAxisType();

		if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
			return "Trace Offset";
		else if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_INDEX)
			return "Trace Offset Index";
		else
			return "Trace Offset Abs Index";
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
			StsMessageFiles.infoMessage(getHorizontalAxisLabel() + "= " + point.v[0] +
										" " + model.getProject().getZDomainString() + "= " + point.v[1]);
		return cursorPoint;
	}

	protected void displayGather(StsGLPanel3d glPanel3d, GL gl, GLU glu)
	{
        if(superGather == null || superGather.nGathersInSuperGather == 0)return;
		if(!lineSetClass.getDisplayWiggles())return;
        
        superGather.checkInitializeSuperGatherGeometry();
        
        // display must be done on event thread so we won't thread this one out
        // should redesign so gather data is computed on worker thread and we then display
        // by calling viewObjectRepaint as we do in viewSemblance and viewCVS. TJL 1/25/07
        synchronized(superGather)
        {
            superGather.displayWiggleTraces(gl, this);
        }

        StsVelocityProfile velocityProfile = null;
        if(lineSet.velocityModel != null)
        {
            velocityProfile = lineSet.velocityModel.getComputeVelocityProfile(superGather.superGatherRow, superGather.superGatherCol);
            if(velocityProfile != null && lineSet.getWiggleDisplayProperties().getDrawNMOCurves())
                superGather.displayNMOCurve(glPanel3d, axisRanges[0]);
            
            superGather.displayPickStretchMute(glPanel3d.getGL(), StsColor.CYAN);
        }
        
        if (lineSet.getWiggleDisplayProperties().getDrawAttributeCurve()) superGather.displayAttributeOnGather(glPanel3d, StsColor.SADDLEBROWN);
        
        StsVelocityCursor.getStsVelocityCursor().displayCursor(gl, this);
        StsGatherCursor.getStsGatherCursor().displayCursor(gl, this);
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
		StsColor.GRAY.setGLColor(gl);
		gl.glLineWidth((float)3.f);
		gl.glBegin(GL.GL_LINES);
		gl.glVertex2f((float)axisRanges[0][0], xyz[2]);
		gl.glVertex2f((float)axisRanges[0][1], xyz[2]);
		gl.glEnd();

	}
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

	private void setTotalAxisRanges()
	{
		totalAxisRanges = new float[][]
			{
			{0.0f, 100.0f},
			{0.0f, 200.0f}
		};
//		totalAxisRanges = StsMath.copyFloatArray(glPanel3d.cursor3d.getTotalAxisRanges());
		axisRanges = StsMath.copyFloatArray(totalAxisRanges);
	}
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
	 * Output the mouse tracking readout to the information panel on the main screen
	 * @param mouse mouse object
	 */
	public void cursor3dLogReadout2d(StsMouse mouse)
	{
		StsCursorPoint cursorPoint = null;
//		StsCursorPoint cursorPoint = glPanel3d.cursor3d.getCursorPoint2d(mouse, axesFlipped);
//		if(cursorPoint == null) return;
		setCursorXOR(glPanel3d, mouse, cursorPoint.point.v);
//		if(mouse.getButtonState(StsMouse.LEFT) == StsMouse.RELEASED)
//			glPanel3d.cursor3d.logReadout2d(cursorPoint, axesFlipped);
	}

	public StsCursorPoint getCursorPoint(StsMouse mouse)
	{
		return null;
//		return glPanel3d.cursor3d.getCursorPoint2d(mouse, axesFlipped);
	}

	/**
	 * Set the crossplot associated with this crossplot view to a new object
	 * @param object the new crossplot object
	 * @return true if successfully changed, false if object sent in was not a crossplot object
	 */
	public boolean setViewObject(Object object)
	{
		if(object instanceof StsPreStackLineSet)
		{
			return volumeChanged((StsPreStackLineSet)object);
		}
		else
			return false;

	}

	/**
	 * Set the volume associated with this gather view to a new volume object
	 * @param lineSet the new volume object
	 * @return true if successfully changed
	 */
	public boolean volumeChanged(StsPreStackLineSet lineSet)
	{
		changed = (this.lineSet != lineSet);
		if(!changed)return false;
        setLineSet(lineSet);
		initialize();
//        setVolumesAndAxisRanges();
		return changed;
	}

	public void showPopupMenu(StsMouse mouse)
	{
		lineSet.gatherPropertiesDialog(glPanel3d.window);
	}

	/*
	 * custom serialization requires versioning to prevent old persisted files from barfing.
	 * if you add/change fields, you need to bump the serialVersionUID and fix the
	 * reader to handle both old & new
	 */
	static final long serialVersionUID = 1l;

    public void stretch(StsMouse mouse, int index, float factor)
    {
        super.stretch(mouse, index, factor, BOTH);
    }

    public void setDefaultView()
    {
        super.setDefaultView(BOTH);
        resetRangeWithScale();
	}

    public void resetRangeWithScale()
    {
        rescaleInchesPerSecond();
        rescaleTracesPerInch();
        computePixelScaling();
    }

    public byte getHorizontalAxisType()
    {
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
		int offsetAxisType = wiggleProperties.getOffsetAxisType();
		if(offsetAxisType == StsWiggleDisplayProperties.OFFSET_AXIS_VALUE)
            return AXIS_TYPE_DISTANCE;
        else
            return AXIS_TYPE_TRACES; 
    }

    public byte getVerticalAxisType() { return AXIS_TYPE_TIME; }
    
    public void viewChanged()
	{
		super.viewChanged();
		deleteTexture = true;
	}

    public double getCursorValue(double x)
    {
        byte axisType = getHorizontalAxisType();
        if(axisType == AXIS_TYPE_DISTANCE) return x;
        int traceNum = (int) x;
        if (traceNum < 0) traceNum = 0;
        if (traceNum >= superGather.nSuperGatherTraces) traceNum = superGather.nSuperGatherTraces - 1;
        StsWiggleDisplayProperties wiggleProperties = lineSet.getWiggleDisplayProperties();
        String attributeName = wiggleProperties.getAttributeName();
        return superGather.getTraceAttribute(attributeName, traceNum);
    }
}