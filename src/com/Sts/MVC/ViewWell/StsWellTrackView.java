package com.Sts.MVC.ViewWell;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.StsMouse;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;
import com.magician.fonts.*;

import javax.media.opengl.*;
import javax.media.opengl.glu.*;
import javax.swing.*;
import java.awt.*;

public class StsWellTrackView extends StsWellView
{
	String units = MDEPTH;
	/** height of GLPanel */
	protected int displayHeight;
	/** pixels in from left for tick */
	protected static int tickStart = 4;
	/** pixels in from left for first label */
	protected static int labelStart = tickStart;
	/** pixels in from left for label tick start */
	protected static int labelTickStart = 60;
	/** pixels in from left for label tick end */
	protected static int labelTickEnd = labelTickStart + 5;
	/** pixels in from left of non-label tick end */
	protected static int tickEnd = labelTickEnd;
	/** pixels in from left for vertical well line left side */
	protected static int wellLineLeft = tickEnd + 2;
	/** well line width in pixels */
	protected static int wellLineWidth = 10;
	/** pixels in from left for vertical well line right side */
	protected static int wellLineRight = wellLineLeft + wellLineWidth;
	/** total well track width in pixels */
	/** displays units at top of well track panel */
    transient StsComboBoxFieldBean unitsBean;
	/** menuBar for general wellWindow operations */
	transient private StsComboBoxFieldBean operationsMenuComboBox = null;
 	/** name of current operation on this well model */
	transient private String currentOperationType;

    static final boolean runTimer = false;
    static StsTimer timer = null;
    {
        if(runTimer && timer == null) timer = new StsTimer();
    }
	static String MDEPTH = StsWellViewModel.MDEPTH;
	static String DEPTH = StsWellViewModel.DEPTH;
	static String TIME = StsWellViewModel.TIME;
	static String MDEPTH_TIME = StsWellViewModel.MDEPTH_TIME;
	static String DEPTH_TIME = StsWellViewModel.DEPTH_TIME;

	static public final String viewTrack = "Well Track View";

	public StsWellTrackView()
	{
    }

    public boolean initializeView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow)
    {
        try
		{
            this.wellViewModel = wellViewModel;
            this.model = model;
            this.well = wellViewModel.well;
            StsWellWindowPanel wellWindowPanel = wellViewModel.wellWindowPanel;
            setDefaultAction();

			if(!model.getProject().supportsTime())
				unitsBean = new StsComboBoxFieldBean(wellViewModel, "units", "", new String[] {MDEPTH, DEPTH});
			else
				unitsBean = new StsComboBoxFieldBean(wellViewModel, "units", "", new String[] {MDEPTH, DEPTH, TIME, MDEPTH_TIME, DEPTH_TIME});
			unitsBean.setBorder(BorderFactory.createRaisedBevelBorder());

			//create new splitpane

			wellWindowPanel.getRootPanel();
			innerPanel = wellWindowPanel.innerPanel();
            innerPanel.gbc.anchor = GridBagConstraints.CENTER;
            innerPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
            innerPanel.gbc.weightx = 0;
            innerPanel.gbc.weighty = 0;
			innerPanel.gbc.gridy = 0;
			unitsBean.setMaximumSize(new Dimension(1000,20));
            innerPanel.addBeanPanel(unitsBean);

			StsComboBoxFieldBean operationsMenuComboBox = getOperationsComboBox();
			innerPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
			innerPanel.gbc.weightx = 0;
			innerPanel.gbc.weighty = 0;
			innerPanel.gbc.gridy = 1;
			operationsMenuComboBox.setMaximumSize(new Dimension(1000,20));
			innerPanel.add(operationsMenuComboBox);

            if (Main.useJPanel)
			   glPanel = new StsGLJPanel(model, actionManager, wellViewModel.curveTrackWidth, wellViewModel.displayHeight, this);
		    else
			   glPanel = new StsGLPanel(model, actionManager, wellViewModel.curveTrackWidth, wellViewModel.displayHeight, this);
            //glPanel.setBorder(BorderFactory.createLoweredBevelBorder());
            // already set in the glPanel constructor
            // glPanel.setBackgroundColor(model.project.getStsWellPanelColor().getColor());
			innerPanel.gbc.fill = GridBagConstraints.BOTH;
            innerPanel.gbc.weightx = 1;
            innerPanel.gbc.weighty = 1;
			innerPanel.gbc.anchor = GridBagConstraints.LAST_LINE_START;
			innerPanel.gbc.gridy = 2;
			innerPanel.add(glPanel);
			innerPanel.gbc.anchor = GridBagConstraints.CENTER;
			innerPanel.gbc.gridy = 0;

            return true;
        }
		catch(Exception e)
		{
			StsException.outputException("StsWellTrackView.constructor() failed.", e, StsException.WARNING);
            return false;
        }
	}

	public boolean viewObjectChanged(Object source, Object object)
	{
        glPanel.repaint();
		return true;
	}
	/** object being viewed is changed. Repaint this view if affected.
	 * Implement as needed in concrete subclasses.
	 */
	public boolean viewObjectRepaint(Object source, Object object)
	{
        glPanel.repaint();
		return true;
	}

    public void display(GLAutoDrawable component)
	{
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }

        computeProjectionMatrix();

        if(Main.isGLDebug) StsException.systemDebug(this, "display");
		// gl.glClearColor(clearColor);
        // gl.glDrawBuffer(GL.GL_BACK);
        clearToBackground(GL.GL_COLOR_BUFFER_BIT);
        // gl.glClear(GL.GL_COLOR_BUFFER_BIT);
        gl.glDisable(GL.GL_DEPTH_TEST);
     //   displayHeight = glPanel.getGLHeight();
	//	wellViewModel.setHeight(displayHeight);
		drawTicksAndLabels(gl, glu);
		drawVerticalWellLine(gl, glu);
		drawMarkers(gl, glu, wellLineRight);
		drawCursor(gl, 0, getWidth(), getHeight(), true);
        gl.glEnable(GL.GL_DEPTH_TEST);
//        swapBuffers();
	}

    protected int getDisplayZHorizPosition(int left, int right, String zLabel, GLBitmapFont font)
    {
        int labelLength = StsGLDraw.getFontStringLength(horizontalFont, zLabel);
        return right - labelLength;
    }

    public void setAxisRanges()
	{

	}

	public void setDefaultView()
	{

	}
/*
	public void init()
	{
		if(isGLInitialized)return;
        initGL();
		gl.glShadeModel(GL.GL_FLAT);
		gl.glEnable(GL.GL_ALPHA_TEST);
		gl.glAlphaFunc(GL.GL_GREATER, 0.1f);
		initializeFont(gl);
		isGLInitialized = true;
	}
*/
//    public void startGL()
//    {
//        wellTrackDisplayPanel.startGL();
//    }

	public void viewChangedRepaint()
	{
		glPanel.repaint();
	}

	public StsGLPanel getWellGLPanel()
	{
		return glPanel;
	}
	public void savePixels(boolean b)
	{
		glPanel.savePixels(b);
    }

    public boolean initializeDefaultAction()
    {
        defaultAction = new StsWellPanelReadout();
        return true;
    }

    public void setDefaultAction()
    {
        defaultAction = new StsWellPanelReadout();
    }

    static public String getViewClassname()
	{
		return viewTrack;
    }
    public void resetToOrigin() {}

    public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }
    public byte getVerticalAxisType() { return AXIS_TYPE_MEASURED_DEPTH; }

    public void reshape(int x, int y, int width, int height)
	{
        super.reshape(x, y, width, height);
        viewPortChanged(x, y, width, height);
//		displayHeight = height;
	}

	public void drawTicksAndLabels(GL gl, GLU glu)
	{
		float convertedZ = 0.0f;
		String zLabel = "None";
		gl.glLineWidth(1.0f);

		double unitsPerLabel = wellViewModel.zoomLevel.unitsPerLabel;
		double unitsPerTick = wellViewModel.zoomLevel.unitsPerTick;
		double unitsPerPixel = wellViewModel.zoomLevel.unitsPerPixel;
		int pixelsPerTick = wellViewModel.zoomLevel.pixelsPerTick;
		double ticksPerLabel = wellViewModel.zoomLevel.ticksPerLabel;

        windowMdepthMax = wellViewModel.getWindowMdepthMax();
		windowMdepthMin = wellViewModel.getWindowMdepthMin(this);

		double firstLabeledZ = StsMath.intervalRoundUp(windowMdepthMin, unitsPerLabel);
		double firstTickZ = StsMath.intervalRoundUp(windowMdepthMin, unitsPerTick);
		int firstLabeledTick = StsMath.roundOffInteger((firstLabeledZ - firstTickZ) / unitsPerTick);
		int y = getHeight() - (int)((firstTickZ - windowMdepthMin) / unitsPerPixel);
		double z = firstTickZ;
		int nTicks = (y / pixelsPerTick) + 1;
        getForegroundColor().setGLColor(gl);
        //StsColor.BLACK.setGLColor(gl);
		wellViewModel.initializeFont(gl); // jbw make sure its constructed
		// SAJ - Add ability to plot time and td or md
        getForegroundColor().setGLColor(gl);
        for(int n = 0; n < nTicks; n++)
		{
			if(n % ticksPerLabel == firstLabeledTick)
			{
				// StsColor.BLACK.setGLColor(gl);
				String typeString = wellViewModel.getUnits();
				if(typeString == wellViewModel.MDEPTH_TIME)
				{
					gl.glRasterPos2i(labelStart, y + (int)(wellViewModel.fontHeight * .35));
					zLabel = wellViewModel.labelFormatter.format(z);
					wellViewModel.font.drawString(zLabel);

					gl.glRasterPos2i(labelStart, y - (int)(wellViewModel.fontHeight * 1.75));
					convertedZ = wellViewModel.well.getValueFromMDepth((float)z, wellViewModel.TIME);
					zLabel = wellViewModel.labelFormatter.format(convertedZ);
					wellViewModel.font.drawString(zLabel);
				}
				else if(typeString == wellViewModel.DEPTH_TIME)
				{
					gl.glRasterPos2i(labelStart, y + (int)(wellViewModel.fontHeight * .35));
					convertedZ = wellViewModel.well.getValueFromMDepth((float)z, wellViewModel.DEPTH);
					zLabel = wellViewModel.labelFormatter.format(convertedZ);
					wellViewModel.font.drawString(zLabel);

					gl.glRasterPos2i(labelStart, y - (int)(wellViewModel.fontHeight * 1.75));
					convertedZ = wellViewModel.well.getValueFromMDepth((float)z, wellViewModel.TIME);
					zLabel = wellViewModel.labelFormatter.format(convertedZ);
					wellViewModel.font.drawString(zLabel);
				}
				else
				{
					gl.glRasterPos2i(labelStart, y - wellViewModel.fontHeight / 2);
					convertedZ = wellViewModel.well.getValueFromMDepth((float)z, typeString);
					zLabel = wellViewModel.labelFormatter.format(convertedZ);
					wellViewModel.font.drawString(zLabel);
				}
				// StsColor.BLACK.setGLColor(gl);
				gl.glBegin(GL.GL_LINES);
				gl.glVertex2i(labelTickStart, y);
				gl.glVertex2i(labelTickEnd, y);
				gl.glEnd();
			}
			else
			{
				// StsColor.BLACK.setGLColor(gl);
				gl.glBegin(GL.GL_LINES);
				gl.glVertex2i(tickStart, y);
				gl.glVertex2i(tickEnd, y);
				gl.glEnd();
			}
			z += unitsPerTick;
			y -= pixelsPerTick;
		}
	}

	public void drawVerticalWellLine(GL gl, GLU glu)
	{
		int yTop = getHeight();
		int yBot = 0;

		double mdepthMax = wellViewModel.mdepthMax;
		double pixelsPerUnit = wellViewModel.zoomLevel.pixelsPerUnit;
        double mdepthMin = wellViewModel.mdepthMin;
		if(windowMdepthMin < mdepthMin)
		{
			yTop -= (mdepthMin - windowMdepthMin) * pixelsPerUnit;
		}
		if(windowMdepthMax > mdepthMax)
		{
			yBot += (windowMdepthMax - mdepthMax) * pixelsPerUnit;
		}
		int wellLineCenter = wellLineLeft + wellLineWidth / 2;
		well.getStsColor().setGLColor(gl);
		gl.glLineWidth(wellLineWidth);
		gl.glLineStipple(1, StsGraphicParameters.dottedLine);
		gl.glEnable(GL.GL_LINE_STIPPLE);
		gl.glBegin(GL.GL_LINES);
		gl.glVertex2i(wellLineCenter, yTop);
		gl.glVertex2i(wellLineCenter, yBot);
		gl.glEnd();
		gl.glDisable(GL.GL_LINE_STIPPLE);
	}

    public StsComboBoxFieldBean getOperationsComboBox()
	{
        if(operationsMenuComboBox != null) return operationsMenuComboBox;
        operationsMenuComboBox =  new StsComboBoxFieldBean(this, "operationType", null, "viewList");
		return operationsMenuComboBox;
	}

    public String getOperationType()
	{
		return currentOperationType;
	}

	public void setOperationType(String operationType)
	{
        if(wellViewModel == null) return;
        wellViewModel.setOperationType(operationType);
    }

	public String[] getViewList()
	{
        if(wellViewModel == null)
            return new String[] { "none" };
        else
            return wellViewModel.getViewList();
    }

     public void showPopupMenu(StsMouse mouse)
     {
     }
}
