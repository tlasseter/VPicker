package com.Sts.MVC.ViewWell;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S development</p>
 * <p/>
 * <p>Description: Integrated seismic to simulation software</p>
 * <p/>
 * <p>Copyright: Copyright (c) 2004</p>
 * <p/>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsTimeLogCurveWellView extends StsWellTextureView implements StsTextureSurfaceFace
{
    /** dataset being isVisible */
    private StsTimeLogCurve dataset;
    /** current position of vertical cursor */
    int cursorX;

    transient private StsJPanel curveNamePanel;
    transient private JLabel currentValueLabel;
    transient StsTimeLogCurveClass datasetClass = null;
    /** cursor is currently being moved */
    transient public boolean verticalCursorPicked = false;
    transient boolean cursorXOffscreen = false;

    static public final String viewCurtain = "Seismic Curtain View";

    public StsTimeLogCurveWellView()
    {
    }

    public StsTimeLogCurveWellView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, StsTimeLogCurve dataset, int nSubWindow)
    {
        this.dataset = dataset;
        initializeView(wellViewModel, model, actionManager, nSubWindow);
    }

    public boolean initializeView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow)
    {
        try
        {
            this.wellViewModel = wellViewModel;
            this.model = model;

            datasetClass = (StsTimeLogCurveClass)dataset.getStsClass();
            initializeRange();

			int width = dataset.getNSurveys() * 4;
            int height = wellViewModel.displayHeight;

			StsWellWindowPanel wellWindowPanel = wellViewModel.getWellWindowPanel();
			wellWindowPanel.setSize(width, height);

            wellWindowPanel.getNewPanel(wellViewModel.curveTrackWidth, wellViewModel.displayHeight);
			innerPanel = wellWindowPanel.innerPanel();
            GridBagConstraints gbc = innerPanel.gbc;
			gbc.weightx = 0.1;
			gbc.weighty = 0;
			gbc.gridy = 0;
			gbc.gridx = nSubWindow;
			gbc.fill = GridBagConstraints.HORIZONTAL;

            curveNameBackPanel = new StsJPanel();
            curveNameBackPanel.setBorder(BorderFactory.createRaisedBevelBorder());
            constructCurveNamePanel();
            curveNameBackPanel.setMaximumSize(new Dimension(1000,40));
			innerPanel.add(curveNameBackPanel, gbc);

			if (Main.useJPanel)
			  glPanel = new StsGLJPanel(model, actionManager, width, height, this);
			else
	  		  glPanel = new StsGLPanel(model, actionManager, width, height, this);

			/*
            gl = glPanel.getGL();
            glu = glPanel.getGLU();
	        */

            // jbw glPanel.setBorder(BorderFactory.createLoweredBevelBorder());
/*
			gbc.gridy = 2;
			gbc.weightx = 1;
			gbc.weighty = 1;
			gbc.fill = GridBagConstraints.BOTH;
*/
			wellWindowPanel.innerPanel.gbc.anchor = GridBagConstraints.LAST_LINE_START;
			wellViewModel.addToWellWindowPanel(glPanel, gbc);
			wellWindowPanel.innerPanel.gbc.anchor = GridBagConstraints.CENTER;

/*
            gbc.gridx = nSubWindow;
            gbc.gridy = 1;
            gbc.fill = GridBagConstraints.BOTH;
            gbc.weightx = 1.0;
            gbc.weighty = 1.0;
			wellWindowPanel.gbc.anchor = GridBagConstraints.LAST_LINE_START;
            wellViewModel.addToWellWindowPanel(glPanel, gbc);
            gbc.weightx = gbc.weighty = 0;
			wellWindowPanel.gbc.anchor = GridBagConstraints.CENTER;
*/
		    setTextureChanged();

            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initializeView", e);
            return false;
        }
    }

    static public StsTimeLogCurveWellView constructor(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, StsTimeLogCurve dataset, int nSubWindow)
    {
        return new StsTimeLogCurveWellView(wellViewModel, model, actionManager, dataset, nSubWindow);
    }

    /** StsView2d */
    public boolean viewObjectChanged(Object source, Object object)
    {
        if (object instanceof StsTimeLogCurve || object instanceof StsTimeLogCurveType)
            setTextureChanged();
        return false;
    }

    public void setTextureChanged() { textureChanged = true; }

    /**
     * object being viewed is changed. Repaint this view if affected.
     * Implement as needed in concrete subclasses.
     */
    public boolean viewObjectRepaint(Object source, Object object)
    {
        if (object instanceof StsTimeLogCurve)
        {
            viewChangedRepaint();
            return true;
        }
        if (object instanceof StsTimeLogCurveType)
        {
            viewChangedRepaint();
            return true;
        }
        return false;
    }

    public void display(GLAutoDrawable component)
    {
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }

        if(glPanel.viewChanged)
        {
            computeProjectionMatrix();
            glPanel.viewChanged = false;
        }
        if (!wellViewModel.isVisible) return;
        //Seismic.SeismicView.init(component);
        // init(component);

        if (!displayValues && wellViewModel.cursorPicked && !pixelsSaved)
            savePixels(true);

        // put a pixel ortho2D projection on stack
         gl.glDrawBuffer(GL.GL_BACK);
         gl.glClearColor(1.f, 1.f, 1.f, 1.f);
         gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
         gl.glDisable(GL.GL_LIGHTING);
         gl.glDisable(GL.GL_DEPTH_TEST);

         gl.glViewport(0, 0, getWidth(), getHeight());

         gl.glMatrixMode(GL.GL_PROJECTION);
         gl.glPushMatrix();
         gl.glLoadIdentity();
         glu.gluOrtho2D(0.f, getWidth(), 0.f, getHeight());
         gl.glMatrixMode(GL.GL_MODELVIEW);

        if (pixelsSaved)
        {
            doRestorePixels();
        }
        else
        {
            gl.glDrawBuffer(GL.GL_BACK);
            gl.glClearColor(1.f, 1.f, 1.f, 1.f);
            gl.glClear(GL.GL_COLOR_BUFFER_BIT);
            gl.glDisable(GL.GL_LIGHTING);
            gl.glDisable(GL.GL_DEPTH_TEST);

            gl.glMatrixMode(GL.GL_PROJECTION);
            gl.glPushMatrix();
            computeProjectionMatrix(gl, glu);

            displayTexture();

            gl.glMatrixMode(GL.GL_PROJECTION);
            gl.glPopMatrix();
            gl.glMatrixMode(GL.GL_MODELVIEW);
            if (doPixelsSaved) doSavePixels();
        }

        if (displayValues)
		{   wellViewModel.push2DOrtho(gl, glu, getWidth(), getHeight());
            displayValues(gl);
			wellViewModel.pop2DOrtho(gl);
        }
        else
        {
            if (wellViewModel.cursorPicked)
			{
				wellViewModel.push2DOrtho(gl, glu, getWidth(), getHeight());
				displayValues(gl);
				wellViewModel.pop2DOrtho(gl);
			}
        }

        //TODO may be easier to do the cursor drawing after the PopMatrix below when we are back in screen coordinates
        //TODO would have to convert screenCoor to globalCoor
		gl.glPopMatrix();
        drawCursor(gl, 0, getWidth(), getHeight(), false);

        drawVerticalCursor(gl);
        // drawVerticalCursor(gl, StsColor.BLACK);
        //gl.glDisable(GL.GL_COLOR_LOGIC_OP);
        //gl.glLogicOp(GL.GL_REPLACE);

        adjustTimeLabel();


        gl.glEnable(GL.GL_DEPTH_TEST);
    }

    public boolean getUseShader() { return datasetClass.getContourColors(); }

    protected boolean getIsPixelMode() { return datasetClass.getIsPixelMode(); }

    protected void setGLColorList(GL gl)
    {
        dataset.setGLColorList(gl, false, StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS);
    }

    protected void clearShader()
    {
        dataset.setGLColorList(gl, false, StsJOGLShader.NONE);
    }

    protected byte[] getData()
    {
        return dataset.getData();
    }


    public String getValueLabel(double xCoordinate)
    {
        long time = (long) (xCoordinate + dataset.bornDate);
        return StsTimeLogCurve.dateTimeFormat.format(time);
    }


    protected void adjustTimeLabel()
    {
        String currentTimeString = dataset.getProjectTimeString();
        currentValueLabel.setText(currentTimeString);
    }

    public void performMouseAction(StsActionManager actionManager, StsMouse mouse)
    {
        try
        {
			int ibut= mouse.getCurrentButton();

            {
                if (ibut == StsMouse.LEFT)
                {
                    int buttonState = mouse.getButtonStateCheckClear(StsMouse.LEFT);
                    if (buttonState == StsMouse.PRESSED)
                    {
                        displayValues = true;
                        if (!wellViewModel.cursorPicked)
                            wellViewModel.checkCursorPicked(wellViewModel.getCursorY(), this); // jbw
                        if (!verticalCursorPicked)
                            checkVerticalCursorPicked(mouse.getMousePoint().x);
                        savePixels(true);
                    }
                    else if (buttonState == StsMouse.DRAGGED)
                    {
                        displayValues = true;
                        if (wellViewModel.cursorPicked)
                            wellViewModel.moveCursor(mouse.getMousePoint().y, this);
                        else if (verticalCursorPicked)
                            moveVerticalCursor(mouse.getMousePoint().x);
                        else
                        {
                            // wellViewModel.setCursorPanel(glPanel);
                            glPanel.setCursor(Cursor.getPredefinedCursor(Cursor.CROSSHAIR_CURSOR));
                            viewChangedRepaint();
                        }
                    }
                    else if (buttonState == StsMouse.RELEASED)
                    {
                        glPanel.setCursor(Cursor.getPredefinedCursor(Cursor.DEFAULT_CURSOR));
                        wellViewModel.cursorPicked = false;
                        verticalCursorPicked = false;
                        displayValues = false;
                        savePixels(false);
                        wellViewModel.display();
                    }
                }

                // If any right mouse action, move view
                else if (ibut == StsMouse.RIGHT)
                {
                    int buttonState = mouse.getButtonStateCheckClear(StsMouse.RIGHT);
                    if ((buttonState == StsMouse.DRAGGED))
                        wellViewModel.moveWindow(mouse);
                }

                // If middle mouse button clicked, terminate any active function.

                // If none active, trigger pop-up menu
                else if (ibut == StsMouse.MIDDLE)
                {
                    mouse.clearButtonState(StsMouse.MIDDLE);
                }
            }
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "performMouseAction", e);
        }
    }

    public void drawVerticalCursor(GL gl)
    {
        wellViewModel.push2DOrtho(gl, glu, getWidth(), getHeight());

        computeCursorX();
		gl.glLineWidth(1.0f);
        getForegroundColor().setGLColor(gl);

        if (cursorXOffscreen)
        {
            gl.glLineStipple(1, StsGraphicParameters.dottedLine);
            gl.glEnable(GL.GL_LINE_STIPPLE);
        }

		gl.glBegin(GL.GL_LINES);
        gl.glVertex2i(cursorX, 0);
        gl.glVertex2i(cursorX, getHeight());
		gl.glEnd();

        if (cursorXOffscreen)
            gl.glDisable(GL.GL_LINE_STIPPLE);

	    wellViewModel.pop2DOrtho(gl);
    }

    private void computeCursorX()
    {
        double currentDeltaT = dataset.getDeltaT();
        int width = getWidth();
        cursorXOffscreen = false;
        double panelDeltaT = axisRanges[0][1] - axisRanges[0][0];
        //System.out.println("Width= " + width + " currentDeltaT= " + currentDeltaT + " panelDeltaT= " + panelDeltaT);
        if (currentDeltaT < 0)
        {
            cursorXOffscreen = true;
            cursorX = 4;
        }
        else if (currentDeltaT > panelDeltaT)
        {
            cursorXOffscreen = true;
            cursorX = width - 4;
        }
        else
        {
            cursorX = (int)(width*currentDeltaT/panelDeltaT);
        }
        //System.out.println("XCursor= " + cursorX);
    }

    public boolean checkVerticalCursorPicked(int mouseX)
    {
        moveVerticalCursor(mouseX);
        verticalCursorPicked = true;
        return verticalCursorPicked;
    }

    public void moveVerticalCursor(int mouseX)
    {
        cursorX = mouseX;
        float deltaT = getXCoordinate(mouseX);
        dataset.setDeltaT(deltaT);
        viewChangedRepaint();
    }

    protected float compute2dValue(GL gl, double xCoordinate, double yCoordinate)
    {
        double fraction = (xCoordinate - axisRanges[0][0]) / (axisRanges[0][1] - axisRanges[0][0]);
        long time = dataset.getTimeFromFraction(fraction);
        float temp = dataset.getValueAtTimeAndMdepth(time, (float) yCoordinate);
        return temp;
    }
/*
    protected void display2dValue(GL gl, int mouseX, int glMouseY)
    {
        double[] pickedPoint = wellTextureGLPanel.getPickedViewPoint();
        long time = (long)pickedPoint[0];
        float mdepth = (float)pickedPoint[1];
        float temp = dataset.getValueAtTimeAndMdepth(time, mdepth);
        if (temp == StsParameters.nullValue) return;
        String tempString = Float.toString(temp);
        StsGLDraw.fontHelvetica12(gl, mouseX, glMouseY, tempString);
    }
*/
    public float getUnitsPerVerticalPixel()
    {
        return (axisRanges[1][1] - axisRanges[1][0]) / getHeight();
    }

    protected void constructCurveNamePanel()
    {
        String curveName = dataset.getName();
        curveNamePanel = new StsJPanel();
        curveNamePanel.setName(curveName);
        curveNamePanel.setBorder(BorderFactory.createRaisedBevelBorder());
        JLabel minValueLabel = new JLabel();
        JLabel curveLabel = new JLabel();
        JLabel maxValueLabel = new JLabel();
        currentValueLabel = new JLabel();

        StsLogCurveType logCurveType = dataset.getLogCurveType();
        Color color = logCurveType.getStsColor().getColor();

        curveNamePanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        curveNamePanel.setBackground(color);

        curveLabel.setText(curveName);
        curveNamePanel.addToRow(curveLabel);

        String currentTimeString = dataset.getProjectTimeString();
        currentValueLabel.setText(currentTimeString);
        curveNamePanel.addEndRow(currentValueLabel);

        String minValueString = dataset.getBornTimeString();
        minValueLabel.setText(minValueString);
        curveNamePanel.addToRow(minValueLabel);

        String maxValueString = dataset.getDeathTimeString();
        maxValueLabel.setText(maxValueString);
        curveNamePanel.addEndRow(maxValueLabel);

        curveNameBackPanel.add(curveNamePanel);
    }

    public void removePanel()
    {

        glPanel.removeAll();
        wellViewModel.removeWellWindowPanel(glPanel);
        wellViewModel.removeView(this);
        // this.delete();
        wellViewModel.rebuild();
    }

    protected void initializeRange()
	{
        float[] mdepthFloats = dataset.getMDepthVectorFloats();
        yCoordinates = StsMath.convertFloatToDoubleArray(mdepthFloats);
        nSamples = yCoordinates.length;
        long[] times = dataset.getTimes();
        nTraces = times.length;
        long t0 = times[0];
        xCoordinates = new double[nTraces];
        for(int n = 0; n < nTraces; n++)
            xCoordinates[n] =  (double)(times[n] - t0);

        totalAxisRanges = new float[2][2];
        totalAxisRanges[0][0] = 0.f;
		totalAxisRanges[0][1] = (float)xCoordinates[nTraces-1];
		totalAxisRanges[1][0] = (float)yCoordinates[0];
		totalAxisRanges[1][1] = (float)yCoordinates[nSamples-1];

        axisRanges = new float[2][2];
        axisRanges[0][0] = 0.f;
		axisRanges[0][1] = (float)xCoordinates[nTraces-1];
		axisRanges[1][0] = (float)yCoordinates[0];
		axisRanges[1][1] = (float)yCoordinates[nSamples-1];
	}

    public void adjustColorscale()
    {
        StsColorscalePanel colorPanel = new StsColorscalePanel(true, StsColorscalePanel.COLORSCALE);
        colorPanel.setColorscale(dataset.getColorscale());
	    StsActionManager actionManager = wellViewModel.actionManager;
	    actionManager.startAction(StsColorscaleAction.class, new Object[] { colorPanel } );
    }

    public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }
    public byte getVerticalAxisType() { return AXIS_TYPE_NONE; }
    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
    public Class getDisplayableClass() { return StsTimeLogCurve.class; }
}
