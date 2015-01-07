package com.Sts.MVC.ViewWell;

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;

/**
 * <p>Title: S2S development</p>
 *
 * <p>Description: Integrated seismic to simulation software</p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author not attributable
 * @version c51c
 */
public class StsWellSeismicCurtainView extends StsWellTextureView
{
	private transient StsSeismicCurtain seismicCurtain;
    StsSeismicVolume seismicVolume = null;
	StsSeismicVolumeClass seismicClass = null;

	public StsWellSeismicCurtainView()
	{
	}

    public StsWellSeismicCurtainView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow, StsSeismicCurtain seismicCurtain)
	{
        this.seismicCurtain = seismicCurtain;
        initializeView(wellViewModel, model, actionManager, nSubWindow);
    }

    public boolean initializeView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow)
    {
        try
        {
            this.wellViewModel = wellViewModel;
            this.model = model;

            this.seismicVolume = seismicCurtain.getSeismicVolume();
            this.seismicClass = seismicCurtain.getSeismicVolumeClass();
            seismicVolume.addActionListener(this);

            int width = seismicCurtain.getNCols() * 2;
            glPanel = new StsGLPanel(model, actionManager, width, wellViewModel.displayHeight, this);
            gl = glPanel.getGL();
            glu = glPanel.getGLU();

            // jbw glPanel.setBorder(BorderFactory.createLoweredBevelBorder());
            StsJPanel innerPanel = wellViewModel.wellWindowPanel.innerPanel();
            GridBagConstraints gbc = innerPanel.gbc;
            gbc.gridx = nSubWindow;
            gbc.gridy = 1;
            gbc.fill = GridBagConstraints.BOTH;
            gbc.weightx = 1.0;
            gbc.weighty = 1.0;
            wellViewModel.addToWellWindowPanel(glPanel, gbc);
            gbc.weightx = gbc.weighty = 0;
            return true;
        }
        catch(Exception e)
        {
            StsException.outputWarningException(this, "initializeView", e);
            return false;
        }
    }

    protected void initializeRange()
	{
        nSamples = seismicCurtain.getNSlices();
		nTraces = seismicCurtain.getNCols();

        float tMin = seismicCurtain.getZMin();
		float tMax = seismicCurtain.getZMax();
        StsWell well = wellViewModel.well;
        float mdMin = well.getMDepthFromTime(tMin);
		float mdMax = well.getMDepthFromTime(tMax);
        float mdInc = (mdMax - mdMin)/(nSamples-1); // hack just to display time data in window
        yCoordinates = new double[nSamples];
        double mdepth = mdMin;
        for(int n = 0; n < nSamples; n++, mdepth += mdInc)
            yCoordinates[n] = mdepth;

        StsGridPoint[] arcGridPoints = seismicCurtain.cellGridCrossingPoints;
        totalAxisRanges = new float[2][2];
        totalAxisRanges[0][0] = 0.f;
		totalAxisRanges[0][1] = arcGridPoints[nTraces-1].point.getF();
		totalAxisRanges[1][0] = (float)yCoordinates[0];
		totalAxisRanges[1][1] = (float)yCoordinates[nSamples-1];

        axisRanges = new float[2][2];
        axisRanges[0][0] = 0.f;
		axisRanges[0][1] = totalAxisRanges[0][1];
		axisRanges[1][0] = totalAxisRanges[1][0];
		axisRanges[1][1] = totalAxisRanges[1][1];
	}

    public void removePanel()
	{
		glPanel.removeAll();
		wellViewModel.removeWellWindowPanel(glPanel);
		wellViewModel.removeView(this);
        wellViewModel.rebuild();
	}

	public void rebuild(int nSubWindow)
	{
		GridBagLayout g = new GridBagLayout();
		GridBagConstraints gbc = g.getConstraints(glPanel);
		gbc.gridx = nSubWindow;
		g.setConstraints(glPanel, gbc);
	}

	/** StsView2d */
	public boolean viewObjectChanged(Object source, Object object)
	{
		if(!(object instanceof StsSeismicVolume)) return false;
        if(object == seismicVolume) return false;
        seismicVolume = (StsSeismicVolume)object;
        textureChanged();
		return true;
	}

	/** object being viewed is changed. Repaint this view if affected.
	 * Implement as needed in concrete subclasses.
	 */
	public boolean viewObjectRepaint(Object source, Object object)
	{
		if(!(object instanceof StsSeismicVolume)) return false;
        if(object == seismicVolume) return false;
        viewChangedRepaint();
		return true;
    }

	/** view2d */
	public void doInitialize()
	{

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
	public void setInsets(boolean axisOn)
	{
		if(axisOn)
		{
			int leftInset = halfWidth + majorTickLength + 2 * verticalFontDimension.width + 2 * fontLineSpace;
			;
			int bottomInset = 0;
			int topInset = halfWidth + majorTickLength + 4 * horizontalFontDimension.height + 3 * fontLineSpace;
			int rightInset = 0;
			insets = new Insets(topInset, leftInset, bottomInset, rightInset);
		}
		else
		{
			insets = new Insets(0, 0, 0, 0);
		}
	}

	public void computeProjectionMatrix()
	{

	}

	public void reshape(GLAutoDrawable component, int x, int y, int width, int height)
	{

	}

    protected StsColor getGridColor()
	{
		return StsColor.BLACK;
	}

	public void setDisplayPanelSize(Dimension size)
	{
		glPanel.setPreferredSize(size);
	}

	public void viewChangedRepaint()
	{
		glPanel.repaint();
	}

	public StsSeismicCurtain getSeismicCurtain()
	{
		return seismicCurtain;
	}

    public void resetToOrigin() {}

    public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }
    public byte getVerticalAxisType() { return AXIS_TYPE_NONE; }
    public boolean getUseShader() { return seismicClass.getContourColors(); }
    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
    protected boolean getIsPixelMode() { return seismicVolume.getIsPixelMode(); }

    protected void setGLColorList(GL gl)
    {
        seismicVolume.seismicColorList.setGLColorList(gl, false, StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS);
    }

    protected void clearShader()
    {
        seismicVolume.seismicColorList.setGLColorList(gl, false, StsJOGLShader.NONE);
    }

    protected byte[] getData()
    {
        return seismicCurtain.getData2D(seismicVolume);
    }

    protected float compute2dValue(GL gl, double x, double y)
    {
        return StsParameters.nullValue;
    }
    protected String getValueLabel(double xCoordinate)
    {
        return "";
    }

    public void adjustColorscale()
    {
        StsColorscalePanel colorPanel = new StsColorscalePanel(true, StsColorscalePanel.COLORSCALE);
        colorPanel.setColorscale(seismicVolume.getColorscale());
	    StsActionManager actionManager = wellViewModel.actionManager;
	    actionManager.startAction(StsColorscaleAction.class, new Object[] { colorPanel } );
    }

    public boolean textureChanged()
    {
        if (seismicVolume == null)
        {
            textureChanged = false;
            return false;
        }
        return textureChanged;
    }

    public Class getDisplayableClass() { return StsSeismicVolume.class; }
}
