package com.Sts.MVC.ViewWell;

/**
 * Title:        Sts Well Viewer Description:  Well Model-Viewer Copyright:    Copyright (c) 2001
 * Company:      4D Systems LLC
 * @author T.J.Lasseter
 * @version 1.0
 */

import com.Sts.Actions.Wizards.Color.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;

public class StsVspView extends StsWellTextureView
{
	private StsVsp vsp;
	transient float tracesPerInch = 25.f;
    transient float inchesPerSecond = 3.f;
    /** indicates GL is initialized and ready to go. */
    transient public boolean initialized = false;
    transient double timeShift = 0.0f;
    transient double userDatum = 0.0f;
    transient boolean reDatum = false;
    transient int width = 100;
    /** measured depths corresponding to each sample time */
	double[] mdepths;
    transient double [] [] auxValues = null;
	transient boolean tryAux = false;

    static public final String viewVsp = "Well VSP View";

    public StsVspView()
	{
    }

	public StsVspView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow, StsVsp vsp)
	{
        this.vsp = vsp;
        initializeView(wellViewModel, model, actionManager, nSubWindow);
    }

    public boolean initializeView(StsWellViewModel wellViewModel, StsModel model, StsActionManager actionManager, int nSubWindow)
	{
        try
		{
            this.model = model;
			this.wellViewModel = wellViewModel;
            this.well = wellViewModel.well;

            //vsp.setVspView(this);
			vsp.setWell(well);
			//containerPanel = new StsJPanel();
			//splitPane = new JSplitPane(JSplitPane.HORIZONTAL_SPLIT, false);
			int nTraces = vsp.getNCols();
            int nSamples = vsp.getNSlices();

			int pixelsPerInch = 96;
			try
			{
				pixelsPerInch = Toolkit.getDefaultToolkit().getScreenResolution();
			}
			catch(Exception e)
			{}

			float tracesPerInch = vsp.getTracesPerInch();
			int width = (int)(pixelsPerInch*nTraces/tracesPerInch); // 96 dpi
			if(width < 100)
			{
				tracesPerInch = (float)nTraces / (100.f / pixelsPerInch);
				vsp.setTracesPerInch(tracesPerInch, false);
				width = (int)(pixelsPerInch*nTraces/tracesPerInch);
			}
			width = StsMath.minMax(width, 10, 3000);
			width += 29; // insets, but don't have font stuff yet

            // Add height which should change viewport only, not window height.
			float inchesPerSecond = vsp.getInchesPerSecond();
			float heightInSeconds = (vsp.getZMax() - vsp.getZMin())/1000.0f;
            float heightInInches = heightInSeconds*inchesPerSecond;
            int heightInPixels = (int)(heightInInches * pixelsPerInch); // 96 dpi
            if(heightInPixels < 100)
            {
				heightInInches = 100.f/pixelsPerInch; // height in inches
				inchesPerSecond = heightInInches/heightInSeconds;
                vsp.setInchesPerSecond(inchesPerSecond);
                heightInPixels = (int)(heightInInches*pixelsPerInch);
            }
            heightInPixels = StsMath.minMax(heightInPixels, 10, 3000);
			heightInPixels += 29; // insets, but don't have font stuff yet

            /* VSP View is not currently using the height constraint.....Needs to be added - SAJ */

			glPanel = new StsGLPanel(model, actionManager, width, wellViewModel.displayHeight, this);
			//jbw glPanel.setBorder(BorderFactory.createLoweredBevelBorder());
			gl = glPanel.getGL();
			glu = glPanel.getGLU();

            calculateTimeShift();
            initializeRange();
    //		initializeData();
            initialized = true;
            glPanel.setMinimumSize(new Dimension(100, 100));

            //wellGLPanel.setBorder(BorderFactory.createLoweredBevelBorder());
			StsWellWindowPanel wellWindowPanel = wellViewModel.getWellWindowPanel();
			innerPanel = wellWindowPanel.getNewPanel(wellViewModel.curveTrackWidth, wellViewModel.displayHeight);
            GridBagConstraints gbc = innerPanel.gbc;
			gbc.gridx = nSubWindow;
			gbc.gridy = 0;
			gbc.weightx = 1.0;
			gbc.weighty = 1.0;
			gbc.fill = GridBagConstraints.BOTH;
			gbc.anchor = GridBagConstraints.LAST_LINE_START;
			gbc.ipadx = 0;
			gbc.ipady = 0;
			innerPanel.add(glPanel);
        	glPanel.initAndStartGL(glPanel.getGraphicsDevice());
			gbc.anchor = GridBagConstraints.CENTER;
            return true;
		}
		catch(Exception e)
		{
			StsException.outputException("StsVspView.constructor() failed.", e, StsException.WARNING);
            return false;
        }
	}

    private void calculateTimeShift()
    {
        double[] datumValues = null;
        double[] velocityValues = null;
        double userVelocity = 0.0f;
        double newTimeShift = 0.0f;
        double newUserDatum = 0.0f;

        try
        {
            String datumAttribute = vsp.datumProperties.getDatumAttribute();
            String velocityAttribute = vsp.datumProperties.getVelocityAttribute();
            if (!datumAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.NONE]))
            {
                reDatum = true;
                if(!datumAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.USER_SPECIFIED]))
                {
                    datumValues = vsp.getAttributeArray(datumAttribute);
                    newUserDatum = (float)datumValues[0]; // Told it would be the same on every trace by Bill
                }
                else
                    newUserDatum = vsp.datumProperties.getDatum();
            }
            if (!velocityAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.NONE]))
            {
                if(!velocityAttribute.equals(StsDatumProperties.attributeStrings[StsDatumProperties.USER_SPECIFIED]))
                {
                    velocityValues = vsp.getAttributeArray(velocityAttribute);
                    userVelocity = (float)velocityValues[0]; // Told it would be the same on every trace by Bill
                }
                else
                    userVelocity = vsp.datumProperties.getVelocity();
            }
            // Currently assumes time gathers, need to add code for depth gathers and depth shift.
            if ((reDatum) && (newUserDatum != 0.0f) && (userVelocity > 0.0f))
            {
                newTimeShift = newUserDatum/(userVelocity / 1000.0f);
            }

            // If it changed, re-calculate mDepths and times
            if((newTimeShift != timeShift) || (newUserDatum != userDatum))
            {
                userDatum = newUserDatum;
                timeShift = newTimeShift;
                initializeRange();
            }
        }
        catch (Exception e)
        {
            StsException.outputException("Error accessing datum information", e, StsException.WARNING);
        }
    }

    public void removePanel()
	{
		//glPanel.removeAll()
		StsWellWindowPanel wellWindowPanel = wellViewModel.getWellWindowPanel();
		wellViewModel.removeWellWindowPanel(glPanel);
        wellViewModel.removeView(this);
        wellWindowPanel.removeInner(innerPanel);
		wellWindowPanel.rebuild();
        wellViewModel.rebuild();
	}

	public void rebuild(int nSubWindow)
	{
		GridBagLayout g = new GridBagLayout();
		GridBagConstraints gbc = g.getConstraints(glPanel);
		gbc.gridx = nSubWindow;
		g.setConstraints(glPanel,gbc);
    }

	/** StsView2d */
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


	public void display(GLAutoDrawable component)
	{
        if(glPanel.panelViewChanged)
        {
            initializeView();
            glPanel.panelViewChanged = false;
        }
        if(Main.isGLDebug) StsException.systemDebug(this, "display");
		if(!wellViewModel.isVisible)return;

		if(!displayValues && wellViewModel.cursorPicked && !pixelsSaved)
			savePixels(true);

        // Calculate time/depth shift if user has specified.
        calculateTimeShift();

		if(pixelsSaved)
		{
			gl.glDrawBuffer(GL.GL_BACK);
			gl.glClearColor(0.f, 0.f, 0.f, 1.f);
			gl.glClear(GL.GL_COLOR_BUFFER_BIT | GL.GL_DEPTH_BUFFER_BIT);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glDisable(GL.GL_DEPTH_TEST);

            gl.glViewport(0, 0, getWidth(), getHeight());

            gl.glMatrixMode(GL.GL_PROJECTION);
			gl.glPushMatrix();
			gl.glLoadIdentity();
			glu.gluOrtho2D(0.f, getWidth(), 0.f, getHeight());
            gl.glMatrixMode(GL.GL_MODELVIEW);

			doRestorePixels();

		}
		else
		{
			gl.glDrawBuffer(GL.GL_BACK);
			gl.glClearColor(1.f, 1.f, 1.f, 1.f);
			gl.glClear(GL.GL_COLOR_BUFFER_BIT);
			gl.glDisable(GL.GL_LIGHTING);
			gl.glDisable(GL.GL_DEPTH_TEST);

            // case 0: too small;
            /*
            if (getWidth() < 100)
            {
				setWidth(100);
				return;
            }
            */

            /*** Add support for inches per second for vertical scaling, plumbing in place just need display SAJ*/

			// case 1: traces/inch state ok, but window resized
			if (tracesPerInch == vsp.getTracesPerInch())
			{
			   if (width != getWidth())
			   {
				   //glc.setSize((int)(width-4), getHeight());
				   float pixelsPerInch = 96;
				   try
				   {
					   pixelsPerInch = Toolkit.getDefaultToolkit().getScreenResolution();
				   }
				   catch(Exception e)
				   {}

				   width = getWidth();
				   setInsets(vsp.getDisplayAxis());
					int insetLeft = insets.left;
					width -= insetLeft;

				   float tpi = (float)nTraces / ((width) / pixelsPerInch);
				   width = getWidth();
				   //tracesPerInch = tpi
				   vsp.setTracesPerInch(tpi, false);


				   //System.out.println(this+"width change set both tpi to "+tpi);
			   }
            }
			// case 2, user has set the traces/inch -- resize the window

			if(tracesPerInch != vsp.getTracesPerInch())
			{
				//System.out.println(this+"tpi change tpi "+tracesPerInch+ " vsp "+vsp.getTracesPerInch());
				float pixelsPerInch = 96;
				try
				{
					pixelsPerInch = Toolkit.getDefaultToolkit().getScreenResolution();
				}
				catch(Exception e)
				{}

				int width = (int)(((nTraces) / vsp.getTracesPerInch()) * pixelsPerInch); // 96 dpi
				tracesPerInch = vsp.getTracesPerInch();
				{
					setInsets(vsp.getDisplayAxis());
					int insetLeft = insets.left;
					width += insetLeft;
					//System.out.println("targe width is " + width + " plus inset " + insetLeft + "=  current " + getWidth() + " " + getWidth());
					if((Math.abs(getWidth() - width) > 20) && width >= 100) // 96 dpi
					{
						//System.out.println("setSize "+width);
						glPanel.setSize(new Dimension((int)(width+8), getHeight()));
						glPanel.setPreferredSize(new Dimension((int)(width+8), getHeight()));
						//setPreferredSize(new Dimension(width + 8, getHeight()));
						//setMinimumSize(new Dimension(width + 8, getHeight()));
						wellViewModel.rebuild();
					}

				}

			}

			boolean axesOn = vsp.getDisplayAxis();
            if(axesOn)
            {
                setInsets(true);
                insetViewPort();
            }

			computeProjectionMatrix(gl, glu);

			if(vsp.getDisplayVAR())
				displayTexture();

			displayWiggles = vsp.getWiggleDisplayProperties().getDisplayWiggles();
			if(displayWiggles)
			{
				float w2p =  (float)getWidth() / (float)nTraces ;
				if ( w2p > vsp.getWiggleToPixelRatio() )
				   displayWiggleTraces(gl);
			}

            // Display requested attributes on traces  -- SAJ
			String attributeName = vsp.getWiggleDisplayProperties().getAttributeName();
            if(!attributeName.equals(StsWiggleDisplayProperties.ATTRIBUTE_NONE))
            {
                if(!displayAttributeOnGather(gl, attributeName, StsColor.BLUE))
                    vsp.getWiggleDisplayProperties().setAttributeName(StsWiggleDisplayProperties.ATTRIBUTE_NONE);
            }

            if(axesOn)
            {
                setInsets(false);
                resetViewPort();
                gl.glViewport(0, 0, getWidth(), getHeight());

                gl.glMatrixMode(GL.GL_PROJECTION);
                gl.glPushMatrix();
                gl.glLoadIdentity();
                /* pixel correctness ?*/
                glu.gluOrtho2D( -.375f, getWidth() - .375f, -.375f, getHeight() - .375f);
                //glu.gluOrtho2D(0.f,getWidth(),0.f,getHeight());
                gl.glMatrixMode(GL.GL_MODELVIEW);
                drawAxes(component);
            }
		}

		if(doPixelsSaved)doSavePixels();

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

		//gl.glEnable(GL.GL_COLOR_LOGIC_OP);

		//gl.glLogicOp(GL.GL_XOR);

		drawCursor(gl, 0, getWidth(), getHeight(), false);

		gl.glDisable(GL.GL_COLOR_LOGIC_OP);

		gl.glLogicOp(GL.GL_REPLACE);

		gl.glPopMatrix();

		gl.glEnable(GL.GL_DEPTH_TEST);
		gl.glEnable(GL.GL_LIGHTING);
	}

   protected void displayTexture()
    {
        textureTiles = vsp.getTextureTiles();
        if(textureTiles == null) return;

        //computeProjectionMatrix(gl, glu); // should call only if view changed
        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);

        setGLColorList(gl);
        textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null, nullByte);

        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);

        if (textureTiles.shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
    }
/*
    protected void displayTexture()
    {
        if(textureTiles.shaderChanged()) textureChanged = true;
        if(textureTiles == null)
        {
            if(!initializeTextureTiles(gl))return;
            textureChanged = true;

        }

        //computeProjectionMatrix(gl, glu); // should call only if view changed
        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glEnable(GL.GL_TEXTURE_2D);
        gl.glEnable(GL.GL_BLEND);
        gl.glTexEnvf(GL.GL_TEXTURE_ENV, GL.GL_TEXTURE_ENV_MODE, GL.GL_MODULATE);
        gl.glBlendFunc(GL.GL_SRC_ALPHA, GL.GL_ONE_MINUS_SRC_ALPHA);
        gl.glShadeModel(GL.GL_FLAT);

//		if(vsp.colorListChanged)textureChanged = true;

        if(textureTiles.shaderChanged()) textureChanged = true;

        if(isPixelMode != vsp.getIsPixelMode())
        {
            textureChanged = true;
            isPixelMode = !isPixelMode;
        }

        // try a shader
        if(textureTiles.shader != StsJOGLShader.NONE)
        {
            vsp.setGLColorList(gl, false, textureTiles.shader);
            StsJOGLShader.enableARBShader(gl, textureTiles.shader);
        }
        else
        {
            vsp.setGLColorList(gl, false, 0);
        }
        if(textureChanged)
        {
            textureTiles.deleteTextures(gl);
            textureTiles.displayTiles(this, gl, isPixelMode, vsp.getScaledFloatData(), vsp.getWiggleDisplayProperties().getWiggleReversePolarity());
            textureChanged = false;
        }
        else
            textureTiles.displayTiles(this, gl, isPixelMode, (byte[])null);

        gl.glDisable(GL.GL_TEXTURE_2D);
        gl.glDisable(GL.GL_BLEND);
        gl.glEnable(GL.GL_LIGHTING);

        if(textureTiles.shader != StsJOGLShader.NONE) StsJOGLShader.disableARBShader(gl);
    }
*/
    public byte[] getData()
    {
        return vsp.getData();
    }

    private void displayWiggleTraces(GL gl)
	{
		float mdMax = totalAxisRanges[1][0];
		float mdMin = totalAxisRanges[1][1];
		float height = mdMax - mdMin;
		float horizInc = 1;
		float[] data = vsp.getScaledFloatData();

		// If density is less than defined traces to pixels, display wiggles
		if(((axisRanges[0][1] - axisRanges[0][0]) / horizInc) > getWidth() / getPixelsPerWiggle())
			return;

		float horizScale = (100.f + vsp.getWiggleDisplayProperties().getWiggleOverlapPercent()) / 100.f * horizInc / 2.f;

        StsColor.BLACK.setGLColor(gl);
		gl.glDisable(GL.GL_LIGHTING);
		gl.glLineWidth(0.5f);
		boolean isDrawing = false;
		int n = 0;
		//System.out.println("data size "+data.length+" "+nTraces+" "+nSamples);
		for(int t = 0; t < nTraces; t++)
		{
			drawFilledWiggleTraces(gl, n, nSamples, t, mdepths[0], horizScale, vsp.getWiggleDisplayProperties().getWiggleReversePolarity(), data);
			if(vsp.getWiggleDisplayProperties().getWiggleDrawLine())
			{
				for(int s = 0; s < nSamples; s++, n++)
				{
					if(data[n] == -128)
					{
						if(isDrawing)
						{
							gl.glEnd();
							isDrawing = false;
						}
						continue;
					}
					else if(!isDrawing)
					{
						isDrawing = true;
						gl.glBegin(GL.GL_LINE_STRIP);
					}
					if(vsp.getWiggleDisplayProperties().getWiggleReversePolarity())
						gl.glVertex2d(t + horizScale * -(data[n]), mdepths[s]);
					else
						gl.glVertex2d(t + horizScale * (data[n]), mdepths[s]);
				}
				if(isDrawing)gl.glEnd();
				isDrawing = false;
			}
			else
				n += nSamples;
//			gl.glTranslatef(horizInc, 0.0f, 0.0f);
		}
//		gl.popMatrix();
		gl.glEnable(GL.GL_LIGHTING);
	}

	public void drawFilledWiggleTraces(GL gl, int startN, int nSlices, float x0, double zMin, float horizScale, boolean reversePolarity, float[] data)
	{

		try
		{
			if(!vsp.getWiggleDisplayProperties().hasFill())return;
			StsColor plusColor = vsp.getWiggleDisplayProperties().getWigglePlusColor();
			StsColor minusColor = vsp.getWiggleDisplayProperties().getWiggleMinusColor();

			gl.glDisable(GL.GL_LIGHTING);
			//		boolean isDrawing = false;
			//		double tStart = drawPoints[1][0];
			//		double tInc = drawPoints[1][0] - drawPoints[0][0];
			//		int i = 0;
			float t1 = (float)zMin;
			float a1 = horizScale * (data[startN]);
			if(reversePolarity)a1 = -a1;
			boolean b1 = a1 >= 0;
			if(b1)
				plusColor.setGLColor(gl);
			else
				minusColor.setGLColor(gl);

			gl.glBegin(GL.GL_QUAD_STRIP);
			gl.glVertex2f(x0, t1);
			gl.glVertex2f(x0 + a1, t1);
			int n = 1;
			for(int s = startN + 1; s < startN + nSlices; s++, n++)
			{
				float t0 = t1;
				t1 = (float)mdepths[n];
				float a0 = a1;
				a1 = horizScale * (data[s]);
				if(reversePolarity)a1 = -a1;
				boolean b0 = b1;
				b1 = a1 >= 0;
				if(b0 && b1 || !b0 && !b1)
				{
					gl.glVertex2f(x0, t1);
					gl.glVertex2f(x0 + a1, t1);
				}
				else
				{
					float tm = t0 + a0 * (t1 - t0) / (a0 - a1);
					gl.glVertex2f(x0, tm);
					gl.glVertex2f(x0, tm);
					if(b1)
						plusColor.setGLColor(gl);
					else
						minusColor.setGLColor(gl);
					gl.glVertex2f(x0, tm);
					gl.glVertex2f(x0, tm);
					gl.glVertex2f(x0, t1);
					gl.glVertex2f(x0 + a1, t1);
				}
			}
			gl.glEnd();
			gl.glEnable(GL.GL_LIGHTING);
		}
		catch(Exception e)
		{
			StsException.outputWarningException(this, "drawFilledWiggleTraces", e);
		}
	}

	private int getPixelsPerWiggle()
	{

		//return vsp.getWiggleToPixelRatio();
		return 1;
	}
   public boolean displayAttributeOnGather(GL gl, String attName, StsColor color)
    {
        double x = 0.0;
        int horizInc = 1;
        float horizScale = (100.f + vsp.getWiggleDisplayProperties().getWiggleOverlapPercent()) / 100.f * horizInc / 2.f;

        try
        {
            gl.glDisable(GL.GL_LIGHTING);

            double[] attribute = null;
            if(StsSEGYFormat.isTimeAttribute(attName))              // Time attributes to MDepth
                attribute = vsp.getAttributeArrayInMdepth(attName, StsParameters.TD_TIME);
            else if (StsSEGYFormat.isDistanceAttribute(attName))    // Depth attributes
                attribute = vsp.getAttributeArrayInMdepth(attName, StsParameters.TD_DEPTH);
            else                                                    // Other attributes ---- assuming in mdepth already
                attribute = vsp.getAttributeArrayInMdepth(attName, StsParameters.TD_DEPTH);

            if(attribute == null)
            {
                StsMessageFiles.logMessage("Requested attribute: " + attName + " cannot be accessed for current VSP");
                return true;
            }

            // Check attribute range against axis range
            double min = Double.MAX_VALUE;
            double max = -Double.MAX_VALUE;
            for(int i=0; i<attribute.length; i++)
            {
                if(min > attribute[i]) min = attribute[i];
                if(max < attribute[i]) max = attribute[i];
            }
            // Selected attribute is outside time axis range
            if(((min > totalAxisRanges[1][0]) && (max > totalAxisRanges[1][0])) ||
               ((min < totalAxisRanges[1][1]) && (max < totalAxisRanges[1][1])))
            {
                new StsMessage(wellViewModel.getParentFrame(), StsMessage.WARNING,
                               "Entire attribute range in time and in depth is outside axis range.\n Verify that "
                               + attName + " is in time or depth units?");
                return false;
            }

            // Draw the attributes on the traces.
            color.setGLColor(gl);
            gl.glLineWidth(0.5f);
            gl.glBegin(GL.GL_LINE_STRIP);
            for(int i=0; i<attribute.length; i++)
            {
                // Selected attribute is outside time axis range
                if((attribute[i] > totalAxisRanges[1][0]) || (attribute[i] < totalAxisRanges[1][1]))
                    continue;
                StsGLDraw.drawPoint2d((float) i,(float)attribute[i] ,StsColor.BLUE, gl, 4);
            }
            gl.glEnd();
            gl.glEnable(GL.GL_LIGHTING);
        }
        catch(Exception e)
        {
            StsException.outputException("Error displaying attribute on gather", e, StsException.WARNING);
            return false;
        }
        return true;
    }

    private void drawAxes(GLAutoDrawable component)
    {
        float[] range = new float[2];
        range[0] = -.5f;
        range[1] = nTraces - .5f;

        gl.glDisable(GL.GL_DEPTH_TEST);
        gl.glDisable(GL.GL_LIGHTING);
        gl.glColor3f(1.f, 1.f, 1.f);
        // blank out pseudo-inset area
        gl.glBegin(GL.GL_POLYGON);
        gl.glVertex2i(0, getHeight());
        gl.glVertex2i(getWidth(), getHeight());
        gl.glVertex2i(getWidth(), getHeight() - insets.top);
        gl.glVertex2i(0, getHeight() - insets.top);
        gl.glEnd();
        gl.glColor3f(0.f, 0.f, 0.f);
        double [] relev = null;
        if (!tryAux)
        {
            try {
                relev = vsp.getAttributeArray("RELEV");
            } catch (Exception e) {}
        }
        tryAux = true;
        if (relev != null)
        {
            auxValues = new double[2][];
            double[] ix = new double[nTraces];
            boolean valid = false;
            for(int i = 0; i < nTraces; i++)
            {
                ix[i] = i;
                if (relev.length > i && relev[i] != 0)
                    valid = true;
            }
            auxValues[0] = ix;
            auxValues[1] = relev;
            if (!valid) auxValues = null;
        }
        drawHorizontalAxis(vsp.getName(), auxValues == null ? "Trace" : "TVD", range, true, true, gl, true, auxValues);
        //range[1] = (float)wellViewModel.getWindowMdepthMin(getHeight()) + (float)(wellViewModel.zoomLevel.unitsPerPixel * insets.top);
        range[0] = (float)wellViewModel.getWindowMdepthMax();
        range[1] = (float)wellViewModel.getWindowMdepthMin(this);
        drawVerticalAxis("MDepth", range, true, gl);

    }

    protected void displayValues(GL gl)
    {
        StsMousePoint mousePoint = glPanel.getMousePoint();
        int mouseX = mousePoint.x;
        int mouseY = wellViewModel.cursorY;
        int width = getWidth();
        int height = getHeight();
        double depth;
        String depthLabel;
        if(wellViewModel.cursorPicked)
        {
            mouseY = wellViewModel.getCursorY();
            depth = wellViewModel.getMdepthFromMouseY(mouseY, glPanel);
            depthLabel = wellViewModel.getMdepthStringFromGLCursorY();
        }
        else
        {
            depth = wellViewModel.getMdepthFromMouseY(mouseY, glPanel);
            depthLabel = wellViewModel.getMdepthStringFromGLCursorY();
            mouseY = height - mouseY;

        }

        int locationY = StsMath.minMax(mouseY, 20, Math.max(0, height - (int)(5.f * wellViewModel.fontHeight)));

        double fraction = (double)(mouseX - insets.left) / (width - insets.left);
        int value = getValueFromPanelXFraction(fraction);
        //String depthLabel = wellViewModel.labelFormatter.format(depth);
        String valueLabel = "Trace: " + wellViewModel.labelFormatter.format(value);
        String mouseXLabel = wellViewModel.labelFormatter.format(mouseX);
        //String mouseYLabel = wellViewModel.labelFormatter.format(glMouseY);
        float timeStep = 0;
        for(int i = 0; i < mdepths.length - 1; i++)
            if(mdepths[i] <= depth && mdepths[i + 1] >= depth)
            {
                timeStep = i;
                break;
            }
        String ampLabel = "AMP: " + vsp.getAmplitudeAt(value, timeStep);

        float w = wellViewModel.fontWidth * depthLabel.length() * 1.8f;
        w = Math.max(w, wellViewModel.fontWidth * valueLabel.length() * 1.8f);
        w = Math.max(w, wellViewModel.fontWidth * ampLabel.length() * 1.8f);
        int mouseGLy = locationY;
        int locationX = StsMath.minMax(mouseX, 0, Math.max(0, width - (int)w));
        gl.glDisable(GL.GL_LIGHTING);
        gl.glColor3f(1.f, 1.f, 1.f);
        gl.glBegin(GL.GL_POLYGON);
        gl.glVertex2f((float)locationX - 1, (float)mouseGLy + 1);
        gl.glVertex2f((float)locationX + w, (float)mouseGLy + 1);
        if(wellViewModel.currentWellView == this)
        {
            gl.glVertex2f((float)locationX + w, (float)mouseGLy + 5.f * wellViewModel.fontHeight);
            gl.glVertex2f((float)locationX - 1, (float)mouseGLy + 5.f * wellViewModel.fontHeight);
        }
        else
        {
            gl.glVertex2f((float)locationX + w, (float)mouseGLy + 1.5f * wellViewModel.fontHeight);
            gl.glVertex2f((float)locationX - 1, (float)mouseGLy + 1.5f * wellViewModel.fontHeight);
        }
        gl.glEnd();

        getForegroundColor().setGLColor(gl);
        //StsColor.BLACK.setGLColor(gl);
        StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy + .5f * wellViewModel.fontHeight, depthLabel);
        if(wellViewModel.currentWellView == this)
        {
            StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy + 2.f * wellViewModel.fontHeight, valueLabel);
            StsGLDraw.fontHelvetica12(gl, locationX, mouseGLy + 3.5 * wellViewModel.fontHeight, ampLabel);
        }
    }

    protected void initializeRange()
	{
		float tMin, tMax, tInc, mdMin, mdMax;

		nTraces = vsp.getNCols();
		nSamples = vsp.getNSlices();

		// Assumes one-way travel time, thus conversions to one-way below.
		if(vsp.getZDomain() == StsParameters.TD_TIME)
		{
			tMin = vsp.getZMin();
			if(tMin != 0)
				tMin = tMin; /// 2.0f;
			tMax = vsp.getZMax(); /// 2.0f;
            if(reDatum)
            {
                tMin =+ (float)timeShift;
                tMax =+ (float)timeShift;
            }
			tInc = vsp.getZInc(); /// 2.0f;
			mdMin = well.getMDepthFromTime(tMin);
			mdMax = well.getMDepthFromTime(tMax);
			mdepths = well.getMDepthsFromTimes(tMin, tInc, nSamples);
		}
		else
		{
			// vsp is in depth
			float zMin = vsp.getZMin();
			float zMax = vsp.getZMax();
			float zInc = vsp.getZInc();
            if(reDatum)
            {
                zMin =+ (float)userDatum;
                zMax =+ (float)userDatum;
            }
			mdMin = well.getMDepthFromDepth(zMin);
			mdMax = well.getMDepthFromDepth(zMax);
			mdepths = well.getMDepthsFromDepths(zMin, zInc, nSamples);
			//           well.getMDepthsFromDepths(
		}

        totalAxisRanges = new float[2][2];
        totalAxisRanges[0][0] = -0.5f;
		totalAxisRanges[0][1] = nTraces - .5f;
		totalAxisRanges[1][0] = mdMax;
		totalAxisRanges[1][1] = mdMin;

        axisRanges = new float[2][2];
        axisRanges[0][0] = -0.5f;
		axisRanges[0][1] = nTraces - 0.5f;
		axisRanges[1][0] = mdMax;
		axisRanges[1][1] = mdMin;
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
		glPanel.setPanelSize(size);
	}

	public boolean viewObjectChanged(Object source, Object object)
	{
		if (!(object instanceof StsVsp)) return false;
        if(object != vsp) return false;
	    textureChanged(); // for e.g polarity
		return true;
	}

	public boolean viewObjectRepaint(Object source, Object object)
	{
	    if (!(object instanceof StsVsp)) return false;
        if(object != vsp) return false;
        viewChangedRepaint();
		return true;
    }

	public void viewChangedRepaint()
	{
		glPanel.repaint();
	}

	public StsVsp getVsp()
	{
		return vsp;
	}

    public String getViewClassname()
	{
		return viewVsp;
	}
    public void resetToOrigin() { }

    public byte getHorizontalAxisType() { return AXIS_TYPE_NONE; }
    public byte getVerticalAxisType() { return AXIS_TYPE_TIME; }
    public boolean getUseShader() { return vsp.getContourColors(); }
    public int getDefaultShader() { return StsJOGLShader.ARB_TLUT_NO_SPECULAR_LIGHTS; }
    protected boolean getIsPixelMode() { return vsp.getIsPixelMode(); }

    public void adjustColorscale()
    {
        StsColorscalePanel colorPanel = new StsColorscalePanel(true, StsColorscalePanel.COLORSCALE);
        colorPanel.setColorscale(vsp.getColorscale());
	    StsActionManager actionManager = wellViewModel.actionManager;
	    actionManager.startAction(StsColorscaleAction.class, new Object[] { colorPanel } );
    }

    protected float compute2dValue(GL gl, double x, double y)
    {
        return StsParameters.nullValue;
    }

    private boolean initializeTextureTiles(GL gl)
    {
        // axisRanges have origin at lower-left (xmin, zmax) and first axis is horizontal, second is vertical up.
        // textureRanges have origin at upper-left (xmin, zmin) and first axis is vertical down and second is horizontal across.
        float[][] textureRanges = new float[2][2];
        textureRanges[0][0] = totalAxisRanges[1][1];
        textureRanges[0][1] = totalAxisRanges[1][0];
        textureRanges[1] = totalAxisRanges[0];
        if(textureTiles != null) deleteTexturesAndDisplayLists(gl);
        textureTiles = StsTextureTiles.constructor(model, this, nTraces, nSamples, true, textureRanges);
        if(textureTiles == null)return false;
        textureChanged = true;
        return true;
    }

    public void setTextureChanged(boolean v)
    {
        textureChanged = v;
    }

    public void drawTextureTileSurface(StsTextureTile tile, GL gl, boolean is3d)
    {
        tile.drawQuadStripSurface2d(gl, axisRanges[0][0], axisRanges[0][1], mdepths);
    }

    public void deleteTexturesAndDisplayLists(GL gl)
    {}

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

    protected void clearShader()
    {
        vsp.seismicColorList.setGLColorList(gl, false, StsJOGLShader.NONE);
    }

    public Class getDisplayableClass() { return StsVsp.class; }

    protected String getValueLabel(double xCoordinate) { return ""; }

    public void setGLColorList(GL gl) { vsp.setGLColorList(gl, false, textureTiles.shader); }
}
