
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.DB.DBCommand.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.io.*;

/** Constructs a set of StsColors with a defined name */

public class StsSpectrum extends StsMainObject implements StsTreeObjectI
{
    private StsColor[] stsColors;	/** @param colors array of spectrum colors 		*/
    protected StsColor[] keys = null;	/** @param colors array of spectrum colors 		*/
//    protected StsColor currentColor = null;     /** @currentColorNo current color being used 	*/
    protected int currentIndex = 0;

    transient protected StsColorGridFrame colorGridDialog = null;
    transient boolean autoIncrement = true;  /** If dialog sets index, don't increment on next selection */

    static protected StsObjectPanel objectPanel = null;
    static protected Color nullColor = new Color(0,0,0,0);

    static final int PLUS = StsParameters.PLUS;
    static final int MINUS = StsParameters.MINUS;

    static public final StsFieldBean[] displayFields = null;
    static public final StsFieldBean[] propertyFields = new StsFieldBean[]
    {
        new StsStringFieldBean(StsSpectrum.class, "name", true, "Name:")
//        new StsSpectrumFieldBean(StsSpectrum.class, "spectrum")
    };

    /** Default Spectrum constructor */
	public StsSpectrum()
	{
    }

	public StsSpectrum(boolean persistent)
	{
		super(persistent);
    }

    /** Constructs and initializes spectrum
      * @param name spectrum name
      */
/*
    public StsSpectrum(String name, int nColors, StsColor[] inputColors) throws StsException
	{
        this(name);
        StsColor[] stsColors = new StsColor[nColors + 1];
        for(int n = 0; n < nColors; n++)
        	stsColors[n] = new StsColor(inputColors[n]);
        currentIndex = 0;
        stsColors[nColors] = new StsColor(nullColor);
        setStsColors(stsColors);
	}

    public StsSpectrum(String name, Color[] inputColors) throws StsException
    {
		setName(name);
        int nColors = inputColors.length;
        StsColor[] stsColors = new StsColor[nColors+1];
        for(int n = 0; n < nColors; n++)
        	stsColors[n] = new StsColor(inputColors[n]);
        currentIndex = 0;
        stsColors[nColors] = new StsColor(nullColor);
        setStsColors(stsColors);
    }

    public StsSpectrum(String name, StsColor[] inputColors)
    {
		setName(name);
        int nColors = inputColors.length;
        StsColor[] stsColors = new StsColor[nColors + 1];
        for(int n = 0; n < nColors; n++)
        	stsColors[n] = inputColors[n];
        currentIndex = 0;
        stsColors[nColors] = new StsColor(nullColor);
        setStsColors(stsColors);
    }
*/
	static public StsSpectrum constructor(String name, StsColor[] newColors)
	{
		try
		{
			return new StsSpectrum(name, newColors);
		}
		catch(Exception e)
		{
			return null;
		}
	}

	private StsSpectrum(String name, StsColor[] newColors) throws Exception
	{
		super(false);
		setName(name);

		if (newColors == null) throw new Exception();
		keys = null;

		int nColors = newColors.length;
		StsColor[] stsColors = new StsColor[nColors + 1];

		for (int n = 0; n < nColors; n++)
			stsColors[n] = newColors[n];

			// Null Value
		stsColors[nColors] = new StsColor(nullColor);
		setStsColors(stsColors);

	    addToModel();
	}	
	static public StsSpectrum constructor(String name, Color[] newColors)
	{
		try
		{
			return new StsSpectrum(name, newColors);
		}
		catch(Exception e)
		{
			return null;
		}
	}

	private StsSpectrum(String name, Color[] newColors) throws Exception
	{
		super(true);
		setName(name);

		if (newColors == null) throw new Exception();
		keys = null;

		int nColors = newColors.length;
		StsColor[] stsColors = new StsColor[nColors + 1];

		for (int n = 0; n < nColors; n++)
			stsColors[n] = new StsColor(newColors[n]);

			// Null Value
		stsColors[nColors] = new StsColor(nullColor);
		setStsColors(stsColors);

	    addToModel();
	}

	static public StsSpectrum constructor(String name, StsColor[] keys, int nColors)
	{
		try
		{
			return new StsSpectrum(name, keys, nColors);
		}
		catch(Exception e)
		{
			return null;
		}
	}

	private StsSpectrum(String name, StsColor[] keys, int nColors) throws Exception
	{
		super(true);
		if (nColors > 255)
		{
			StsMessageFiles.errorMessage("StsColorWizard.createSpectrum() failed. Maximum number of colors is 255");
			throw new Exception();
		}
		if (keys[0].idx == -1)
			setInterpColors(keys, nColors);
		else
			setInterpColorsWithIdx(keys, nColors);
        setName(name);
        addToModel();
	}

    public StsSpectrum(String name)
	{
		super(false);
		setName(name);
	}
/*
    public StsSpectrum(String name, int nColors, StsColor[] inputColors,  boolean persistent) throws StsException
    {
        this(name, persistent);

        StsColor[] stsColors = new StsColor[nColors+1];
        for(int n = 0; n < nColors; n++)
        	stsColors[n] = new StsColor(inputColors[n]);
        currentIndex = 0;
        stsColors[nColors] = new StsColor(nullColor);
        setStsColors(stsColors);
    }

    public StsSpectrum(String name, boolean persistent)
    {
        super(persistent);
		setName(name);
    }
*/
    public boolean initialize(StsModel model)
    {
        return true;
    }

    /** return the array of Sts colors
      * @return colors
      */
    public StsColor[] getStsColors() { return stsColors; }

    public void setStsColors(StsColor[] stsColors)
    {
        this.stsColors = stsColors;
        if(getIndex() < 0) return;
        if(currentModel == null ) return;
        StsActionManager actionManager = currentModel.mainWindowActionManager;
        if(actionManager == null) return;
        StsChangeCmd cmd = new StsChangeCmd(this, this.stsColors, "stsColors", false);
        currentModel.getCreateTransactionAddCmd("setStsColors", cmd);
    }
    
    public Color[] getColors()
    {
        if(stsColors == null) return null;
        int nColors = stsColors.length - 1;
        Color[] colors = new Color[nColors];
        for( int i=0; i<nColors; i++ )
            colors[i] = stsColors[i].getColor();

        return colors;
    }

    public StsSpectrum getSpectrum() { return this; }
    public int getNColors()
    {
        if(stsColors == null) return 0;
        return stsColors.length - 1;
    }
    public StsColor[] getKeys() { return keys; }

    public void setCurrentColorIndex(int index)
    {
        autoIncrement = false;
        currentIndex = index;
    }

    public int getCurrentColorIndex()
    {
        return currentIndex;
    }

    /** increment the current colorNo
      * if user has selected a color (autoIncrement = false), don't increment
      * this time, but set autoIncrement = true for next time.
      * In the applications, we should delay color requests until the very
      * last moment, so the user as the opportunity to select the one they
      * want.
      */
    public StsColor incrementCurrentColor()
    {
        incrementColorIndex();
        return getCurrentColor();
    }

    public int incrementColorIndex()
    {
        if(autoIncrement)
        {
            int nColors = stsColors.length;

            currentIndex++;
            if(currentIndex >= nColors) currentIndex = 0;

            if(colorGridDialog != null)
            {
                try{ colorGridDialog.buttonSelected(currentIndex); }
                catch(Exception e) { return 0; }
            }
        }
        else
            autoIncrement = true;

        return currentIndex;
    }

/*  Not used anymore:  expensive as each add requires writing all current colors to db.
    Better solution is to add all the colors in one go and then write to db which
    is now accomplished by setStsColors(StsColor[] stsColors)

    public void addColor(StsColor color)
    {
        stsColors = (StsColor[])StsMath.arrayAddElement(stsColors, color);
    }
*/
    /** return the current color
      * @return the current color
      */
    public StsColor getCurrentColor()
    {
        return stsColors[currentIndex];
    }

    public StsColor getColor(int index, StsColor defaultColor)
    {
        if(index < 0 || index > stsColors.length) return defaultColor;
        return stsColors[index];
    }

    public StsColor getTransparentColor()
    {
        return stsColors[stsColors.length-1];
    }

    public StsColor getColor(int index)
    {
        int nColors = stsColors.length;
        if(index < 0) index = -index;
        return stsColors[index%nColors];
    }

    public void setTransparentColor(StsColor stsColor)
    {
        stsColors[stsColors.length-1] = new StsColor(stsColor.getColor());
    }

    public void setTransparentColor(Color color)
    {
        stsColors[stsColors.length-1] = new StsColor(color);
    }
//    public void setCurrentColor(Color color) {
//        StsColor[] stsColors = getStsColors();
//        stsColors[currentIndex].setBeachballColors(color);
//    }
//
    /** display this spectrum */

    public void createColorDialog(int nRows, int nCols)
    {
        colorGridDialog = new StsColorGridFrame(nRows, nCols, this, 0);
    }

    public void display()
    {
        if(colorGridDialog == null)
        {
            int nColors = stsColors.length;
            int nRows = 1 + nColors/32;
            int nCols = nColors/nRows;
    	    colorGridDialog = new StsColorGridFrame(nRows, nCols, this, 0);
        }
        else
            colorGridDialog.setVisible(true);
    }

    /** Close the spectrum display */
    public void closeDisplay()
    {
        if (colorGridDialog != null) colorGridDialog.setVisible(false);
    }

    /** Create the default Rainbow Spectrum */
    public void createRainbowSpectrum()
    {
//        StsColor[] stsColors = new StsColor[65];
        interpolateHue(600, 300, 1.0f, 1.0f, 64);
//        setStsColors(stsColors);
    }

    static public StsSpectrum createRainbowSpectrum(String name, int nColors)
    {
		StsSpectrum spectrum = new StsSpectrum(name);
        spectrum.interpolateHue(300.0f, 0.0f, 1.0f, 1.0f, nColors);
		return spectrum;
    }

    /** Create the default Grayscale Spectrum*/
    public void setGrayScale(int nColors)
    {
        float h = 0.0f;
        float s = 0.0f;
        float b = 0.0f;
        float db;

        try
        {
            StsColor[] stsColors = new StsColor[nColors+1];
            db = 1.0f/(nColors-1);

            for(int n = 0; n < nColors; n++, b = b + db)
                stsColors[n] = new StsColor(Color.getHSBColor(h, s, b));

            stsColors[nColors] = new StsColor(nullColor);
            setStsColors(stsColors);

        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrum.interpolateHue() failed." +
                " couldn't add HSB color - h: " + h + " s: " + s + " b: " + b,
                e, StsException.WARNING);
        }
    }

   /** Sets this color array to a set of standard 32 colors */
    public void setBasic32Colors()
    {
    	float h, s, b, dh;
        int i, n, nc;

       	i = 0;

        StsColor[] stsColors = new StsColor[32];
    /**  full saturation and brightness, hue 0 to 300 deg every 60 deg */

       	s = 1.0f; b = 1.0f; h = 0.0f; nc = 6;
        h = h/360.f;
        dh = 1.0f/(float)nc;
     	for(n = 0; n < nc; n++, i++, h = h + dh)
             stsColors[i] = new StsColor(Color.getHSBColor(h, s, b));

    /**  full saturation and brightness, hue 30 to 330 deg every 60 deg */

       	s = 1.0f; b = 1.0f; h = 30.0f; nc = 6;
        h = h/360.f;
        dh = 1.0f/(float)nc;
     	for(n = 0; n < nc; n++, i++, h = h + dh)
             stsColors[i] = new StsColor(Color.getHSBColor(h, s, b));


    /**  67% saturation and brightness, hue 0 to 300 deg every 60 deg */

       	s = 0.67f; b = 0.67f; h = 0.0f; nc = 6;
        h = h/360.f;
        dh = 1.0f/(float)nc;
     	for(n = 0; n < nc; n++, i++, h = h + dh)
             stsColors[i] = new StsColor(Color.getHSBColor(h, s, b));


    /**  67% saturation and brightness, hue 30 to 330 deg every 60 deg */

       	s = 0.67f; b = 0.67f; h = 30.0f; nc = 6;
        h = h/360.f;
        dh = 1.0f/(float)nc;
     	for(n = 0; n < nc; n++, i++, h = h + dh)
        	stsColors[i] = new StsColor(Color.getHSBColor(h, s, b));

    /**  33% saturation 67% brightness, hue 0 to 300 deg every 60 deg */

       	s = 0.33f; b = 0.67f; h = 0.0f; nc = 6;
        h = h/360.f;
        dh = 1.0f/(float)nc;
     	for(n = 0; n < nc; n++, i++, h = h + dh)
        	stsColors[i] = new StsColor(Color.getHSBColor(h, s, b));

        stsColors[i++] = new StsColor(StsColor.WHITE);
        stsColors[i++] = new StsColor(StsColor.GRAY);

//        stsColors[32] = new StsColor(nullColor);
        setStsColors(stsColors);
    }

   /** Sets this color array to a set of standard 10 colors */
    public void setBasic10Colors()
    {
        StsColor[] stsColors = new StsColor[12];
        stsColors[0] = new StsColor(StsColor.RED);
        stsColors[1] = new StsColor(255, 64, 0);  // red orange
        stsColors[2] = new StsColor(StsColor.ORANGE);
        stsColors[3] = new StsColor(255, 192, 0); // yellow orange
        stsColors[4] = new StsColor(StsColor.YELLOW);
        stsColors[5] = new StsColor(128, 255, 0); // yellow green
        stsColors[6] = new StsColor(StsColor.GREEN);
        stsColors[7] = new StsColor(0, 128, 128); // blue green
        stsColors[8] = new StsColor(StsColor.BLUE);
        stsColors[9] = new StsColor(64, 0, 192); // blue violet
        stsColors[10] = new StsColor(nullColor);
        setStsColors(stsColors);
    }

    public boolean delete()
    {
        stsColors = null;
        super.delete();
        return true;
    }

   /** Interpolates hue specified in degrees with a constant saturation and brightness.
    *  h1 must be > h0
    *  @param h0 starting hue (between 0 and n*360) - you can wrap colors around circle
    *  @param h1 ending hue (between 0 and n*360)
    *  @param s  saturation (between 0.0f and 1.0f)
    *  @param b  brightness (between 0.0f and 1.0f)
    *  @param nColors number of interpolated colors
    */
    public void interpolateHue(float h0, float h1, float s, float b, int nColors)
    {
        float h = -1.0f;
        try
        {
            StsColor[] stsColors = new StsColor[nColors+1];

            h0 /= 360;
            h1 /= 360;

            float dh = (h1 - h0)/(float)nColors;
            h = h0;
            for(int n = 0; n < nColors; n++, h = h + dh)
            {
                float hh = h%1;
                stsColors[n] = new StsColor(Color.getHSBColor(h, s, b));
            }
            stsColors[nColors] = new StsColor(nullColor);
            setStsColors(stsColors);
        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrum.interpolateHue() failed." +
                " couldn't add HSB color - h: " + h + " s: " + s + " b: " + b,
                e, StsException.WARNING);
        }
    }

    public StsSpectrum copySpectrum()
    {
        StsSpectrum spectrumTwo = new StsSpectrum(this.getName() + 1);
        int nColors = stsColors.length;
        StsColor[] newColors = new StsColor[stsColors.length];
        for(int n = 0; n < nColors; n++ )
            newColors[n] = new StsColor(stsColors[n]);
        spectrumTwo.setStsColors(newColors);
        return spectrumTwo;
    }

    public void flipSpectrum()
    {
        int min = 0;
        int max = stsColors.length-1;
        while(min < max)
        {
            StsColor tempColor = stsColors[min];
            stsColors[min] = stsColors[max];
            stsColors[max] = tempColor;
            min++;
            max--;
        }
    }

    /** set colors as an interpolation of n evenly-spaced colors in spectrum */
    public void setInterpColors(StsColor[] inputColors, int nColors) throws StsException
    {
        if(inputColors == null)
        {
            StsException.systemError("StsSpectrum.setInterpColors() failed: inputColors are null.");
            return;
        }

        int nInputColors = inputColors.length;
        if(nInputColors < 2)
        {
            StsException.systemError("StsSpectrum.setInterpColors() failed: number of inputColors < 2.");
            return;
        }

        if(nColors < 2)
        {
            StsException.systemError("StsSpectrum.setInterpColors() failed: number of spectrum colors < 2.");
            return;
        }

        keys = new StsColor[inputColors.length];
        for(int i=0; i<inputColors.length; i++) {
            keys[i] = new StsColor(inputColors[i].getRGBA());
            keys[i].idx = -1;
        }

        StsColor[] stsColors = new StsColor[nColors+1];

        int nColorIntervals = nInputColors - 1;
        float nColorsPerIntervalF = (float)nColors/nColorIntervals;
        float nextColorIndexF = 0.0f;
        int nextColorIndex = 0;
        float[] nextColorRGB = inputColors[0].getRGBA();
        float[] scale = new float[3];
        int nn = 0;
        for(int n = 0; n < nColorIntervals; n++)
        {
            int prevColorIndex = nextColorIndex;
            nextColorIndexF += nColorsPerIntervalF;
            nextColorIndex = Math.round(nextColorIndexF);
            nextColorIndex = Math.min(nextColorIndex, nColors-1);
            float[] prevColorRGB = nextColorRGB;
            nextColorRGB = inputColors[n+1].getRGBA();

            for(nn = prevColorIndex; nn < nextColorIndex; nn++)
            {
                float[] interpolatedColor = StsMath.interpolate(prevColorRGB, nextColorRGB, (nn-prevColorIndex)/nColorsPerIntervalF);
                stsColors[nn] = new StsColor(interpolatedColor);
            }
        }
        for( ; nn < nColors; nn++)
            stsColors[nn] = new StsColor(nextColorRGB);

        currentIndex = 0;
        stsColors[nColors] = new StsColor(nullColor);
        setStsColors(stsColors);
	}

        public void setInterpColorsWithIdx(StsColor[] inputColors, int nColors) throws StsException
        {
        if(inputColors == null)
        {
            StsException.systemError("StsSpectrum.setInterpColors() failed: inputColors are null.");
            return;
        }

        int nInputColors = inputColors.length;
        if(nInputColors < 2)
        {
            StsException.systemError("StsSpectrum.setInterpColors() failed: number of inputColors < 2.");
            return;
        }

        if(inputColors.length < 2)
        {
            StsException.systemError("StsSpectrum.setInterpColors() failed: number of spectrum colors < 2.");
            return;
        }
        int iColors = inputColors.length;
        if(inputColors[inputColors.length-1].idx == nColors)
            iColors--;

        keys = new StsColor[iColors];
        for(int i=0; i<iColors; i++)
        {
            if(inputColors[i].idx > (nColors - 1))
                inputColors[i].idx = nColors - 1;
            keys[i] = new StsColor(inputColors[i].getIRGBA());
        }

        StsColor[] stsColors = new StsColor[nColors+1];

        int nColorIntervals = nInputColors - 1;

        int nextColorIndex = 0;
        float[] nextColorRGB = inputColors[0].getRGBA();
        float[] scale = new float[3];
        int nn = 0;
        for(int n = 0; n < nColorIntervals; n++)
        {
            int nColorsInInterval = inputColors[n+1].idx - inputColors[n].idx;
            if(nColorsInInterval <= 1)
            {
                stsColors[nn] = new StsColor(inputColors[n].getRGBA());
                nextColorIndex = inputColors[n+1].idx;
                nextColorRGB = inputColors[n+1].getRGBA();
                continue;
            }
            else
            {
                int prevColorIndex = nextColorIndex;
                nextColorIndex += nColorsInInterval;
                float[] prevColorRGB = nextColorRGB;
                nextColorRGB = inputColors[n+1].getRGBA();

                for(nn = prevColorIndex; nn < nextColorIndex; nn++)
                {
                    float temp = (float)(nn-prevColorIndex)/(float)nColorsInInterval;
                    float[] interpolatedColor = StsMath.interpolate(prevColorRGB, nextColorRGB, temp);
                    stsColors[nn] = new StsColor(interpolatedColor);
                }
            }
        }
        for( ; nn < nColors; nn++)
            stsColors[nn] = new StsColor(nextColorRGB);

        currentIndex = 0;

        // If the input keys included a key beyond the number of colors, assume it is the null value
        if(inputColors[inputColors.length-1].idx == nColors)
            stsColors[nColors] = new StsColor(inputColors[inputColors.length-1].getRGBA());
        else
            stsColors[nColors] = new StsColor(nullColor);
        setStsColors(stsColors);
	}

    public void addBlackContours(int interval)
    {
        int nColors = stsColors.length;
        for(int n = interval; n < nColors; n += interval)
            stsColors[n] = StsColor.BLACK;
    }

    final public float[][] getArrayRGBA()
    {
        int nColors = stsColors.length;
        float[][] arrayRGBA = new float[4][nColors];
        for(int n = 0; n < nColors; n++)
        {
            StsColor color = stsColors[n];
            arrayRGBA[0][n] = color.red;
            arrayRGBA[1][n] = color.green;
            arrayRGBA[2][n] = color.blue;
            if(color.alpha >= 0.0f)
                arrayRGBA[3][n] = color.alpha;
            else
                arrayRGBA[3][n] = 1.0f;
        }
        return arrayRGBA;
    }


/*
    public StsColorscalePanel spectrumConfigure(JPanel palPanel)
    {
        StsColorscalePanel colorscalePanel = null;

        if(colorscalePanel != null) palPanel.remove(colorscalePanel);
        colorscalePanel = new StsColorscalePanel();
        colorscalePanel.setColors(this.getColors());
        int w = (palPanel.getWidth()*9)/10;
        colorscalePanel.setSize(w, palPanel.getHeight());
        colorscalePanel.setBorder(BorderFactory.createLineBorder(Color.black));
        palPanel.add(colorscalePanel);

        colorscalePanel.repaint();
        return colorscalePanel;
    }
*/
    public StsColorscalePanel spectrumConfigure(JPanel palPanel)
    {
        StsColorscalePanel colorscalePanel = null;

        if(colorscalePanel != null) palPanel.remove(colorscalePanel);
        colorscalePanel = new StsColorscalePanel(this);
        int w = (palPanel.getWidth()*9)/10;
        colorscalePanel.setSize(w, palPanel.getHeight());
        colorscalePanel.setBorder(BorderFactory.createLineBorder(Color.black));
        palPanel.add(colorscalePanel);

        colorscalePanel.repaint();
        return colorscalePanel;
    }

    public void replaceColors(Color[] newColors)
    {
        int nColors = newColors.length + 1;
        StsColor[] stsColors = new StsColor[nColors];

        for(int n = 0; n < nColors; n++)
            stsColors[n] = new StsColor(newColors[n]);
        setStsColors(stsColors);
    }

    public StsFieldBean[] getDisplayFields() { return displayFields; }
    public StsFieldBean[] getPropertyFields() { return propertyFields; }
    public Object[] getChildren() { return new Object[0]; }
    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
//        else objectPanel.setViewObject(this);
        return objectPanel;
    }

    public void treeObjectSelected() { }

    public boolean anyDependencies()
    {
        StsModel model = StsObject.getCurrentModel();
        StsCrossplot[] cp = (StsCrossplot[])model.getCastObjectList(StsCrossplot.class);
        for(int n = 0; n < cp.length; n++)
        {
            if(this == cp[n].getColorscale().getSpectrum())
            {
                StsMessageFiles.infoMessage("Spectrum " + getName() + " by Crossplot " + cp[n].getName());
                return true;
            }
        }

        StsSeismicVolume[] sv = (StsSeismicVolume[])model.getCastObjectList(StsSeismicVolume.class);
        for(int n = 0; n < sv.length; n++)
        {
            if(this == sv[n].getSpectrum())
            {
                StsMessageFiles.infoMessage("Spectrum " + getName() + " by Seismic PostStack3d " + sv[n].getName());
                return true;
            }
        }
        return false;
    }
    public boolean canExport() { return true; }
    public boolean export()
    {
        PrintWriter printWriter;

        StsModel model = StsObject.getCurrentModel();
        String name =  StsStringUtils.cleanString(getName());
        String filename = model.getProject().getRootDirString() + "palette.txt." + name;

        try
        {
            printWriter = new PrintWriter(new FileWriter(filename, false));

            for(int i=0; i<stsColors.length; i++)
            {
                StsColor color = stsColors[i];
                printWriter.println((int)(color.red * 255.0f) + " " +
                                    (int)(color.green * 255.0f) + " " +
                                    (int)(color.blue * 255.0f) + " " +
                                    (int)(color.alpha * 255.0f));
            }
            printWriter.close();

            StsMessageFiles.infoMessage("Spectrum exported to file: " + model.getProject().getRootDirString()
                                    + "palette.txt." + name);
            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSpectrum.export() failed.",
                e, StsException.WARNING);
            return false;
        }
    }

    /** test program */
    public static void main(String[] args)
    {
        try
        {
            //StsColor.initStsColors();
            new StsModel();
            StsSpectrum s = new StsSpectrum(StsSpectrumClass.SPECTRUM_RWB);
            StsColor[] colors = new StsColor[] { StsColor.BLUE, StsColor.WHITE, StsColor.RED };
            s.setInterpColors(colors, 255);
            s.createColorDialog(2, 128);
            s.setTransparentColor(nullColor);
            s.display();
        }
        catch (Exception e) { e.printStackTrace(); }
    }
}



