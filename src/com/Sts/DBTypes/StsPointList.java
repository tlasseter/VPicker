package com.Sts.DBTypes;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.awt.event.*;
import java.io.*;

public class StsPointList extends StsMainObject implements StsCustomSerializable, StsInstance3dDisplayable, StsTreeObjectI, StsViewSelectable
{
    transient StsPoint[] points = null;
//    String name = "pointSet";  // already exists in superclass
    int nPoints = 0;
    int initialSize;
    boolean showAllOn2d = true;
    int minVolume = 999;
    int maxVolume = -999;
    double dataMin = StsParameters.largeDouble;
    double dataMax = -StsParameters.largeDouble;
    byte sizeBy = AMPLITUDE;
    byte colorBy = AMPLITUDE;
    int numberBins = 10;
    String[] volumes = null;
    byte setType = VALUE;
    String spectrumName = "Basic";
    byte pointType = POINT;
    transient protected StsColorscale colorscale;

    /** Histogram of the data distribution */
	public float[] dataHist = new float[255];

    /** Total samples in each of 255 steps */
    transient private int dataCnt[] = new int[255];
    transient private int ttlHistogramSamples = 0;

    transient protected StsModel model = null;
    transient protected String currentVolumeName = "none";
    transient protected int startVol = 0;
    transient protected int endVol = 10000;
    transient protected int currentVol = 0;
    transient protected boolean animated = false;
    static StsObjectPanel objectPanel = null;

    transient public static final byte POINT = 0;
    transient public static final byte SPHERE = 1;
    transient public static final byte DISK = 2;

    transient public static final byte NONE = 0;
    transient public static final byte VOLUME = 1;
    transient public static final byte AMPLITUDE = 2;
    transient public static final byte RED = 3;
    transient public static final byte GREEN = 4;
    transient public static final byte BLUE = 5;
    transient public static final byte YELLOW = 6;

    transient public static final byte VALUE = 0;
    transient public static final byte PERCENTAGE = 1;
    transient public static String[] SETTYPES = {"Data Range", "Percentage"};
//    transient private int size = 2;
    transient private int colorIdx = 0;
    transient StsMethodPick pointPicker;

	transient StsProgressBarDialog progressBarDialog = null;

    static public final String[] SYMBOL_TYPE_STRINGS = new String[] { "Point", "Sphere", "Cylinder"};
    static public final byte[] SYMBOL_TYPES = new byte[] { POINT, SPHERE, DISK };

    static public final String[] SIZE_TYPE_STRINGS = new String[] { "None", "PostStack3d", "Amplitude"};
    static public final String[] COLOR_TYPE_STRINGS = new String[] { "PostStack3d", "Amplitude", "Red", "Green", "Blue", "Yellow" };
    static public final byte[] SIZE_TYPES = new byte[] { NONE, VOLUME, AMPLITUDE };
    static public final byte[] COLOR_TYPES = new byte[] { VOLUME, AMPLITUDE, RED, GREEN, BLUE, YELLOW };

    static public StsComboBoxFieldBean symbolListBean ;
    static public StsComboBoxFieldBean sizeListBean;
    static public StsComboBoxFieldBean colorListBean;
    static public StsIntFieldBean numBinsBean;
    static public StsIntFieldBean startVolumeBean;
    static public StsIntFieldBean endVolumeBean;
    static public StsStringFieldBean volumeNameBean;
    static StsEditableColorscaleFieldBean colorscaleBean;

    static public StsFieldBean[] displayFields = null;

    public StsPointList(boolean persistent)
    {
        super(persistent);
    }

    public StsPointList()
    {
        this(10);
    }

    public StsPointList(int size)
    {
        initialSize = Math.max(size, 2);
        points = new StsPoint[initialSize];
        model = currentModel;
		setName("pointSet");
    }

    public StsPointList(int size, double min, double max, boolean persistant)
    {
        super(persistant);
        initialSize = Math.max(size, 2);
        dataMin = min;
        dataMax = max;
        points = new StsPoint[initialSize];
        model = currentModel;
        initializeColorscale();
        numberBins = getPointListClass().getDefaultNumberBins();
		setName("pointSet");
   }

    public StsPointList(String name, StsPoint[] points, String[] names, boolean persistant)
    {
        super(persistant);
        initialSize = 2;
        setName(name);
        volumes = names;
        setPoints(points);
        model = currentModel;
        initializeColorscale();
        buildHistogram();
        numberBins = getPointListClass().getDefaultNumberBins();
    }

    static public StsPointListClass getPointListClass()
    {
        return (StsPointListClass)currentModel.getCreateStsClass(StsPointList.class);
    }

    public boolean initialize(StsModel model)
    {
        try
        {
            this.model = model;
            startVol = minVolume;
            endVol = maxVolume;
            if(volumes != null) setCurrentVolumeName(volumes[endVol]);
            super.isVisible = isVisible;

            if(colorscale == null) initializeColorscale();
            buildHistogram();
            setDataHistogram();

            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPointList.classInitialize() failed.", e, StsException.WARNING);
            return false;
        }
    }

	private void initializeColorscale()
	{
		spectrumName = getPointListClass().getDefaultSpectrumName();
		colorscale = new StsColorscale("PointList", model.getSpectrum(spectrumName),(float) dataMin, (float) dataMax);
		colorscale.addActionListener(this);
	}

    private void resetColorscale()
    {
        if(model == null)
            return;
        if(colorscale != null)
        {
            if (colorBy == AMPLITUDE)
            {
                colorscale.setRange((float) dataMin, (float) dataMax);
                colorscale.resetOpacity();
            }
            else
            {
                colorscale.setRange((float) minVolume, (float) maxVolume);
                colorscale.resetOpacity();
            }
        }
    }

	public StsColorscale getColorscale() { return colorscale; }
	public void setColorscale(StsColorscale colorscale) { this.colorscale = colorscale; }

    public void pick(GL gl, StsGLPanel glPanel)
    {
        float[] xyz;
        if (points == null) return;

        for (int i = 0; i < points.length; i++)
        {
            try
            {
                if(points[i] == null) continue;
                xyz = points[i].getPointValues();

                // Which volumes are being isVisible
                if ( ( (int) xyz[3] < startVol) || ( (int) xyz[3] > endVol))
                    continue;

                gl.glPushName(i);
                StsGLDraw.drawPoint(xyz, gl, computeSize(xyz) + 2);
                gl.glPopName();
            }
            catch(Exception e)
            {
                StsException.outputException("StsPointList.pick() failed.", e, StsException.WARNING);
            }
        }
    }

    public void mouseSelectedEdit(StsMouse mouse)
    {
      StsPickItem items = StsJOGLPick.pickItems[0];
      int pointIdx = items.names[1];
      logMessage(pointIdx);
    }

    public void showPopupMenu(StsGLPanel glPanel, StsMouse mouse) { }

    public void accumulateHistogram(int bindex)
    {
        if (bindex > 254)
        {
            bindex = 254;
        }
        if (bindex < 0)
        {
            bindex = 0;
        }
        dataCnt[bindex] = dataCnt[bindex] + 1;
        ttlHistogramSamples++;
    }

    private void accumulateHistogram(float value)
    {
        float scaledFloat = 254.0f * (float)(((double)value - dataMin) / (dataMax - dataMin));
        int scaledInt = StsMath.minMax(Math.round(scaledFloat), 0, 254);
        accumulateHistogram(scaledInt);
    }

    public void calculateHistogram()
    {
        for (int i = 0; i < 255; i++)
        {
            dataHist[i] = (float) ( (float)dataCnt[i] / (float)ttlHistogramSamples) * 100.0f;
        }
    }

    public void clearHistogram()
    {
        for (int i = 0; i < 255; i++)
        {
            dataCnt[i] = 0;
            dataHist[i] = 0.0f;
        }
        ttlHistogramSamples = 0;
    }

    public void setDataHistogram()
    {
        colorscaleBean.setHistogram(dataHist);
    }

    private int computeSize(float[] xyz)
    {
        double interval = (dataMax-dataMin)/(double)numberBins;
        int size = 2;
        if (sizeBy == VOLUME)
            size = (int) xyz[3] + 2;
        else if (sizeBy == AMPLITUDE)
        {
            int bin = (int) ((xyz[4] - dataMin) / interval);
//            int bin = (int) (Math.abs(Math.log(xyz[4]) - Math.log(dataMin)) / interval);
            size = bin + 1; // Set based on amplitude
        }
        else
            size = 3;
        return size;
    }

    public double getDataMin() { return dataMin; }
    public double getDataMax() { return dataMax; }
    public void setDataMin(double min) { dataMin = min; }
    public void setDataMax(double max) { dataMax = max; }

    public void logMessage(int pointIdx)
    {
        float[] point = points[pointIdx].getPointValues();
        double[] xy = model.getProject().getAbsoluteXYCoordinates(point);
        StsMessageFiles.infoMessage("X=" + xy[0] + " Y=" + xy[1] + " Z=" + point[2] + " PostStack3d=" + (int)point[3] + " Magnitude=" + point[4]);

        StsGLDraw.drawPoint(point, StsColor.BLACK, model.win3d.getGlPanel3d(), computeSize(point)+4);
        StsGLDraw.drawPoint(point, colorscale.getStsColor(colorIdx), model.win3d.getGlPanel3d(), computeSize(point), 2.0);
    }

    public boolean addPoint(StsPoint point)
    {
        try
        {
            if (nPoints >= points.length)
            {
                int newLength = points.length + initialSize;
                points = (StsPoint[]) StsMath.arraycopy(points, newLength, StsPoint.class);
            }
            points[nPoints++] = point;
            if(point.getPointValues()[3] > maxVolume)
                maxVolume = (int) point.getPointValues()[3];
            if(point.getPointValues()[3] < minVolume)
                minVolume = (int) point.getPointValues()[3];
            if(point.getPointValues()[4] > dataMax)
                dataMax = (int) point.getPointValues()[4];
            if(point.getPointValues()[4] < dataMin)
                dataMin = (int) point.getPointValues()[4];

            return true;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPointList.addPoint() failed.", e, StsException.WARNING);
            return false;
        }
    }

    private void setPoints(StsPoint[] points)
    {
        this.points = points;
        nPoints = points.length;
        for(int i=0; i<nPoints; i++)
        {
            float[] point = points[i].getPointValues();
            if(minVolume > point[3])
                minVolume = (int)point[3];
            if(maxVolume < point[3])
                maxVolume = (int)point[3];
            if(dataMin > point[4])
                dataMin = point[4];
            if(dataMax < point[4])
                dataMax = point[4];
        }
        if(volumes == null)
        {
            volumes = new String[maxVolume-minVolume+1];
            for(int i=0; i<maxVolume-minVolume+1; i++)
                volumes[i] = new String("Volume_" + i);
        }
        startVol = minVolume;
        endVol = maxVolume;
    }

    public void buildHistogram()
    {
        clearHistogram();
        for(int i=0; i<nPoints; i++)
        {
            float[] point = points[i].getPointValues();
            if(colorBy == VOLUME)
                accumulateHistogram(point[3]);
            else
                accumulateHistogram(point[4]);
        }
        calculateHistogram();
        resetColorscale();
        colorscaleBean.setHistogram(dataHist);
    }
    public void trimPointsArray()
    {
        points = (StsPoint[])StsMath.trimArray(points, nPoints);
        dbFieldChanged("maxVolume", maxVolume);
        dbFieldChanged("minVolume", minVolume);
        dbFieldChanged("nPoints", nPoints);
        dbFieldChanged("points", points);
    }

    public void display2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        if (glPanel3d == null) { return; }
        if (!isVisible) { return; }

        if(points.length < 1) { return; }

        GL gl = glPanel3d.getGL();
        if (gl == null) { return; }

        gl.glDisable(GL.GL_LIGHTING);
        glPanel3d.setViewShift(gl, StsGraphicParameters.edgeShift);
        displayLine2d(gl, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        glPanel3d.resetViewShift(gl);

        gl.glEnable(GL.GL_LIGHTING);
    }

    private void displayLine2d(GL gl, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed,
                               boolean yAxisReversed)
    {
        int bin = 0;
        StsColor color = new StsColor();
 //       double interval = (Math.log(dataMax) - Math.log(dataMin))/(double)numberBins;
        double interval = (dataMax - dataMin)/(double)(colorscale.getNColors()-1);

        if(isVisible == false)
            return;
        if(colorscale == null) initializeColorscale();
        try
        {
            for (int i = 0; i < nPoints; i++)
            {
                float[] xyz = points[i].getPointValues();
//                System.out.println("dirCoordinate=" + dirCoordinate + " XYZ[0]=" + xyz[0] + " XYZ[1]=" + xyz[1] + " XYZ[2]=" + xyz[2]);
                if(!showAllOn2d)
                {
                    switch (dirNo)
                    {
                        case 0:
                            if(dirCoordinate != xyz[0])
                                continue;
                            break;
                        case 1:
                            if(dirCoordinate != xyz[1])
                                continue;
                            break;
                        case 2:                    // Z
                            if(dirCoordinate != xyz[2])
                                continue;
                            break;
                    }
                }
                // If behind current slice make semi-transparent
                switch (dirNo)
                {
                    case 0:
                        if((!xAxisReversed) && (dirCoordinate > xyz[0]))
                            continue;
                        else if((xAxisReversed) && (dirCoordinate < xyz[0]))
                            continue;
                        break;
                    case 1:
                        if((!yAxisReversed) && (dirCoordinate < xyz[1]))
                            continue;
                        else if((yAxisReversed) && (dirCoordinate > xyz[1]))
                            continue;
                        break;
                    case 2:                    // Z
                        if(dirCoordinate < xyz[2])
                            continue;

                        if (axesFlipped)
                        {
                            float temp = xyz[0];
                            xyz[0] = xyz[1];
                            xyz[1] = temp;
                        }

                        break;
                }
//                System.out.println("i=" + i + " dirCoordinate=" + dirCoordinate + " XYZ[0]=" + xyz[0] + " XYZ[1]=" + xyz[1] + " XYZ[2]=" + xyz[2]);
                // Scale the amplitude
//                bin = (int)((Math.log(xyz[4]) - Math.log(dataMin))/interval);
                bin = (int)((xyz[4] - dataMin)/interval);

                // Which volumes are being isVisible
                if ( ( (int) xyz[3] < startVol) || ( (int) xyz[3] > endVol))
                    continue;

                // Colored by volume id or amplitude
                if(colorBy == VOLUME)
                    colorIdx = (int)xyz[3] - ((int)xyz[3]/colorscale.getNColors()) * colorscale.getNColors();
                else if(colorBy == AMPLITUDE)
                    colorIdx = bin;             // Set based on amplitude
                else
                    colorIdx = 0;

                int size = computeSize(xyz)+2;
                float[] pt2d = new float[2];
                switch (dirNo)
                {
                    case 1:                    // Y
                        pt2d[0] = xyz[0];
                        pt2d[1] = xyz[2];
                        break;
                    case 0:                    // X
                        pt2d[0] = xyz[1];
                        pt2d[1] = xyz[2];
                        break;
                    case 2:
                        pt2d[0] = xyz[0];
                        pt2d[1] = xyz[1];
                        break;
                }
                switch(colorBy)
                {
                    case VOLUME:
                        colorIdx = (int)xyz[3] - ((int)xyz[3]/colorscale.getNColors()) * colorscale.getNColors();
                        colorIdx = colorIdx % colorscale.getNColors();
                        color = colorscale.getStsColor(colorIdx);
                        break;
                    case AMPLITUDE:
                        colorIdx = bin;
                        colorIdx = colorIdx % colorscale.getNColors();
                        color = colorscale.getStsColor(colorIdx);
                        break;
                    case RED:
                        color = new StsColor(Color.RED);
                        break;
                    case GREEN:
                        color = new StsColor(Color.GREEN);
                        break;
                    case BLUE:
                        color = new StsColor(Color.BLUE);
                        break;
                    case YELLOW:
                      color = new StsColor(Color.YELLOW);
                      break;
                }
                StsGLDraw.drawPoint2d(pt2d, StsColor.BLACK, gl, computeSize(xyz) + 2);
                StsGLDraw.drawPoint2d(pt2d, color, gl, computeSize(xyz));
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPointList.display failed.", e, StsException.WARNING);
            return;
        }
    }

    public void display(StsGLPanel3d glPanel3d)
    {
        int bin = 0;
        StsColor color = new StsColor();
        float[] blackFloats = new float[] {0.0f, 0.0f, 0.0f, 1.0f};
//        StsColor black = new StsColor(blackFloats);
//        double interval = (Math.log(dataMax) - Math.log(dataMin))/(double)numberBins;

        if(isVisible == false)
            return;
        if(colorscale == null) initializeColorscale();
//        double interval = (Math.log(dataMax) - Math.log(dataMin))/(double)colorscale.getNColors();
        double interval = (dataMax- dataMin)/(double)(colorscale.getNColors()-1);

        try
        {
			float zscale = glPanel3d.getZScale();
            for (int i = 0; i < nPoints; i++)
            {
                float[] xyz = points[i].getPointValues();

                // Scale the amplitude
 //               bin = (int)((Math.log(xyz[4]) - Math.log(dataMin))/interval);
                bin = (int)((xyz[4] - dataMin)/interval);

                // Which volumes are being isVisible
                if ( ( (int) xyz[3] < startVol) || ( (int) xyz[3] > endVol))
                    continue;

                // Colored by volume id or amplitude
                switch(colorBy)
                {
                    case VOLUME:
                        colorIdx = (int)xyz[3] - ((int)xyz[3]/colorscale.getNColors()) * colorscale.getNColors();
                        colorIdx = colorIdx % colorscale.getNColors();
                        color = colorscale.getStsColor(colorIdx);
                        blackFloats[3] = color.alpha;
                        break;
                    case AMPLITUDE:
                        colorIdx = bin;
                        colorIdx = colorIdx % colorscale.getNColors();
                        color = colorscale.getStsColor(colorIdx);
                        blackFloats[3] = color.alpha;
                        break;
                    case RED:
                        color = new StsColor(Color.RED);
                        break;
                    case GREEN:
                        color = new StsColor(Color.GREEN);
                        break;
                    case BLUE:
                        color = new StsColor(Color.BLUE);
                        break;
                    case YELLOW:
                      color = new StsColor(Color.YELLOW);
                      break;
                }
                if(pointType == POINT)
                {
                    StsColor black = new StsColor(blackFloats);
                    StsGLDraw.drawPoint(xyz, black, glPanel3d, computeSize(xyz) + 2);
                    StsGLDraw.drawPoint(xyz, color, glPanel3d, computeSize(xyz), 2.0);
                }
                else if (pointType == SPHERE)
                    StsGLDraw.drawSphere(glPanel3d, xyz, color, computeSize(xyz) * 5);
                else if (pointType == DISK)
                    StsGLDraw.drawCylinder(glPanel3d, xyz, color, computeSize(xyz) * 5, 5);
            }
        }
        catch(Exception e)
        {
            StsException.outputException("StsPointList.display failed.", e, StsException.WARNING);
            return;
        }
    }
	public void actionPerformed(ActionEvent e)
	{
		if (e.getSource() instanceof StsColorscale)
		{
			currentModel.win3dDisplayAll();
		}
	}

    public void setVolumeNames(String[] names)
    {
        volumes = names;
    }
    public void addVolumeNames(String[] names)
    {
        volumes = (String[])StsMath.arrayAddArray(volumes, names);
		dbFieldChanged("volumes", volumes);
    }

    public int getNumberBins() { return numberBins; }
    public void setNumberBins(int bins)
    {
        numberBins = bins;
        dbFieldChanged("numberBins", numberBins);
        model.win3dDisplayAll();
    }

    public void setPointSetType(byte val)
    {
        setType = val;
    }

    public byte getPointSetType()
    {
        return setType;
    }
    public String getSetTypeAsString()
    {
        return SETTYPES[setType];
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields != null) return displayFields;
        symbolListBean = new StsComboBoxFieldBean(StsPointList.class, "symbolString", "Symbol:", SYMBOL_TYPE_STRINGS);
        sizeListBean = new StsComboBoxFieldBean(StsPointList.class, "sizeByString", "Size By:", SIZE_TYPE_STRINGS);
        colorListBean = new StsComboBoxFieldBean(StsPointList.class, "colorByString", "Color By:", COLOR_TYPE_STRINGS);
        numBinsBean = new StsIntFieldBean(StsPointList.class, "numberBins", 10, 200, "Amplitude Scale:");
        startVolumeBean = new StsIntFieldBean(StsPointList.class, "startVolume", 0, 10000, "Start PostStack3d:", true);
        endVolumeBean = new StsIntFieldBean(StsPointList.class, "endVolume", 0, 10000, "End PostStack3d:", true);
        volumeNameBean = new StsStringFieldBean(StsPointList.class, "currentVolumeName", true, "End PostStack3d Name:");
        colorscaleBean = new StsEditableColorscaleFieldBean(StsPointList.class, "colorscale");

        displayFields = new StsFieldBean[]
        {
            new StsStringFieldBean(StsPointList.class, "name", true, "Name"),
            new StsDoubleFieldBean(StsPointList.class, "dataMin", false, "Data Min"),
            new StsDoubleFieldBean(StsPointList.class, "dataMax", false, "Data Max"),
            new StsIntFieldBean(StsPointList.class, "numPoints", false, "Number of Points"),
            new StsStringFieldBean(StsPointList.class,"setTypeAsString", false, "Set Type:"),
            new StsBooleanFieldBean(StsPointList.class, "isVisible", "Enable:"),
            startVolumeBean,
            endVolumeBean,
            volumeNameBean,
            symbolListBean,
            colorListBean,
            sizeListBean,
            numBinsBean,
            colorscaleBean
        };
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        return null;
    }

    public Object[] getChildren()
    {
        return new Object[0];
    }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null)
        {
            objectPanel = StsObjectPanel.constructor(this, true);
        }
        return objectPanel;
    }

    public boolean anyDependencies()
    {
        return false;
    }

    public void treeObjectSelected()
    {
//        startVolumeBean.setMaximum(maxVolume);
//        startVolumeBean.setMinimum(minVolume);
//        endVolumeBean.setMaximum(maxVolume);
//        endVolumeBean.setMinimum(minVolume);
    }

    public void setIsVisible(boolean b)
    {
        if (b == isVisible)
        {
            return;
        }
        isVisible = b;
        super.setIsVisible(isVisible);
        currentModel.win3dDisplayAll();
    }
    public boolean getIsVisible() { return isVisible; }

    public void setStartVolume(int vol)
    {
        if(vol > volumes.length)
            vol = volumes.length -1;
        startVol = vol;
        model.win3dDisplayAll();
    }
    public int getStartVolume() { return startVol; }
    public void setEndVolume(int vol)
    {
        if(vol > volumes.length)
            vol = volumes.length -1;
        endVol = vol;
        if(volumes != null) setCurrentVolumeName(volumes[endVol]);
        model.win3dDisplayAll();
    }
    public int getEndVolume() { return endVol; }
    public void setCurrentVolume(int vol)
    {
        startVol = vol;
        endVol = vol;
        currentVol = vol;
        if(volumes != null) setCurrentVolumeName(volumes[endVol]);
        model.win3dDisplayAll();
    }
    public int getCurrentVol() { return currentVol; }

    public String getCurrentVolumeName() { return currentVolumeName; }
    public void setCurrentVolumeName(String name)
    {
        currentVolumeName = name;
        StsObjectPanel panel = model.getProject().getObjectPanel();
        if(panel != null)
        {
            if(model.getCurrentObject(StsPointList.class) == this)
                volumeNameBean.setValue(currentVolumeName);
        }
    }

    public boolean getIsAnimated() { return animated; }
    public void setIsAnimated(boolean bool)
    {
        animated = bool;
    }
    public int getMinVolume() { return minVolume; }
    public int getMaxVolume() { return maxVolume; }
    public int getNumberVolumes() { return volumes.length; }
    public int getNumPoints() { return nPoints; }
    public boolean getShowAllOn2d() { return showAllOn2d; }
    public void setShowAllOn2d(boolean showAll)
    {
        showAllOn2d = showAll;
        model.win3dDisplayAll();
    }

    public void incrementCurrentVolume()
    {
        currentVol++;
        if(currentVol > maxVolume)
            currentVol = minVolume;
    }
    public void setNextCummlativeVolume(int vol)
    {
        if(vol == minVolume)
        {
            setStartVolume(vol);
            setEndVolume(vol);
        }
        else if(vol > maxVolume)
        {
            setStartVolume(minVolume);
            setEndVolume(maxVolume);
        }
        else
            setEndVolume(vol);

        updatePanel();
        model.win3dDisplayAll();
    }
    public void updatePanel()
    {
        startVolumeBean.getValueFromPanelObject();
        endVolumeBean.getValueFromPanelObject();
        endVolumeBean.getValueFromPanelObject();
    }
    public String getSizeByString()
    {
        return SIZE_TYPE_STRINGS[sizeBy];
    }
    public void setSizeByString(String sizeString)
    {
        for(int i=0; i<SIZE_TYPE_STRINGS.length; i++)
        {
            if (sizeString.equals(SIZE_TYPE_STRINGS[i]))
            {
                sizeBy = SIZE_TYPES[i];
                break;
            }
        }
        dbFieldChanged("sizeBy", sizeBy);
        model.win3dDisplayAll();
        return;
    }

    public String getSymbolString()
    {
        return SYMBOL_TYPE_STRINGS[pointType];
    }

    public void setSymbolString(String symbolString)
    {
        for(int i=0; i<SYMBOL_TYPE_STRINGS.length; i++)
        {
            if (symbolString.equals(SYMBOL_TYPE_STRINGS[i]))
            {
                pointType = SYMBOL_TYPES[i];
                break;
            }
        }
        dbFieldChanged("pointType", pointType);
        model.win3dDisplayAll();
        return;
    }

    public String getColorByString()
    {
        return COLOR_TYPE_STRINGS[colorBy-1];
    }

    public void setColorByString(String color)
    {
        for(int i=0; i<COLOR_TYPE_STRINGS.length; i++)
        {
            if (color.equals(COLOR_TYPE_STRINGS[i]))
            {
                colorBy = COLOR_TYPES[i];
                buildHistogram();
                break;
            }
        }
        model.win3dDisplayAll();
        dbFieldChanged("colorBy", colorBy);
        return;
    }

    public boolean canExport() { return true; }
    public boolean export()
	{
        progressBarDialog = StsProgressBarDialog.constructor(currentModel.win3d, "Pointset Export", false);
        progressBarDialog.setLabelText("Exporting point set with name" + getName());
        Runnable runExport = new Runnable()
        {
            public void run()
            {
                exportPointset();
            }
        };

        Thread exportThread = new Thread(runExport);

        exportThread.start();
        return true;
	}

    public void exportPointset()
    {
        int maxProgress = nPoints - 1;
        PrintWriter printWriter = null;
        progressBarDialog.setProgressMax(maxProgress);
        try
        {
            StsProject project = currentModel.getProject();
            String filename = project.getRootDirString() + getName() + ".csv";
            File file = new File(filename);
            if (file.exists())
            {
                boolean overWrite = StsYesNoDialog.questionValue(currentModel.win3d,
                    "File " + filename + " already exists. Do you wish to overwrite it?");
                if (!overWrite)
                {
                    return;
                }
            }
            progressBarDialog.setLabelText("Exporting Pointset to File:    " + filename);
            progressBarDialog.pack();
            progressBarDialog.setVisible(true);

            // Output file.
            printWriter = new PrintWriter(new FileWriter(filename, false));
            int pointSize = points[0].getPointValues().length;
            switch(pointSize)
            {
                case 3:
                    printWriter.println(" X,Y,Z ");
                    break;
                case 4:
                    printWriter.println(" X,Y,Z,PostStack3d");
                    break;
                case 5:
                    printWriter.println(" X,Y,Z,PostStack3d,Amplitude");
                    break;
                default:
                    break;
            }
            for (int i = 0; i < nPoints; i++)
            {
                float[] point = points[i].getPointValues();
                switch(pointSize)
                {
                    case 3:
                        printWriter.println(point[0] + "," + point[1] + "," + point[2]);
                        break;
                    case 4:
                        printWriter.println(point[0] + "," + point[1] + "," + point[2] + ","
                                            + point[3]);
                        break;
                    case 5:
                        printWriter.println(point[0] + "," + point[1] + "," + point[2] + ","
                                            + point[3] + "," + point[4]);
                        break;
                    default:
                        break;
                }
                progressBarDialog.setProgress(i);
            }
        }
        catch (Exception e)
        {
            StsException.outputException("StsPointSet.export() failed.", e, StsException.WARNING);
            return;
        }
        finally
        {
            if (printWriter != null)
            {
                printWriter.flush();
                printWriter.close();
            }
            progressBarDialog.setProgress(maxProgress);
//            progressBarDialog.dispose();
        }
    }

	public void writeObject(StsDBOutputStream out) throws IllegalAccessException, IOException
	{
		for (int i = 0; i < nPoints; i++)
		{
			float[] v = points[i].v;
			out.writeInt(v.length);
			for (int j = 0; j < v.length; j++)
			{
				out.writeFloat(v[j]);
			}
		}
	}

	public void readObject(StsDBInputStream in) throws IllegalAccessException, IOException
	{
		if (nPoints > 0)
		{
			points = new StsPoint[nPoints];
			for (int i = 0; i < nPoints; i++)
			{
				int len = in.readInt();
				points[i] = new StsPoint(len);
				float[] v = points[i].v;
				for (int j = 0; j < len; j++)
				{
					v[j] = in.readFloat();
				}
			}
		}

		updatePanel();
	}

    public void exportObject(StsDBFileObjectTrader objectTrader) { } 
}
