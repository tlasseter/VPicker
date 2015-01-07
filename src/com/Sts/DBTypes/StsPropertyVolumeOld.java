
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

public class StsPropertyVolumeOld extends StsObject
{
    // instance fields
    // propertyVolume is dimensioned [nRows-1][nCols-1][nLayers]
    protected String propertyName;
//    protected File hdfPropertyFile = null;
    protected String scatteredFileName;
	protected int nRows, nCols;
    protected int nLayers;
    protected float angle = 0.0f;
    protected float xInc, yInc;
    protected float xMin, yMin;
    protected StsSpectrum propertySpectrum = null;
    protected float colorScale = 0.0f;
    protected boolean gridAligned = true;                 /** true if grid same as S2S grid */
//    protected float gridRowStart, gridColStart; // If not grid aligned, this maps global grid to prop grid
    protected float gridRowInc, gridColInc;

    protected boolean inMemory = true;
    static int maxValuesInMemory = 1000000;

    // property values
    transient private float[] propertyValues = null;
    transient private float[] valueRange = null;
//    transient private HDF3dModelFile hdfPropertyFileObject = null;

    transient StsWin3d win3d;  // for mouse & message handling
    transient StsModel model;

    static boolean debug = false;
    static public final float nullValue = StsParameters.nullValue;

    static public final StsFieldBean[] fieldBeans =
    {
        new StsBooleanFieldBean(StsPropertyVolumeOld.class, "isVisible", "Enable")
    };

    static public final StsFieldBean[] propertyFields =
    {
        new StsFloatFieldBean(StsPropertyVolumeOld.class, "colorScale"),
        new StsFloatFieldBean(StsPropertyVolumeOld.class, "rangeMax", "Max"),
        new StsFloatFieldBean(StsPropertyVolumeOld.class, "rangeMin", "Min")
    };

    /** constructor for property volume with default x/y sizes & increments */
    public StsPropertyVolumeOld(StsModel model, String propertyName,
                             StsGridDefinition gridDefinition, int nLayers)
    {
        // save inputs
        this.model = model;
        this.win3d = model.win3d;
        this.xInc = gridDefinition.getXInc();
        this.nCols = gridDefinition.getNCols();
        this.yInc = gridDefinition.getYInc();
        this.nRows = gridDefinition.getNRows();
        this.nLayers = nLayers;

        setPropertyName(propertyName);

        int nValues = (nRows-1)*(nCols-1)*nLayers;
        inMemory = nValues < maxValuesInMemory;
    }

    /** constructor for DB */
	public StsPropertyVolumeOld()
	{
	}

    public StsPropertyVolumeOld(StsSurface surface, StsModel model, String propertyName)
    {
        // save inputs
        this.model = model;
        this.win3d = model.win3d;
        this.xInc = surface.getXInc();
        this.nCols = surface.getNCols();
        this.yInc = surface.getYInc();
        this.nRows = surface.getNRows();
        this.nLayers = 1;

        setPropertyName(propertyName);

        int nValues = (nRows-1)*(nCols-1)*nLayers;
        inMemory = nValues < maxValuesInMemory;
    }

    /** set values from database */
    public boolean initialize(StsModel model)
    {
        try
        {
            this.model = model;
            this.win3d = model.win3d;
//            if (model.getProperty() == this) readPropertyVolumeFromHDF();
        }
        catch(Exception e)
        {
            StsMessageFiles.logMessage("Unable to classInitialize property volume.");
            return false;
        }

        return true;
    }

   /** accessors */
    public int getNX() { return nRows; }
    public int getNY() { return nCols; }
    public int getNLayers() { return nLayers; }
    public float getXInc() { return xInc; }
    public float getYInc() { return yInc; }
    public void setPropertyName(String name)
    {
        propertyName = name;
        scatteredFileName = getScatteredFileName(model.getName(), propertyName);
    }
    static public String getScatteredFileName(String modelName, String propertyName)
    {
        return "DB." + modelName + ".P." + propertyName + ".scat";
    }
    public String getPropertyName() { return propertyName; }
/*
    public void setHDFPropertyFile(File hdfPropertyFile, HDF3dModelFile hdfPropertyFileObject)
    {
        this.hdfPropertyFile = hdfPropertyFile;
        this.hdfPropertyFileObject = hdfPropertyFileObject;
    }
    public File getHDFPropertyFile() { return hdfPropertyFile; }
*/
    public String getScatteredFileName() { return scatteredFileName; }


    public void setIsVisible(boolean b)
    {
//		if( currentModel != null )
//        	currentModel.setProperty(b ? this : null);
//        if(b && getValues() == null) readPropertyVolumeFromHDF();
    }

    public boolean getIsVisible()
    {
//		if( currentModel != null ) return currentModel.getProperty() == this;
        return false;
    }

    public void setName(String name) { setPropertyName(name); }
    public String getName() { return getPropertyName(); }

	public float getColorScale() { return this.colorScale; }
	public void setColorScale(float colorScale) { this.colorScale = colorScale; }
/*
    public StsColorscale getColorscale()
    {
    	return new StsColorscale(getPropertySpectrum(), valueRange);
    }
*/
    public void setRangeMin(float min)
    {
    	valueRange[0] = min;
        setSpectrumScale();
//		if( currentModel != null ) currentModel.setProperty(this);
    }
    public void setRangeMax(float max)
    {
    	valueRange[1] = max;
        setSpectrumScale();
//		if( currentModel != null ) currentModel.setProperty(this);
    }
    public float getRangeMin() { return valueRange[0]; }
    public float getRangeMax() { return valueRange[1]; }

    /** set/get property spectrum to use in property value display */
    public void setPropertySpectrum(String name)
    {
//        setPropertySpectrum(model.getSpectrum(name));
    }
    public void setPropertySpectrum(StsSpectrum spectrum)
    {
        propertySpectrum = spectrum;
//        setSpectrumScale();
    }
    public StsSpectrum getPropertySpectrum() { return propertySpectrum; }

    /** get property values from the volume */
    public float[] getValues() { return propertyValues; }  /** all the values */

   /** set values */
/*
    public void setValues(float[] values)
    {
        propertyValues = values;
        calcMinMax();
    }
*/
    public float[] getValueRange() { return valueRange; }

    public void setGridAligned(boolean aligned)
    {
        gridAligned = aligned;
        if(!aligned) computeRowColMapping();
    }

    private void computeRowColMapping()
    {
        StsGridDefinition gridDef = currentModel.getGridDefinition();
        int gridNRows = gridDef.getNRows();
        int gridNCols = gridDef.getNCols();
        float gridXInc = gridDef.getXInc();
        float gridYInc = gridDef.getYInc();
//        gridRowStart = (gridYMin - yMin)/yInc;
//        gridColStart = (gridXMin - xMin)/xInc;
        gridRowInc = gridYInc/yInc;
        gridColInc = gridXInc/xInc;
    }

    public boolean getGridAligned() { return gridAligned; }

    /** read the property values into a 1-D array
        Z is fastest dimension, then X, then Y */
    public boolean readPropertyVolumeFromHDF()
    {
        try
        {
/*
            if (hdfPropertyFileObject == null)
            {
          	    // open the HDF file & build the model object
                SDSFile hdfFile = new SDSFile(hdfPropertyFile, SDSFile.ACCESS_READ);
                hdfPropertyFileObject = new HDF3dModelFile(hdfFile, propertyName);
            }

            if(!inMemory)
                propertyValues = hdfPropertyFileObject.getSampleValues();
            else
                propertyValues = hdfPropertyFileObject.getValues();

            if(propertyValues == null) return false;

            valueRange = computeValueRange(propertyValues);
            if(!inMemory) propertyValues = null;

            if(valueRange == null) return false;
            setSpectrumScale();
*/
            return true;
        }
        catch (Exception e) { return false; }
    }

    public float[][] getValuesArray(float [][]rowColFs, int nFirstLayer, int nLayers)
    {
        int nCols = rowColFs.length;
        float[][] values = new float[nCols][];
//        for(int n = 0; n < nCols; n++)
//            values[n] = hdfPropertyFileObject.getColValues((int)rowColFs[n][0], (int)rowColFs[n][1], nFirstLayer, nLayers);

        return values;
    }

    public float[][] getValuesArray(int [][]rowCols, int nFirstLayer, int nLayers)
    {
        int nCols = rowCols.length;
        float[][] values = new float[nCols][];
//        for(int n = 0; n < nCols; n++)
//            values[n] = hdfPropertyFileObject.getColValues(rowCols[n][0], rowCols[n][1], nFirstLayer, nLayers);

        return values;
    }

    public StsColor[][] getSliceColors(StsGridBoundingBox boundingBox, int nSlice)
    {
        try
        {
            int rowMin = boundingBox.rowMin;
            int rowMax = boundingBox.rowMax;
            int colMin = boundingBox.colMin;
            int colMax = boundingBox.colMax;
            float[] values = null;
//            float[] values = hdfPropertyFileObject.getSliceValues(rowMin, rowMax, colMin, colMax, nSlice, nLayers);
            int nRows = rowMax-rowMin;
            int nCols = colMax-colMin;
            StsColor[][] colors = new StsColor[nRows][nCols];
            int n = 0;
            for(int row = 0; row < nRows; row++)
                for(int col = 0; col < nCols; col++, n++)
                    colors[row][col] = getScaledColor(values[n]);

            return colors;
        }
        catch(Exception e)
        {
            StsException.outputException("StsPropertyVolumeOld.getSliceColors() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    /** set property values to null for garbage collection */
    public void doneWithPropertyValues() { propertyValues = null; }

    /* calculate min and max for all the values (assumes no nulls) */

    private float[] computeValueRange(float[] values)
    {
        float value;
        int i;

        if (values==null) return null;
        float[] valueRange = new float[2];

        int nValues = propertyValues.length;

        // find first non-null and set starting range values
        for (i=0; i<nValues; i++)
        {
            value = propertyValues[i];
        /*
            if(value != StsParameters.HDFnullValue)
            {
                valueRange[0] = value;
                valueRange[1] = value;
                break;
            }
        */
        }
        // now compare the remaining non-null values to expand range
        for (i++; i<nValues; i++)
        {
            value = propertyValues[i];
//            if(value == StsParameters.HDFnullValue) continue;
            if(value < valueRange[0]) valueRange[0] = value;
            else if(value > valueRange[1]) valueRange[1] = value;
        }
        return valueRange;
    }

    /** create scale for property spectrum-value correspondance */
    private void setSpectrumScale()
    {
        if (propertySpectrum==null || valueRange==null) return;
        colorScale = (valueRange[1] - valueRange[0]) / (float)propertySpectrum.getNColors();
    }

    /** check ability to return a scaled color for a property value */
    public boolean scaledColorsAreValid()
    {
        if (colorScale==0.0f || propertySpectrum==null) return false;
        StsColor[] colors = propertySpectrum.getStsColors();
        if (colors==null) return false;
        return true;
    }

    public StsColor getScaledColor(int row, int col, int nLayer)
    {
        float value = getValue(row, col, nLayer);
        return getScaledColor(value);
    }

    public float getValue(int row, int col, int nLayer)
    {
        try
        {
            if(!gridAligned)
            {
//                row = (int)((row - gridRowStart)*gridRowInc);
//                col = (int)((col - gridColStart)*gridColInc);
            }

            if(!insideGrid(row, col)) return nullValue;

            // nLayer may be == nLayers (when we are drawing the very bottom slice);
            // in this case, limit it to nLayers-1.
            nLayer = StsMath.minMax(nLayer, 0, nLayers-1);

            if(propertyValues != null)
                return propertyValues[(row*(nCols-1) + col)*nLayers + nLayer];
            else
            {
                float[] value = null;
//                float[] value = hdfPropertyFileObject.getValue(row, col, nLayer);
//                if(value == null) return nullValue;
                return value[0];
            }
        }
        catch (Exception e) { return nullValue; }
    }

    public StsColor getScaledColor(float value)
    {
        if (value == StsParameters.nullValue) return StsColor.BLACK;

        StsColor[] colors = propertySpectrum.getStsColors();
        int c = (int)((value - valueRange[0])/colorScale);
        if (c < 0 || c >= colors.length) return StsColor.BLACK;
        return colors[c];
    }

    public float getValue(float xValue, float yValue, int nLayer)
    {
        int col = (int)((xValue-xMin)/xInc);
        int row = (int)((yValue-yMin)/yInc);
        return getValue(row, col, nLayer);
    }

    public boolean insideGrid(int row, int col)
    {
        return row >= 0 && row < nRows && col >= 0 && col < nCols;
    }

    /** test program */
    /*
    public static void main(String[] args)
    {
        try
        {
            // build model with read project file
            File f = new File("../../../../texas_data/proj.texas");
            String path = f.getAbsolutePath();
		    StsModel model = new StsModel(new StsProject(path));
            String up = ".." + File.separator;
            String up4 = up+up+up+up;
            model.getProject().setDataDirectory(new File(up4 + "texas_data"));

            // read surface files
            float nullZValue = model.project.getMapGenericNull();
            f = new File(up4 + "texas_data" + File.separator + "map.generic.TOP_10000sand");
            StsImportAsciiSurfaces.createSurfaceFromAscii(model.win3d, f.getParent(), f.getName(), nullZValue);
            StsSurface s1 = model.getSurface("TOP_10000sand");
            //f = new File("..\\..\\..\\..\\texas_data\\map.generic.BASE_10000");
            //StsImportAsciiSurfaces.createSurfaceFromAscii(model.win3d, f.getParent(), f.getName(), nullZValue);
            //StsSurface s2 = model.getSurface("BASE_10000");

            // build property volume
            StsSurface s = s1;
            StsGrid grid = s.getGrid();
            StsGridDefinition gridDefinition = new StsGridDefinition(grid.getXMin(), grid.getYMin(),
                grid.getXMax(), grid.getYMax(), grid.getXInc(), grid.getYInc());

            //String pn = "TestProperty";
            String pn = "ILD";
            StsPropertyVolumeOld pv = new StsPropertyVolumeOld(model, pn, gridDefinition, 5);

            String hdfFilename = "vol.hdf";
            //if (!pv.readPropertyVolumeFromHDF(hdfFilename))
            if (!pv.readPropertyVolumeFromHDF())
            {
                System.out.println("Unable to read: " + pv.volumeFileName);
                return;
            }

            // set property values in surfaces
            try { pv.setPropertySpectrum("Basic"); }
            catch (StsException e) { }
            s1.setProperty(pv, 1);
//          pv.setValuesOnSurface(s1, 1);
            //s2.setPropertySpectrum(model.getSpectrum("Basic"));
            //s2.setPropertyName(pn + "_1");
            //s2.setPropertyValues(pv.getValues(s2, 1));

            // build model with read project file
            StsModel model = new StsModel();
            File f = new File("..\\..\\..\\..\\texas_data\\proj.texas");
            String path = f.getAbsolutePath();
            model.readProjectParameters(path);
            model.setName("test");

            // create a property volume object
            float xMin = 1872974.5f;
            float xMax = 1946949.5f;
            int nX = 100;
            float xInc = (xMax-xMin)/(float)nX;
            float yMin = 252215.81f;
            float yMax = 306390.8f;
            int nY = 100;
            float yInc = (yMax-yMin)/(float)nY;
            int nSubZones = 5;
            StsPropertyVolumeOld vol = new StsPropertyVolumeOld(model, "ILD", xMin,
                    xInc, nX, yMin, yInc, nY, nSubZones);

            // read property values from HDF
            f = new File ("..\\..\\..\\..\\texas_data\\vol.hdf.test.NS.5.P.ILD");
            path = f.getAbsolutePath();
            System.out.println("HDFSDSFile:  " + path);
            vol.readPropertyVolumeFromHDF();
            float[] values = vol.getValues();
            if (values==null) { System.out.println("Values = null"); return; }
            else System.out.println("Values length = " + values.length);
            float value = vol.getValue(10);
            System.out.println("value(10) = " + value);
            value = vol.getValue(1, 1, 1);
            System.out.println("value(1,1,1) = " + value);
            value = vol.getValue(1900000.0f, 275000.0f, 3);
            System.out.println("value(1900000.0f, 275000.0f, 3) = " + value);

        }
        catch (Exception e)
        {
            StsException.outputException("StsPropertyVolumeOld.main: ", e, StsException.FATAL);
        }
    }
    */
}
