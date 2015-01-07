
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;
import com.Sts.IO.*;

import java.util.*;

public class StsPropertyVolume extends StsGridBoundingBox
{
    transient StsPropertyType propertyType;
    transient protected StsGridBoundingBox boundingBox;
    transient private byte distributionType = PROP_NONE;
    transient protected String filename;
    transient protected boolean inMemory = false;

    // property values
    transient protected float[] values = null;
    transient protected float valueMin = StsParameters.largeFloat;
    transient protected float valueMax = -StsParameters.largeFloat;

    static int maxValuesInMemory = 1000000;

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

    static public final byte PROP_NONE = 0;
    static public final byte PROP_CONSTANT = 1;
    static public final byte PROP_VARIABLE_Z = 2;
    static public final byte PROP_VARIABLE_XY = 3;
    static public final byte PROP_VARIABLE_XYZ = 4;

    /** constructor for DB */
	public StsPropertyVolume()
	{
	}

    /** constructor for temporary property volume files (not saved to disk). */
    public StsPropertyVolume(StsPropertyTypeDefinition propertyTypeDefinition, StsGridBoundingBox zone)
    {
        super(false);
        this.propertyType = new StsPropertyType(propertyTypeDefinition);
        setName(propertyType.name);
        initializeCellGridVolume(zone);
        this.boundingBox = zone;
    }

    public StsPropertyVolume(StsPropertyType propertyType, StsBlock block)
    {
         super(false);
         initializeCellGridVolume(block);
         setName(propertyType.name);
         this.propertyType = propertyType;
         this.boundingBox = block;
    }

    /** constructor for property volume files which have been saved to disk. Parameters will be read from file. */
     public StsPropertyVolume(StsPropertyType propertyType, StsBlock block, String outputPathname)
     {
         this(propertyType, block);
         setFilename(propertyType, block, outputPathname);
     }

     /** New property volume */
     public StsPropertyVolume(StsPropertyType propertyType, StsBlock block, byte distributionType, String outputPathname)
     {
         this(propertyType, block);
         setDistributionType(distributionType);
         setFilename(propertyType, block, outputPathname);
     }
    
    /** constructor when reloading a propertyVolume file; parameters will be gotten from file. */
    public StsPropertyVolume(StsPropertyType propertyType)
    {
        this.propertyType = propertyType;
    }

    protected void setFilename(StsPropertyType propertyType, StsBlock block, String outputPathname)
    {
        filename = outputPathname + propertyType.eclipseName + ".block." + block.getIndex();
    }

    /** set values from database */
    public boolean initialize(StsModel model)
    {
        return true;
    }

    public String getName() { return propertyType.name; }

    public int getNPoints()
    {
        if(distributionType == PROP_CONSTANT)
            return 1;
        else if(distributionType == PROP_VARIABLE_Z)
            return getNSlices();
        else if(distributionType == PROP_VARIABLE_XY)
            return getNRows()*getNCols();
        else if(distributionType == PROP_VARIABLE_XYZ)
            return  getNRows()*getNCols()*getNSlices();
        else
        {
            StsException.systemError(this, "getNPoints", "Undefined property distributionType " + distributionType + " for property " + toDetailString());
            return 1;
        }
    }

    public void setDistributionType(byte distributionType)
    {
        this.distributionType = distributionType;
        int nValues = getNPoints();
        values = new float[nValues];
    }

    public void setConstant(float value)
    {
        distributionType = PROP_CONSTANT;
        values = new float[] { value };
    }

    public void setVariableZ(float[] values)
    {
        distributionType = PROP_VARIABLE_Z;
        int nValues = getNPoints();
        inMemory = nValues < maxValuesInMemory;
        if(nValues != values.length)
            StsException.systemError(this, "setVariableZ", "Number of input Values (" + values.length + ") not equal to required " + nValues);
        setValues(values);
    }

    public void setVariableZ(float[] values, int offset, int length)
    {
        distributionType = PROP_VARIABLE_Z;
        int nValues = getNPoints();
        inMemory = nValues < maxValuesInMemory;
        if(this.values == null) this.values = new float[values.length];
        setValues(values, offset, length);
    }

    public void setVariableXY(float[] values)
    {
        distributionType = PROP_VARIABLE_XY;
        int nValues = getNPoints();
        inMemory = nValues < maxValuesInMemory;
        if(nValues != values.length)
            StsException.systemError(this, "setVariableXY", "Number of input Values (" + values.length + ") not equal to required " + nValues);
        setValues(values);
    }

    public void setVariableXYZ(float[] values)
    {
        setVariableXYZ();
        setValues(values);
    }

    public void setVariableXYZ()
    {
        distributionType = PROP_VARIABLE_XYZ;
        int nValues = getNPoints();
        inMemory = nValues < maxValuesInMemory;
        if(nValues != values.length)
            StsException.systemError(this, "setVariableXYZ", "Number of input Values (" + values.length + ") not equal to required " + nValues);
    }

    private void setValues(float[] values)
    {
        this.values = values;
        setRange();
    }

    private void setValues(float[] values, int offset, int length)
    {
        System.arraycopy(values, offset, this.values, 0, length);
        setRange();
    }

    public void initializeValues(float value)
    {
        Arrays.fill(values, value);
    }

    static public byte getDistributionType(byte currentType, byte newType)
    {
        byte type;
        if(newType > currentType)
        {
            if(newType == PROP_VARIABLE_XY && currentType == PROP_VARIABLE_Z)
                type = PROP_VARIABLE_XYZ;
            else
                type = newType;
        }
        else if(newType < currentType)
        {
            if(currentType == PROP_VARIABLE_XY && newType == PROP_VARIABLE_Z)
                type = PROP_VARIABLE_XYZ;
            else
                type = currentType;
        }
        else
            type = currentType;

        return type;
    }

   /** accessors */
    public int getNX() { return getNRows(); }
    public int getNY() { return getNCols(); }
    public int getNLayers() { return getNSlices(); }
    public float getXInc() { return xInc; }
    public float getYInc() { return yInc; }

    static public String getScatteredFileName(String modelName, String propertyName)
    {
        return "DB." + modelName + ".P." + propertyName + ".scat";
    }

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

    public void setRangeMin(float min)
    {
    	valueMin = min;
    }
    public void setRangeMax(float max)
    {
    	valueMax = max;
    }
    public float getRangeMin() { return valueMin; }
    public float getRangeMax() { return valueMax; }

    public void setRange(float min, float max)
    {
        valueMin = min;
        valueMax = max;
    }

    /** get property values from the volume */
    public float[] getValues()
    {
        if(!checkLoadData()) return null;
        return values;
    }

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
        if(!checkLoadData()) return null;
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
    public void doneWithPropertyValues() { values = null; }

    /* calculate min and max for all the values (assumes no nulls) */

    private float[] computeValueRange(float[] values)
    {
        float value;
        int i;

        if (values==null) return null;
        float[] valueRange = new float[2];

        int nValues = this.values.length;

        // find first non-null and set starting range values
        for (i=0; i<nValues; i++)
        {
            value = this.values[i];
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
            value = this.values[i];
//            if(value == StsParameters.HDFnullValue) continue;
            if(value < valueRange[0]) valueRange[0] = value;
            else if(value > valueRange[1]) valueRange[1] = value;
        }
        return valueRange;
    }

    public StsColor getScaledColor(int row, int col, int nLayer)
    {
        float value = getValue(row, col, nLayer);
        return getScaledColor(value);
    }

    public float getValue(int row, int col, int layer)
    {
        if(!checkLoadData())
            return nullValue;
        try
        {
            if(distributionType == PROP_CONSTANT)
                return values[0];
            else if(distributionType == PROP_VARIABLE_Z)
                return values[layer - sliceMin];
            if(!insideGrid(row, col)) return nullValue;
            if(distributionType == PROP_VARIABLE_XY)
            {
                return values[getIndex2d(row, col)];
            }
            else //  if(distributionType == PROP_VARIABLE_XYZ)
                return values[getIndex3d(row, col, layer)];
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getValue", "row " + row + " col " + col + " layer " + layer, e);
            return nullValue;
        }
    }

    public float getBlockValue(int blockRow, int blockCol, int layer)
    {
        int index;
        try
        {
            if(distributionType == PROP_CONSTANT)
                index = 0;
            else if(distributionType == PROP_VARIABLE_Z)
                index = getIndexZ(layer);
            else if(distributionType == PROP_VARIABLE_XY)
                index = getBlockIndex2d(blockRow, blockCol);
            else if(distributionType == PROP_VARIABLE_XYZ)
                index = getBlockIndex3d(blockRow, blockCol, layer);
            else
            {
                StsException.systemError(this, "getBlockValue", "Unknown distribution type: " + distributionType);
                return nullValue;
            }
            return getBlockValue(index);
        }
        catch (Exception e)
        {
            StsException.outputWarningException(this, "getValue", "blockRow " + blockRow + " blockCol " + blockCol + " layer " + layer, e);
            return nullValue;
        }
    }

    private int getIndexZ(int layer)
    {
        return layer - sliceMin;
    }

    public float getBlockValue(int index)
    {
        if(!checkLoadData())
            return nullValue;
        if(index < 0 || index >= values.length)
            return nullValue;
        else
            return values[index];
    }

    public void setBlockValue(int[] blockIJK, float value)
    {
        setBlockValue(blockIJK[0], blockIJK[1], blockIJK[2], value);
    }
    public void setBlockValue(int blockRow, int blockCol, int slice, float value)
    {
        if(StsBlock.debugIJK(blockRow + rowMin, blockCol + colMin, slice))
            StsException.systemDebug(this, "setBlockValue",  "volume " + getName() + " row " + (blockRow + rowMin) + " col " + (blockCol + colMin) + " layer " + slice + " value: " + value);
        int index = getBlockIndex3d(blockRow, blockCol, slice);
        values[index] = value;
        adjustRange(value);
    }

    public void computeValueRange()
    {
        if(values == null) return;
        for(int n = 0; n < values.length; n++)
        {
            if(values[n] == nullValue) continue;
            valueMin = Math.min(valueMin, values[n]);
            valueMax = Math.max(valueMax, values[n]);
        }
    }

    private void adjustRange(float value)
    {
        if(value == nullValue) return;
        if(value < valueMin) valueMin = value;
        if(value > valueMax) valueMax = value;
    }

    protected void setRange()
    {
        if(values == null) return;
        for(int n = 0; n < values.length; n++)
        {
            float value = values[n];
            if(value == nullValue) continue;
            if(value < valueMin)
                valueMin = value;
            if(value > valueMax)
                valueMax = value;
        }

    }

    public void setValue(int row, int col, int slice, float value)
    {
        int index = getIndex3d(row, col, slice);
        if(index == 37684)
            StsException.systemDebug(this, "setBlockValue");
        values[index] = value;
        adjustRange(value);
     }

    public void addBlockValue(int blockRow, int blockCol, int slice, float value)
    {
        int index = getBlockIndex3d(blockRow, blockCol, slice);
        values[index] += value;
        adjustRange(value);
    }

    public void addValue(int row, int col, int slice, float value)
    {
        int index = getIndex3d(row, col, slice);
        values[index] += value;
        adjustRange(value);
    }

    public void setValues(StsPropertyVolume inputVolume, int topLayerNumber, int nLayers)
    {
        if(inputVolume == null) return;
        
        if(inputVolume.distributionType == distributionType)
        {
            if(inputVolume.sameAs(this))
                values = inputVolume.values;
            else
            {
                if(distributionType == PROP_CONSTANT)
                    values = inputVolume.values;
                else if(distributionType == PROP_VARIABLE_Z)
                    System.arraycopy(inputVolume.values, 0, values, topLayerNumber, nLayers);
                else if(distributionType == PROP_VARIABLE_XY)
                {
                    int nRows = getNRows();
                    int nCols = getNCols();
                    int nInputRows = inputVolume.getNRows();
                    int nInputCols = inputVolume.getNCols();
                    if(nInputRows == nRows && nInputCols == nCols)
                        values = inputVolume.values;
                    else
                        StsException.systemError(this, "setValues", "Haven't implemented different XY sizes yet.");
                }
                else if(distributionType == PROP_VARIABLE_XYZ)
                {
                    int nRows = getNRows();
                    int nCols = getNCols();
                    int nInputRows = inputVolume.getNRows();
                    int nInputCols = inputVolume.getNCols();
                    if(nInputRows == nRows && nInputCols == nCols)
                        System.arraycopy(inputVolume.values, 0, values, topLayerNumber*nRows*nCols, nRows*nCols*nLayers);
                    else
                        StsException.systemError(this, "setValues", "Haven't implemented different XYZ sizes yet.");
                }
                valueMin = Math.min(valueMin, inputVolume.valueMin);
                valueMax = Math.max(valueMax, inputVolume.valueMax);
            }
        }
    }

    public StsColor getScaledColor(float value)
    {
        return propertyType.getColor(value);
    }

    public float getValue(float xValue, float yValue, int nLayer)
    {
        int col = (int)((xValue-xMin)/xInc);
        int row = (int)((yValue-yMin)/yInc);
        return getValue(row, col, nLayer);
    }

    public boolean insideBlockGrid(int row, int col)
    {
        return row >= 0 && row < getNRows() && col >= 0 && col < getNCols();
    }

    public boolean insideGrid(int row, int col)
    {
        return row >= rowMin && row <= rowMax && col >= colMin && col <= colMax;
    }

    public StsFloatVector computeFloatVector()
    {
        StsFloatVector floatVector = new StsFloatVector(values);
        floatVector.setMinValue(valueMin);
        floatVector.setMaxValue(valueMax);
        values = null;
        return floatVector;
    }

    public void outputDataFile(StsFloatVector floatVector)
    {
        values = floatVector.getValues();
        StsFile file = StsFile.constructor(filename);
        if(file.exists()) file.delete();
        file.createNewFile();
        StsBinaryFile binaryFile = new StsBinaryFile(file);
        binaryFile.openWrite();
        binaryFile.setByteValue(distributionType);
        binaryFile.setFloatVector(floatVector);
        binaryFile.close();
        values = null;
    }

    public boolean loadDataFile()
    {
        StsFile file = StsFile.constructor(filename);
        if(!file.exists())
        {
            filename = replaceDriveLetterHack("E", "F");
            file = StsFile.constructor(filename);
        }
        if(!file.exists())
        {
            return false;
        }
        StsBinaryFile binaryFile = new StsBinaryFile(file);
        binaryFile.openRead();
        distributionType = binaryFile.getByteValue();
        int nValues = getNPoints();
        StsFloatVector floatVector = new StsFloatVector(nValues);
        binaryFile.getFloatVector(floatVector, true);
        binaryFile.close();
        values = floatVector.getValues();
        // valueMin = floatVector.getMinValue();
        // valueMax = floatVector.getMaxValue();
        // propertyType.setColorscaleRange(valueMin, valueMax);
        return true;
    }

    private String replaceDriveLetterHack(String oldLetter, String newLetter)
    {
        return filename.replaceFirst(oldLetter + ":", newLetter + ":");
    }

    public StsColor getColor(int row, int col, int layer)
    {
        if(!checkLoadData()) return null;
        if(distributionType == PROP_NONE)
            return StsColor.GREY;
        else
        {
            int index = getDistributionIndex(row, col, layer);
            if(index < 0 || index >= values.length)
                return StsColor.BLACK;
            float value = values[index];
            if(value == nullValue)
                return StsColor.GREY;
            return propertyType.getColor(values[index]);
        }
    }

    private boolean checkLoadData()
    {
        return values != null || loadDataFile();
    }

    protected int getDistributionIndex(int row, int col, int layer)
    {
        switch(distributionType)
        {
            case PROP_CONSTANT:
                return 0;
            case PROP_VARIABLE_Z:
                return getIndexZ(layer);
            case PROP_VARIABLE_XY:
                return getIndex2d(row, col);
            case PROP_VARIABLE_XYZ:
                return getIndex3d(row, col, layer);
            default:
                return 0;
        }
    }

    public boolean matches(StsBlock block, StsPropertyType propertyType)
    {
        return boundingBox == block && this.propertyType.eclipseName.equals(propertyType.eclipseName);
    }

    public String toString()
    {
        return name + " boundingBox " + boundingBox.name + " rowMin " + rowMin + " rowMax " + rowMax + " colMin " + colMin + " colMax: " + colMax + " sliceMin " + sliceMin + " sliceMax " + sliceMax;
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

    public byte getDistributionType()
    {
        return distributionType;
    }
}