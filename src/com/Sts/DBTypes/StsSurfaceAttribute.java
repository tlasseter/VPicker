//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.IO.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.*;
import java.nio.*;

public class StsSurfaceAttribute extends StsMainObject
{
	public StsSurface surface = null;
	protected StsColorscale colorscale = null;
    public float dataMin = StsParameters.largeFloat;
    public float dataMax = -StsParameters.largeFloat;
    public String binaryFilename = null;
    protected boolean useFilter;

    transient public SurfaceAttributeTexture surfaceTexture;

    private static final long serialVersionUID = 1L;

    /** DB constructor */
    public StsSurfaceAttribute()
    {
    }

    public StsSurfaceAttribute(boolean persistence)
    {
        super(persistence);
    }

    static public StsSurfaceAttribute constructor(StsSurface surface, String name, float aMin, float aMax, float[][] values)
    {
        try
        {
            return new StsSurfaceAttribute(surface, name, aMin, aMax, values);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSurfaceAttribute.constructSurfaceAttribute() failed.", e, StsException.WARNING);
            return null;
        }
    }
    
    static public StsSurfaceAttribute constructor(StsSurface surface, String name, float aMin, float aMax,  float[][] values, StsSurfaceTexture texture)
    {
        try
        {
            return new StsSurfaceAttribute(surface, name, aMin, aMax, values, texture);
        }
        catch (Exception e)
        {
            StsException.outputException("StsSurfaceAttribute.constructSurfaceAttribute() failed.", e, StsException.WARNING);
            return null;
        }
    }

    private StsSurfaceAttribute(StsSurface surface, String name, float dataMin, float aMax,
    		float[][] values, StsSurfaceTexture texture) throws StsException
    {
        setName(name);
        this.surface = surface;
        this.dataMax = aMax;
        this.dataMin = dataMin;
        //System.out.println("Name:" + name + " Min=" + aMin + " Max=" + aMax);
        surfaceTexture = new SurfaceAttributeTexture(this, values);
        surfaceTexture.setTextureData(texture.getTextureData());
    }

    private StsSurfaceAttribute(StsSurface surface, String name, float dataMin, float dataMax, float[][] values) throws StsException
    {
        setName(name);
        this.surface = surface;
        initialize(dataMin, dataMax, values);
    }

    public void initialize(float dataMin, float dataMax, float[][] values)
    {
        this.dataMin = dataMin;
        this.dataMax = dataMax;
        surfaceTexture = new SurfaceAttributeTexture(this, values);
    }

    public void initialize(float dataMin, float dataMax, float[][] values, byte[] textureData)
    {
        this.dataMin = dataMin;
        this.dataMax = dataMax;
        surfaceTexture = new SurfaceAttributeTexture(surface, values, textureData);
    }

    public boolean initialize(StsModel model)
    {
        surfaceTexture = new SurfaceAttributeTexture(this, null);
        return true;
    }

    public String toString()
    {
        return getName();
    }

    public StsColorscale getColorscale()
    {
        if (colorscale == null) createColorscale();
        return colorscale;
    }
    
    public void setColorscale(StsColorscale c)
    {
        colorscale = c;
     }

    public float getDataMin()
    {
        return dataMin;
    }

    public float getDataMax()
    {
        return dataMax;
    }

    public void setDataMin(float value)
    {
        dataMin = value;
        surfaceTexture.setTextureData(null);
    }

    public void setDataMax(float value)
    {
        dataMax = value;
        surfaceTexture.setTextureData(null);
    }
    
    public StsSurfaceTexture getSurfaceTexture()
    {
        if(surfaceTexture == null)
            surfaceTexture = new SurfaceAttributeTexture(this, null);
        return surfaceTexture;
    }
    
    public boolean isFiltered() { return useFilter; }
    public void setFiltered(boolean val) { useFilter = val; }

    public class SurfaceAttributeTexture extends StsSurfaceTexture
    {
        byte[] textureData = null;
        int colorDisplayListNum = 0;
        float[][] values = null;

        public SurfaceAttributeTexture(StsSurfaceAttribute surfaceAttribute, float[][] values)
        {
            super(surfaceAttribute.surface);
            this.values = values;
            initializeColorscaleActionListener();
        }

        public SurfaceAttributeTexture(StsSurface surface, float[][] values, byte[] textureData)
        {
            super(surface);
            this.values = values;
            this.textureData = textureData;
            initializeColorscaleActionListener();
        }

        public void setTextureData(byte[] data) { textureData = data; }

        public byte[] getTextureData()
        {
            if (textureData != null) return textureData;
            checkLoadAttribute();
            byte[][] surfacePointsNull = surface.getPointsNull();
            int nRows = surface.getNRows();
            int nCols = surface.getNCols();
            byte[] textureData = new byte[nRows * nCols];
            float[] editRange = getColorscale().getEditRange();
            float editMin = editRange[0];
            float editMax = editRange[1];
            float scale = StsMath.floatToUnsignedByteScale(editMin, editMax);
            float scaleOffset = StsMath.floatToUnsignedByteScaleOffset(scale, editMin);
            int n = 0;
            if (useFilter)
            {
            	float range = dataMax - dataMin;
                for (int row = 0; row < nRows; row++)
                {
                    for (int col = 0; col < nCols; col++)
                    {
                    	float value = values[row][col] ;
                    	if (value > -StsParameters.largeFloat)
                        {
                   		
    	                    if (value > dataMax) value = dataMax;
    	                    if (value < dataMin) value = dataMin;
                        	float val =  255 *  (value- dataMin) /range;
                        	textureData[n++] = (byte)val;
                        }
                        else
                        {
                            textureData[n++] = -1;
                        }
                    }
                }
            }
            else
            {
	            for (int row = 0; row < nRows; row++)
	            {
	                for (int col = 0; col < nCols; col++)
	                {
	                    float value = values[row][col];
	                    if (value > dataMax) value = dataMax;
	                	if (value < dataMin) value = dataMin;
	                	if (surfacePointsNull[row][col] != StsSurface.SURF_PNT)
	                        textureData[n++] = -1;
	                    else
	                        textureData[n++] = StsMath.floatToUnsignedByte254WithScale(value, scale, scaleOffset);
	                }
	            }
            }
            if(colorscale == null) createColorscale();
            colorscale.setRange(dataMin, dataMax);
            return textureData;
        }

        public void checkLoadAttribute()
        {
            if (values != null) return;
            StsBinaryFile bFile = getBinaryFile(binaryFilename);
            bFile.openRead();
            readBinaryHeader(bFile);
            values = new float[surface.nRows][];
            for (int row = 0; row < surface.nRows; row++)
                values[row] = bFile.getFloatValues();
        }

        // Used to skip the binary header so the data can be loaded.
        private boolean readBinaryHeader(StsBinaryFile binaryFile)
        {
            try
            {
                binaryFile.getByteValues();
                binaryFile.getIntegerValues();
                binaryFile.getDoubleValues();
                binaryFile.getFloatValues();
                binaryFile.getBooleanValues();
                binaryFile.getByteValues();
                return true;
            }
            catch (Exception ex)
            {
                StsException.outputException("StsSurface.writeBinaryFile() failed writing header for surface (" + getName() + ")", ex, StsException.WARNING);
                return false;
            }
        }

        public StsBinaryFile getBinaryFile(String binFilename)
        {
            StsAbstractFile file = null;
            try
		    {
                if(Main.isJarDB)
                {
                    String[] filenames = Main.jar.getFilenames();
                    for(int i=0; i<filenames.length; i++)
                    {
                        if(filenames[i].endsWith(binFilename))
                        {
                            file = Main.jar.getFile(i);
                            break;
                        }
                    }
                    if(file == null) return null;
                }
                else
                {
                    String binaryFilename = StsSurfaceAttribute.currentModel.getProject().getBinaryFullDirString() + binFilename;
                    file = StsFile.constructor(binaryFilename);
                    if (file == null)
                    {
                        StsMessageFiles.infoMessage("Failed to load surface attribute named: " + getName());
                        return null;
                    }
                }
            }
		    catch(Exception e)
		    {
			    StsException.outputException("StsSurfaceAttribute.getBinaryFile() failed.",	e, StsException.WARNING);
			    return null;
		    }
            return new StsBinaryFile(file);
        }

        public String getName()
        {
            return name;
        }

        public void deleteColorDisplayList(GL gl)
        {
            if (colorDisplayListNum == 0)
            {
                return;
            }
            gl.glDeleteLists(colorDisplayListNum, 1);
            colorDisplayListNum = 0;
        }

        public int getColorDisplayListNum(GL gl, boolean nullsFilled)
        {
            deleteColorDisplayList(gl);

            colorDisplayListNum = gl.glGenLists(1);
            if (colorDisplayListNum == 0)
            {
                StsMessageFiles.logMessage(
                    "System Error in StsSurfaceAttribute.getColorListNum(): Failed to allocate a display list");
                return 0;
            }

            gl.glNewList(colorDisplayListNum, GL.GL_COMPILE);
            createColorList(gl, nullsFilled);
            gl.glEndList();

            return colorDisplayListNum;
        }

        private void createColorList(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            int nColors = arrayRGBA[0].length;
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_R, nColors, arrayRGBA[0], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_G, nColors, arrayRGBA[1], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_B, nColors, arrayRGBA[2], 0);
            gl.glPixelMapfv(GL.GL_PIXEL_MAP_I_TO_A, nColors, arrayRGBA[3], 0);
            gl.glPixelTransferf(GL.GL_MAP_COLOR, 1);
            arrayRGBA = null;
        }
/*
        public void createColorTLUT(GL gl, boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            StsJOGLShader.createLoadARBColormap(gl, arrayRGBA);
        }
*/
        public FloatBuffer getComputeColormapBuffer(boolean nullsFilled)
        {
            float[][] arrayRGBA = computeRGBAArray(nullsFilled);
            return StsJOGLShader.computeColormapBuffer(arrayRGBA, 256);
        }
        private float[][] computeRGBAArray(boolean nullsFilled)
        {
            Color[] colors = getColorscale().getNewColorsInclTransparency();
            int nColors = colors.length;
            float[][] arrayRGBA = new float[4][nColors];
            float[] rgba = new float[4];
            for (int n = 0; n < nColors; n++)
            {
                colors[n].getComponents(rgba);
                for (int i = 0; i < 4; i++)
                {
                    arrayRGBA[i][n] = rgba[i];
                }
            }
            return arrayRGBA;
        }

        public StsColorscale getColorscale()
        {
            return StsSurfaceAttribute.this.getColorscale();
        }


        public float[] getHistogram()
        {
            return StsToolkit.buildHistogram(getTextureData(), dataMin, dataMax);
        }


        public boolean isDisplayable()
        {
            return surface.isPersistent();
        }

        public void initializeColorscaleActionListener()
        {
            StsColorscale colorscale = getColorscale();
            colorscale.addActionListener(surface);
        }

        public boolean writeBinaryFile(String filename)
        {
            StsBinaryFile binaryFile;
            try
            {
                binaryFile = getBinaryFile(filename);
                binaryFile.openWrite(false);
                surface.writeBinaryHeader(binaryFile);

                for (int row = 0; row < surface.nRows; row++)
                    binaryFile.setFloatValues(values[row], surface.pointsNull[row], StsSurface.SURF_PNT, surface.nullZValue);
                binaryFile.close();
                return true;
            }
            catch (Exception e)
            {
                StsException.outputException("StsSurface.writeBinaryFile() failed.", e, StsException.WARNING);
                return false;
            }
        }

        public String toString()
        {
            return StsSurfaceAttribute.this.getName();
        }

        public void selected()
        {
            attributeSelected();
        }
    }

    protected void attributeSelected() { }

    protected void createColorscale()
    {
        StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);

        //double[] scale = StsMath.niceScale( (double)aMin,(double)aMax, 32, true);
        double[] scale = StsMath.niceScale((double) dataMin, (double) dataMax, spectrum.getNColors(), true);
        dataMin = (float) scale[0];
        dataMax = (float) scale[1];
        System.out.println("Scaled(" + getName() + "): Min=" + dataMin + " Max=" + dataMax);

        colorscale = new StsColorscale(name, spectrum, dataMin, dataMax);
        colorscale.addActionListener(surface);
    }

    public boolean writeBinaryFile(String filename)
    {
        this.binaryFilename = filename;
        return surfaceTexture.writeBinaryFile(filename);
    }
}