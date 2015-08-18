//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.event.*;

public class StsAttributeSurface extends StsSurface implements StsSelectable, StsXYSurfaceGridable, StsTreeObjectI,
                  ItemListener, StsTextureSurfaceFace, StsDistanceTransformFace, StsCultureDisplayable
{
    /** colorscale for depth display */
    protected StsColorscale[] attributeColorscale = null;
    protected String[] attributeNames = null;
    protected float[][][] attributes = null;

    /** DB constructor */
    public StsAttributeSurface()
    {
    }
    
    static public StsAttributeSurface constructSurface(String name, StsColor stsColor, byte type, int nCols, int nRows,
                             double xOrigin, double yOrigin, float xInc, float yInc, float xMin, float yMin, 
                             float angle, float[][] pointsZ, String[] attributeNames, float[][][] attributes,
                             boolean hasNulls, float nullZValue, byte zDomain, byte vUnits, byte hUnits, StsProgressPanel progressPanel)
    {
        try
        {
            return new StsAttributeSurface(name, stsColor, type, nCols, nRows, xOrigin,
                                  yOrigin, xInc, yInc, xMin, yMin, angle,
                                  pointsZ, attributeNames, attributes, hasNulls, nullZValue, zDomain, progressPanel);
        }
        catch (Exception e)
        {
            StsException.outputException("StsAttributeSurface.constructSurface() failed.", e, StsException.WARNING);
            return null;
        }
    }
    
    private StsAttributeSurface(String name, StsColor stsColor, byte type, int nCols, int nRows, double xOrigin, 
    		              double yOrigin, float xInc, float yInc, float xMin, float yMin, float angle, 
    		              float[][] pointsZ, String[] attributeNames, float[][][] attributes, boolean hasNulls, float nullZValue,
    		              byte zDomain, StsProgressPanel progressPanel) throws StsException
    {
    	initialize(name, stsColor, type, nCols, nRows, xOrigin, yOrigin, xInc,
    			yInc, xMin, yMin, angle, pointsZ, hasNulls, nullZValue,
    			zDomain, zDomain, progressPanel);
        initSurfaceTextureList();
    	this.attributeNames = attributeNames;
    	this.attributes = attributes;
    	this.attributeColorscale = new StsColorscale[attributes.length];
    }
    
    public boolean initialize(StsModel model)
    {
    	super.initialize(model);
        return true;
    }

    /** if surface is deleted, it's index will now be -1; classes display this texture should delete it from their list */
    public boolean isDisplayable()
    {
        return isPersistent();
    }

    public StsColorscale getAttributeColorscale(String name, float min, float max)
    {
        int attIndex = getAttributeIndex(name);  
        if(attIndex == -1)
        	StsMessageFiles.infoMessage("Unable to locate colorscale for selected attribute (" + name + ").");
        if (attributeColorscale[attIndex] == null)
        	createAttributeColorscale(name, min, max);
        return attributeColorscale[attIndex];
    }
    
    public void createAttributeColorscale(String name, float min, float max)
    {
        double[] scale = StsMath.niceScale(min, max, 32, true);
        float smin = (float) scale[0];
        float smax = (float) scale[1];

        StsSpectrum spectrum = currentModel.getSpectrum(StsSpectrumClass.SPECTRUM_RAINBOW);
        int attIndex = getAttributeIndex(name);
        if(attIndex == -1)
        	return;
        attributeColorscale[attIndex] = new StsColorscale(name, spectrum, smin, smax);
        attributeColorscale[attIndex].addActionListener(this);
    }
    
    public float[] getAttributeHistogram()
    {
    	return StsToolkit.buildHistogram(getTextureData(), getZMin(), getZMax());
    }
    
    public int getAttributeIndex(String name)
    {
    	for(int i=0; i<attributes.length; i++)
    	{
    		if(attributeNames[i].equals(name))
    			return i;
    	}
    	return -1;
    }
}