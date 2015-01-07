package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.nio.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 2, 2009
 * Time: 10:14:39 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class StsDisplaySeismicAttribute implements ActionListener
{
    protected StsModel model;
    StsSeismicBoundingBox seismicBoundingBox;
    public StsColorscale colorscale;
    StsColorList colorList;
    StsTextureTiles textureTiles = null;
    boolean textureChanged = true;
    public ByteBuffer byteBuffer = null;

    abstract public boolean computeBytes();
    abstract public String getName();

    public StsColorscale getColorscaleWithName(String name)
    {
        StsException.notImplemented(this, "getColorscaleWithName");
        return null;
    }

    public void addColorscale(StsColorscale colorscale)
    {
        StsException.notImplemented(this, "addColorscale");
    }

    public String toString() { return getName(); }


    public void initialize(StsSeismicBoundingBox seismicBoundingBox, StsModel model)
	{
        this.model = model;
        this.seismicBoundingBox = seismicBoundingBox;
    }

    public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
		if(colorList == null)
			return false;
		else
		    return colorList.setGLColorList( gl, nullsFilled, shader);
	}

    protected void initializeColorscale()
	{
        if(colorscale == null) return;
        colorscale.addActionListener(this);
		colorList = new StsColorList(colorscale);
    }

    public void actionPerformed(ActionEvent e)
	{
		if (! (e.getSource() instanceof StsColorscale)) return;
		colorListChanged();
        model.viewObjectRepaint(this, colorscale);
	}

    public void colorListChanged()
    {
        if(colorList != null) colorList.setColorListChanged(true);
        byteBuffer = null;
    }

    public StsColorscale getInitializeColorscale(String spectrumName, float rangeMin, float rangeMax)
    {
        if(colorscale != null) return colorscale;
        // if we have reloaded an existing db, we may already have this colorscale
        colorscale = getColorscaleWithName(getName());
        // if we don't have this colorscale, build it
        if(colorscale == null)
        {
            StsSpectrum spectrum = model.getSpectrum(spectrumName);
            colorscale = new StsColorscale(getName(), spectrum, rangeMin, rangeMax);
            // add it to the lineset colorscales; it will now be persisted
            addColorscale(colorscale);
        }
        // having recovered or built this colorscale, initialize it
        initializeColorscale();
        return colorscale;
    }
}
