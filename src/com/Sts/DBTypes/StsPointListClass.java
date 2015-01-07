package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.DB.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.View3d.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

import java.util.*;

public class StsPointListClass extends StsClass implements StsSerializable, StsClassDisplayable, StsClassCursorDisplayable
{
    protected boolean displayProperties = false;

    private String defaultSpectrumName = StsSpectrumClass.SPECTRUM_RWB;
    private int defaultNumberBins = 10;

    public StsPointListClass()
    {
    }

    public void initializeDefaultFields()
    {
        defaultFields = new StsFieldBean[]
        {
            new StsIntFieldBean(this, "defaultNumberBins", 2, 50, "Amplitude Scale:"),
            new StsComboBoxFieldBean(this, "defaultSpectrumName", "Spectrum:", StsSpectrumClass.cannedSpectrums)
        };
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        int nPointLists = getSize();
        for(int n = 0; n < nPointLists; n++)
        {
            StsPointList pl = (StsPointList)getElement(n);
            pl.display(glPanel3d);
        }
    }

    public void drawOnCursor3d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, StsPoint[] planePoints, boolean isDragging)
    {
    }

    public void drawOn3dCurtain(StsGLPanel3d glPanel3d, StsSeismicCurtain seismicCurtain)
    {
    }

    public void drawOnCursor2d(StsGLPanel3d glPanel3d, int dirNo, float dirCoordinate, boolean axesFlipped, boolean xAxisReversed, boolean yAxisReversed)
    {
        Iterator iter = getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsPointList pset = (StsPointList)iter.next();
            pset.display2d(glPanel3d, dirNo, dirCoordinate, axesFlipped, xAxisReversed, yAxisReversed);
        }
    }

    public StsPointList getPointSet(String name) { return (StsPointList)getObjectWithName(name); }

    public void setDisplayProperties(boolean b) { displayProperties = b; }
    public boolean getDisplayProperties() { return displayProperties; }
    public boolean setCurrentObject(StsObject object)
    {
        boolean changed = super.setCurrentObject(object);
        if (changed)
            if (object != null) ( (StsPointList) object).treeObjectSelected();
        return changed;
    }

    public String getDefaultSpectrumName()  { return defaultSpectrumName; }

    public void setDefaultSpectrumName(String value)
    {
        if(this.defaultSpectrumName.equals(value)) return;
        this.defaultSpectrumName = value;
//        setDisplayField("defaultSpectrumName", defaultSpectrumName);
    }

    public void setDefaultNumberBins(int value)
    {
        if(this.defaultNumberBins == value) return;
        this.defaultNumberBins = value;
//        setDisplayField("defaultNumberBins", defaultNumberBins);
    }

    public int getDefaultNumberBins()  { return defaultNumberBins; }
}
