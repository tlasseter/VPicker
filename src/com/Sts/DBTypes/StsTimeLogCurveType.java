package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

import javax.media.opengl.*;
import java.awt.event.*;
import java.util.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Feb 25, 2009
 * Time: 7:37:03 AM
 * To change this template use File | Settings | File Templates.
 */
public class StsTimeLogCurveType extends StsLogCurveType
{
    protected StsColorscale colorscale;
    transient private StsSpectrum spectrum;
    transient public StsColorList colorList = null;

    static StsEditableColorscaleFieldBean colorscaleBean;
    
    public StsTimeLogCurveType()
    {
    }

    public StsTimeLogCurveType(StsTimeLogCurve logCurve)
    {
        super(logCurve);
    }

    public StsTimeLogCurveType(String name, float min, float max, StsSpectrum spectrum)
    {
        super(name, min, max);
        this.spectrum = spectrum;
        initializeColorscale();
    }

    private void initialize()
    {
        logCurveTypeClass = (StsTimeLogCurveTypeClass)getStsClass();
    }

    public boolean initialize(StsModel model)
    {
        initialize();
        this.reinitializeScale();
        return true;
    }

    static public String[] getLogCurveTypeStrings()
    {
        StsObject[] logCurveTypes = currentModel.getObjectList(StsTimeLogCurveType.class);
        int nTypes = logCurveTypes.length;
        String[] names = new String[nTypes];
        for(int n = 0; n < nTypes; n++)
        {
            StsLogCurveType logCurveType = (StsLogCurveType)logCurveTypes[n];
            names[n] = logCurveType.getName();
        }
        return names;
    }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields == null)
        {
            colorscaleBean = new StsEditableColorscaleFieldBean(StsTimeLogCurveType.class, "colorscale");

            displayFields = new StsFieldBean[]
            {
                new StsStringFieldBean(StsTimeLogCurveType.class, "name", false, "Name:"),
                new StsComboBoxFieldBean(StsTimeLogCurveType.class, "aliasToString", "Alias to:", "logCurveTypeStrings"),
                new StsFloatFieldBean(StsTimeLogCurveType.class, "displayCurveMin", true, "Display min:"),
                new StsFloatFieldBean(StsTimeLogCurveType.class, "displayCurveMax", true, "Display max:"),
                new StsFloatFieldBean(StsTimeLogCurveType.class, "curveMin", false, "Curve min:"),
                new StsFloatFieldBean(StsTimeLogCurveType.class, "curveMax", false, "Curve max:"),
                new StsColorComboBoxFieldBean(StsTimeLogCurveType.class, "stsColor", "Color:", StsLogCurveTypeClass.logCurveTypeColors),
                new StsComboBoxFieldBean(StsTimeLogCurveType.class, "scaleTypeString", "Curve Type:", scaleTypeStrings),
                colorscaleBean
            };
        }
        return displayFields;
    }

    public void initializeColorscale()
    {
        try
        {
            if(colorscale == null)
            {
                constructColorscale();
                dbFieldChanged("colorscale", colorscale);
            }
            colorList = new StsColorList(colorscale);
            colorscale.addActionListener(this);
        }
        catch (Exception e)
        {
            StsException.outputException("StsTimeLogCurve.initializeColorscale() failed.", e, StsException.WARNING);
        }
    }

    private void constructColorscale()
    {
        colorscale = new StsColorscale(name, spectrum, curveMin, curveMax);
    }

    public boolean setGLColorList(GL gl, boolean nullsFilled, int shader)
	{
		if(colorList == null) initializeColorscale();
	    return colorList.setGLColorList(gl, nullsFilled, shader);
	}
    
    private void colorscaleChanged()
    {
        Iterator iter = logCurves.getVisibleObjectIterator();
        while(iter.hasNext())
        {
            StsTimeLogCurve timeLogCurve = (StsTimeLogCurve)iter.next();
            timeLogCurve.colorscaleChanged();
        }
        if(colorList != null) colorList.setColorListChanged(true);
        currentModel.viewObjectChanged(this, this);
    }

    public void actionPerformed(ActionEvent e)
    {
        if (e.getSource() instanceof StsColorscale)
        {
            colorscaleChanged();
            currentModel.viewObjectRepaint(this, this);
        }
    }

    public StsColorscale getColorscale() { return colorscale; }
    public void setColorscale(StsColorscale colorscale) { this.colorscale = colorscale; }

    public void setDisplayCurveMin(float min)
    {
        colorscale.setEditMin(min);
        colorscaleChanged();
        super.setDisplayCurveMin(min);
    }

    public void setDisplayCurveMax(float max)
    {
        colorscale.setEditMax(max);
        colorscaleChanged();
        super.setDisplayCurveMax(max);
    }
}
