package com.Sts.DBTypes;

import com.Sts.DB.*;
import com.Sts.Utilities.*;
import com.Sts.Types.*;
import com.Sts.MVC.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Apr 27, 2010
 * Time: 9:26:52 AM
 * To change this template use File | Settings | File Templates.
 */

public class StsPropertyType extends StsMainObject
{
    public String eclipseName;
    protected StsColorscale colorscale;
    private float valueMin = StsParameters.largeFloat;
    private float valueMax = -StsParameters.largeFloat;

    static public final StsPropertyTypeDefinition none = new StsPropertyTypeDefinition("None", "NONE");
    static public final StsPropertyTypeDefinition porosity = new StsPropertyTypeDefinition ("Porosity", "PORO", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition permX = new StsPropertyTypeDefinition ("Perm X", "PERMX", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition permY = new StsPropertyTypeDefinition ("Perm Y", "PERMY", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition permZ = new StsPropertyTypeDefinition ("Perm Z", "PERMZ", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition poreVolume = new StsPropertyTypeDefinition ("Pore Volume", "PORV", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition actnum = new StsPropertyTypeDefinition ("Active Flag", "ACTNUM", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition depth = new StsPropertyTypeDefinition ("Depth", "DEPTH", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition tranX = new StsPropertyTypeDefinition ("Tran X", "TRANX", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition tranY = new StsPropertyTypeDefinition ("Tran Y", "TRANY", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition tranZ = new StsPropertyTypeDefinition ("Tran Z", "TRANZ", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition layerColor = new StsPropertyTypeDefinition ("Layer Color", "LAYER_COLOR");
    static public final StsPropertyTypeDefinition pressure = new StsPropertyTypeDefinition ("Pressure", "PRESSURE", StsSpectrumClass.SPECTRUM_SEMBLANCE);
    static public final StsPropertyTypeDefinition wsat = new StsPropertyTypeDefinition ("Water Sat", "SWAT", StsSpectrumClass.SPECTRUM_WATER_SAT);
    static public final StsPropertyTypeDefinition indexMap = new StsPropertyTypeDefinition ("S2S->ECL index map", "INDEX_MAP");

    static public final StsPropertyTypeDefinition[] propertyTypeNames = new StsPropertyTypeDefinition[] { none, porosity, permX, permY, permZ, poreVolume,
        actnum, depth, tranX, tranY, tranZ, layerColor, pressure, wsat, indexMap };

    static public StsPropertyType propertyTypeNone = new StsPropertyType(none, false);

    public StsPropertyType()
    {
    }
    
    public StsPropertyType(StsPropertyTypeDefinition propertyTypeName)
    {
        this(propertyTypeName, true);
    }

    public StsPropertyType(StsPropertyTypeDefinition propertyTypeName, boolean persistent)
    {
        super(persistent);
        setName(propertyTypeName.name);
        this.eclipseName = propertyTypeName.eclipseName;
        if(propertyTypeName.spectrumName == null) return;
        StsSpectrumClass spectrumClass = (StsSpectrumClass)currentModel.getStsClass(StsSpectrum.class);
        StsSpectrum propertySpectrum = spectrumClass.getSpectrum(propertyTypeName.spectrumName);
        colorscale = new StsColorscale(propertySpectrum);
        if(persistent) addToModel();
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    static public StsPropertyType constructor(StsPropertyTypeDefinition propertyTypeDefinition)
    {
        return new StsPropertyType(propertyTypeDefinition);
    }

    static public StsPropertyType constructor(StsPropertyTypeDefinition propertyTypeDefinition, boolean persistent)
    {
        return new StsPropertyType(propertyTypeDefinition, persistent);
    }

    public void adjustRange(float valueMin, float valueMax)
    {
        this.valueMin = Math.min(this.valueMin, valueMin);
        this.valueMax = Math.max(this.valueMax, valueMax);
    }

    public void rangeChanged()
    {
        dbFieldChanged("valueMin", valueMin);
        dbFieldChanged("valueMax", valueMax);
        if(colorscale != null)
        {
            colorscale.setRange(valueMin, valueMax);
            colorscale.dbFieldChanged("editRange", colorscale.editRange);
            colorscale.dbFieldChanged("range", colorscale.range);
        }
    }

    public void setColorscaleRange(float valueMin, float valueMax)
    {
        if(colorscale == null) return;
        colorscale.setRange(valueMin, valueMax);
    }

    protected StsColor getColor(float value)
    {
        return colorscale.getStsColorFromValue(value);
    }

    public float getValueMin()
    {
        return valueMin;
    }

    public void setValueMin(float valueMin)
    {
        this.valueMin = valueMin;
    }

    public float getValueMax()
    {
        return valueMax;
    }

    public void setValueMax(float valueMax)
    {
        this.valueMax = valueMax;
    }
}