
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.DBTypes;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

public class StsSubZone extends StsObject implements StsSelectable
{
    protected String name;
    protected StsColor stsColor;
    protected StsZone zone;
    protected int subZoneNumber;
    protected float topF, botF;  // fractional position of top and bottom between surfaces  (SUBZONE_UNIFORM)
                                 // or actual distances from top or bottom (other layering styles)
    protected StsModelSurface botSurface, topSurface;
    protected int status = UNDEFINED;
    protected int nLayers = 1;
    protected int subZoneType = SUBZONE_UNIFORM;
    protected float layerThickness = 0.0f;
    protected StsWellZoneSet wellSubZoneSet = null;

    // *** comment out types not yet implemented! ***
    static public final String[] subZoneTypeStrings =
    {
    	"Stratigraphic",
        "Onlap",
        "Offlap"
    };

    static public StsFieldBean[] displayFields = null;
    static public StsFieldBean[] propertyFields = null;

    // constants
    static public final int SUBZONE_UNIFORM = 0;
    static public final int SUBZONE_ONLAP = 1;
    static public final int SUBZONE_OFFLAP = 2;
    static public final int SUBZONE_TRUNCATED = 3;

    static public final int UNDEFINED = StsParameters.UNDEFINED;

    public StsSubZone()
    {
    }

    public StsSubZone(int subZoneType, float topF, float botF, StsZone zone, int subZoneNumber)
    {
        this(subZoneType, topF, botF, zone, subZoneNumber, false);
    }

    public StsSubZone(int subZoneType, float topF, float botF, StsZone zone, int subZoneNumber, boolean persistent)
    {
        super(persistent);
        this.subZoneType = subZoneType;
        this.topF = topF;
        this.botF = botF;
        this.zone = zone;
        this.subZoneNumber = subZoneNumber;
        this.name = new String(zone.getName() + "-" + subZoneNumber);
        currentModel.incrementSpectrumColor("Basic");
        stsColor = currentModel.getCurrentSpectrumColor("Basic");
    }

    public StsSubZone(StsModelSurface botSurface, StsModelSurface topSurface, StsZone zone, int subZoneNumber)
    {
        this.botSurface = botSurface;
        this.topSurface = topSurface;
        this.zone = zone;
        this.subZoneNumber = subZoneNumber;
        this.name = new String(zone.getName() + "-" + subZoneNumber);
        currentModel.incrementSpectrumColor("Basic");
        stsColor = currentModel.getCurrentSpectrumColor("Basic");
    }

    public boolean initialize(StsModel model) { return true; }

    public String getName() { return name; }

    public void setNLayers(int nLayers)
    {
        if(this.nLayers == nLayers) return;
        this.nLayers = nLayers;
        zone.subZoneChanged(this);
    }

    public int getNLayers() { return nLayers; }
    public void setSubZoneType(int subZoneType)
    {
        if(this.subZoneType == subZoneType) return;
        this.subZoneType = subZoneType;
        zone.subZoneChanged(this);
    }

    public int getSubZoneType() { return subZoneType; }
    public void setLayerThickness(float layerThickness) { this.layerThickness = layerThickness; }
    public float getLayerThickness() { return layerThickness; }
    public int getSubZoneNumber() { return subZoneNumber; }
    public StsColor getStsColor() { return stsColor; }

    public StsFieldBean[] getDisplayFields()
    {
        if(displayFields != null) return displayFields;
        displayFields = new StsFieldBean[]
        {
            new StsColorComboBoxFieldBean(StsSubZone.class, "stsColor", "Color")
        };
        return displayFields;
    }

    public StsFieldBean[] getPropertyFields()
    {
        if(propertyFields == null)
        {
            propertyFields = new StsFieldBean[]
            {
                new StsIntFieldBean(StsSubZone.class, "nLayers", "Layers"),
                new StsIntFieldBean(StsSubZone.class, "subZoneType", "subUnitType"),
                new StsFloatFieldBean(StsSubZone.class, "layerThickness", "Thickness")
            };
        }
        return propertyFields;
    }
}
