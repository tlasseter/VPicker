package com.Sts.DBTypes;


import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;
import com.Sts.Utilities.*;

public class StsFaultLine extends StsLine implements StsTreeObjectI, StsTimeDepthDisplayable
{
    transient protected StsList faultZones = null;
    transient public StsFaultStickSet faultStickSet = null;
    // display fields
    static public StsFieldBean[] faultDisplayFields = null;

    public StsFaultLine()
    {
    }

    public StsFaultLine(boolean persistent)
    {
        super(persistent);
    }

    static public StsFaultLine buildFault()
    {
        StsFaultLine faultLine = new StsFaultLine();
        faultLine.setZDomainOriginal(currentModel.getProject().getZDomain());
        faultLine.initializeVertices();
        return faultLine;
    }

    static public StsFaultLine buildImportedFault()
    {
        StsFaultLine faultLine = buildFault();
        faultLine.setVerticesRotated(false);
        return faultLine;
    }

    static public StsFaultLine buildImportedFault(StsFaultStickSet faultStickSet)
    {
        StsFaultLine faultLine = buildImportedFault();
        faultLine.faultStickSet = faultStickSet;
        return faultLine;
    }

    static public StsFaultLine buildVerticalFault(StsGridPoint gridPoint)
    {
        try
        {
            StsFaultLine line = StsFaultLine.buildFault();
            line.constructVertical(gridPoint, StsParameters.FAULT);
            return line;
        }
        catch (Exception e)
        {
            StsException.systemError(StsFaultLine.class, "buildVerticalFault");
            return null;
        }
    }

    static public StsFaultLine buildFault(StsSurfaceVertex[] vertices)
    {
        try
        {
            StsFaultLine line = StsFaultLine.buildFault();
            line.construct(vertices);
            addMDepthToVertices(vertices);
            return line;
        }
        catch (Exception e)
        {
            StsException.systemError("StsLine.buildFault(gridPoints) failed.");
            return null;
        }
    }


    public boolean initialize(StsModel model)
    {
        if (initialized) return true;
        if (!extendEnds()) return false;
        initialized = initializeSection();
        checkAdjustFromVelocityModel();
        return initialized;
    }

    public void deleteTransients()
    {
        zones = null;
        faultZones = null;
    }

    public void addToSet(StsFaultStickSet stickSet)
    {
        faultStickSet = stickSet;
        faultStickSet.faultSticks.add(this);
    }
    
    public String getName()
    {
        if (name != null)
            return name;
        else
            return "FaultLine-" + getIndex();
    }

    public String getLabel()
    {
        return getName() + " ";
    }

    public String toString()
    {
        return getName();
    }

    public StsLineZone getLineZone(StsSurfaceVertex topVertex, StsSurfaceVertex botVertex, StsSection section, int side)
    {
        boolean isFault;
        StsLineZone zone;

        if (side == RIGHT)
            zone = getLineZone(topVertex, botVertex, section, isFault = false);
        else
            zone = getLineZone(topVertex, botVertex, section, isFault = true);

        if (zone == null)
        {
            StsException.outputException(new StsException(StsException.WARNING,
                "StsLine.getLineZone() failed.",
                "Tried to construct a zone for vertices: " + topVertex.getLabel() +
                    " and " + botVertex.getLabel() + " on: " + getLabel()));
        }
        return zone;
    }

    public boolean delete()
    {
        if (!super.delete()) return false;
        return true;
    }

    public StsLineZone getLineZone(StsSurfaceVertex topVertex, StsSurfaceVertex botVertex, StsSection section, boolean isFault)
    {
        StsList lineZones;
        StsLineZone zone;

        if (isFault)
        {
            if (faultZones == null) faultZones = new StsList(2, 2);
            lineZones = faultZones;
        }
        else
        {
            if (zones == null) zones = new StsList(2, 2);
            lineZones = zones;
        }

        int nZones = lineZones.getSize();
        for (int n = 0; n < nZones; n++)
        {
            zone = (StsLineZone) lineZones.getElement(n);
            if (zone.getTop() == topVertex) return zone;
        }

        /** Couldn't find zone: make one and store it in zones */
        zone = new StsLineZone(currentModel, topVertex, botVertex, this);
        lineZones.add(zone);
        return zone;
    }

    public StsFieldBean[] getDisplayFields()
    {
        if (faultDisplayFields == null)
        {
            faultDisplayFields = new StsFieldBean[]
                {
                    new StsBooleanFieldBean(StsFaultLine.class, "isVisible", "Enable"),
                    new StsBooleanFieldBean(StsFaultLine.class, "drawLabels", "Draw Labels"),
                    new StsColorComboBoxFieldBean(StsFaultLine.class, "stsColor", "Color", currentModel.getSpectrum("Basic").getStsColors()),
                    new StsFloatFieldBean(StsFaultLine.class, "topZ", false, "Min Depth"),
                    new StsFloatFieldBean(StsFaultLine.class, "botZ", false, "Max Depth"),
                    new StsDoubleFieldBean(StsFaultLine.class, "xOrigin", false, "X Origin"),
                    new StsDoubleFieldBean(StsFaultLine.class, "yOrigin", false, "Y Origin")
                };
        }
        return faultDisplayFields;
    }

    public StsFieldBean[] getPropertyFields() { return null; }

    public Object[] getChildren() { return new Object[0]; }

    public boolean anyDependencies() { return false; }

    static public StsFieldBean[] getStaticDisplayFields() { return faultDisplayFields; }

    public StsObjectPanel getObjectPanel()
    {
        if (objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsFaultLine").selected(this);
    }
}
