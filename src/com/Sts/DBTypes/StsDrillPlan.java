package com.Sts.DBTypes;

import com.Sts.DB.DBCommand.*;
import com.Sts.Interfaces.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.ObjectPanel.*;

import java.awt.*;

public class StsDrillPlan extends StsMainObject implements StsTreeObjectI
{
    public final static int PLATFORM = 0;
    public final static int ARBITRARY = 1;
    public final static int SATURATION = 2;

    static final public String[] dpTypes = new String[]
		{"Platform Slots", "Arbitrary Vertical Paths", "Vertical Saturation Set"};

    static StsObjectPanel objectPanel = null;

    protected StsObjectRefList plannedWells = null;
    protected int planType = PLATFORM;
    protected StsColor defaultColor = null;
    protected double xOrigin, yOrigin;

    // display fields
    static public final StsFieldBean[] planDisplayFields = null;

    public StsDrillPlan()
    {
    }

    private void initialize()
    {
  	    plannedWells = StsObjectRefList.constructor(2, 2, "plannedWells", this);
        StsProject project = currentModel.getProject();
        xOrigin = project.getXOrigin();
        yOrigin = project.getYOrigin();
    }
/*
    static public StsDrillPlan constructor(String name)
    {
        try
        {
            StsDrillPlan drillPlan = new StsDrillPlan();
            drillPlan.setName(name);
            return drillPlan;
        }
        catch(Exception e)
        {
            StsException.outputException("StsDrillPlan.constructor(win3d) failed.", e, StsException.WARNING);
            return null;
        }
    }

    static public StsDrillPlan constructor()
    {
        String name = new String("DrillPlan");
        StsDrillPlan drillPlan = StsDrillPlan.constructor(name);
        return drillPlan;
    }

    static public StsDrillPlan constructor(int type)
    {
        String name = new String("DrillPlan");
        StsDrillPlan drillPlan = StsDrillPlan.constructor(name);
        drillPlan.planType = type;
        return drillPlan;
    }

    static public StsDrillPlan constructor(String name, int type)
    {
        StsDrillPlan drillPlan = StsDrillPlan.constructor(name);
        drillPlan.planType = type;
        return drillPlan;
    }
*/
    /**
     * constructor
     */
    static public StsDrillPlan constructor(int type, String name, StsColor color)
    {
        StsDrillPlan drillPlan = new StsDrillPlan();
        drillPlan.initialize();
        drillPlan.setName(name);
        drillPlan.planType = type;
        drillPlan.defaultColor = color;
        return drillPlan;
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public boolean delete()
    {
        if(!super.delete()) return false;
        return true;
    }

    public double getXOrigin() { return xOrigin; }
    public void setXOrigin(double xOrigin) { this.xOrigin = xOrigin; }
    public double getYOrigin() { return yOrigin; }
    public void setYOrigin(double yOrigin) { this.yOrigin = yOrigin; }

    public StsObjectRefList getPaths()
    {
        return plannedWells;
    }

    public void addWell(StsWell well)
    {
        plannedWells.add(well);
    }

    public void setStsColor(Color color)
    {
        defaultColor = new StsColor(color);
        /*
        if(defaultColor == null)
            defaultColor = new StsColor(color);
        else
            defaultColor.setBeachballColors(color);
        */
        currentModel.addTransactionCmd("default well color change", new StsChangeCmd(this, defaultColor, "stsColor", false));
    }
    public StsColor getStsColor() { return defaultColor; }
    public StsColor getDefaultWellColor() { return defaultColor; }

    public StsFieldBean[] getDisplayFields() { return planDisplayFields; }
    public StsFieldBean[] getPropertyFields() { return null; }
    public Object[] getChildren() { return new Object[0]; }
    public boolean anyDependencies() { return false; }

    public void display(StsGLPanel3d glPanel)
    {
        if (glPanel == null)
        {
            return;
        }
        super.display(glPanel);
    }

    public StsObjectPanel getObjectPanel()
    {
        if(objectPanel == null) objectPanel = StsObjectPanel.constructor(this, true);
        return objectPanel;
    }

    public void treeObjectSelected()
    {
        currentModel.getCreateStsClass("com.Sts.DBTypes.StsDrillPlan").selected(this);
    }
}