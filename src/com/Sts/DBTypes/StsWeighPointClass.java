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
import com.Sts.UI.Beans.*;

public class StsWeighPointClass extends StsClass implements StsSerializable, StsClassDisplayable
{
    boolean displayWeighPointDialog = false;
    boolean enableTimeJumps = true;

    public StsWeighPointClass()
    {
        userName = "Waypoints";
    }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
        {
            new StsBooleanFieldBean(this, "displayWeighPointDialog", "Enable Notes"),
            new StsBooleanFieldBean(this, "enableTimeJumps", "Enable Time Jumps")
        };
    }

    public void selected(StsWeighPoint weighPoint)
    {
        super.selected(weighPoint);
        setCurrentObject(weighPoint);
    }

    public StsWeighPoint getCurrentWeighPoint()
    {
        return (StsWeighPoint)currentObject;
    }

    public boolean getDisplayWeighPointDialog() { return displayWeighPointDialog; }
    public void setDisplayWeighPointDialog(boolean enable)
    {
        displayWeighPointDialog = enable;
    }

    public boolean getEnableTimeJumps() { return enableTimeJumps; }
    public void setEnableTimeJumps(boolean enable)
    {
    	enableTimeJumps = enable;
    }
    public boolean setCurrentObject(StsObject object)
    {
        boolean changed = super.setCurrentObject(object);
        if(object != null) ((StsWeighPoint)object).treeObjectSelected();
        return changed;
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {

    }

    public boolean setCurrentWeighPointName(String name)
    {
        StsWeighPoint newWeighPoint = (StsWeighPoint)getObjectWithName(name);
        return setCurrentObject(newWeighPoint);
    }

    public void close()
    {
        list.forEach("close");
    }
}
