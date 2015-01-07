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
import com.Sts.UI.Beans.*;

public class StsBlockClass extends StsClass implements StsSerializable
{
    protected boolean displayBlockColors = false;

    public StsBlockClass()
    {
    }

    public void initializeDisplayFields()
    {
        displayFields = new StsFieldBean[]
        {
           // new StsBooleanFieldBean(this, "displayBlockColors", "Color Blocks")
        };
    }

    public void setDisplayBlockColors(boolean b) { displayBlockColors = b; }
	public boolean getDisplayBlockColors() { return displayBlockColors; }
}
