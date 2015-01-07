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

public class StsArbitraryLineClass extends StsClass implements StsSerializable
{

    public StsArbitraryLineClass()
    {
        userName = "Arbitrary 3D Line";
    }

    public void selected(StsArbitraryLine aline)
    {
        super.selected(aline);
        setCurrentObject(aline);
    }

    public StsArbitraryLine getCurrentArbitraryLine()
    {
        return (StsArbitraryLine)currentObject;
    }
/*
    public boolean setCurrentObject(StsObject object)
    {
        return super.setCurrentObject(object);
    }
*/
    public boolean setCurrentArbitraryLineName(String name)
    {
        StsArbitraryLine newLine = (StsArbitraryLine)getObjectWithName(name);
        return setCurrentObject(newLine);
    }
}
