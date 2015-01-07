package com.Sts.DBTypes;

import com.Sts.DB.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */


public class StsAttributeSurfaceClass extends StsSurfaceClass implements StsSerializable
{
    public StsAttributeSurfaceClass()
    {
        userName = "Attribute Surfaces";
    }

    public boolean textureChanged(StsObject object)
    {
        return false;
    }
}