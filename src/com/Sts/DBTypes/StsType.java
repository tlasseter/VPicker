package com.Sts.DBTypes;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

import com.Sts.MVC.*;

public class StsType extends StsMainObject
{
    protected StsColor stsColor;

    public StsType()
    {

    }

   public StsType(String name, StsColor color)
    {
        super(false);
        setName(name);
        stsColor = color;
		addToModel();
    }

    public boolean initialize(StsModel model)
    {
        return true;
    }

    public StsColor getStsColor()
    {
        return stsColor;
    }
}
