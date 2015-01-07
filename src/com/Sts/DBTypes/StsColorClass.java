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

public class StsColorClass extends StsClass implements StsSerializable
{
    protected StsColor defaultLabelColor = new StsColor(StsColor.WHITE);
    protected StsColor defaultWellColor = new StsColor(StsColor.RED);
    protected StsColor defaultSurfaceColor = new StsColor(StsColor.BLUE);

    public StsColorClass()
    {
        userName = name;
    }

    public StsColor getDefaultLabelColor() { return defaultLabelColor; }
    public StsColor getDefaultWellColor() { return defaultWellColor; }
    public StsColor getDefaultSurfaceColor() { return defaultSurfaceColor; }

    public void setDefaultLabelColor(StsColor color) { defaultLabelColor = color; }
    public void setDefaultWellColor(StsColor color) { defaultWellColor = color; }
    public void setDefaultSurfaceColor(StsColor color) { defaultSurfaceColor = color; }
}
