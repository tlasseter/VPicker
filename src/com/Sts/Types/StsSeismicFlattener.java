package com.Sts.Types;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Toolbars.*;
import com.Sts.Utilities.*;

public class StsSeismicFlattener
{
    protected StsSeismicVolume seismicVolume = null;
    protected StsSurface surface = null;
    protected StsModel model = null;

    transient public StsToolbar toolbar = null;

    public StsSeismicFlattener()
    {
    }

    public StsSeismicFlattener(StsModel model)
    {
        this.model = model;
    }

    static public StsSeismicFlattener constructor(StsModel model)
    {
        try
        {
            StsSeismicFlattener flattener = new StsSeismicFlattener(model);
            return flattener;
        }
        catch(Exception e)
        {
            StsException.outputException("StsSeismicFlattener.constructTraceAnalyzer() failed.",
                e, StsException.WARNING);
            return null;
        }
    }

    public void initialize(StsSeismicVolume seismicVolumeIn, StsSurface surfaceIn) throws StsException
    {
        this.seismicVolume = seismicVolumeIn;
        this.surface = surfaceIn;
    }

    public StsSeismicVolume getSeismicVolume()
    {
        return seismicVolume;
    }

    public StsSurface getSurface()
    {
        return surface;
    }
    public StsSeismicVolume setSeismicVolume(StsSeismicVolume seismicVolumeIn)
    {
        seismicVolume = seismicVolumeIn;
        return seismicVolume;
    }

    public StsSurface setSurface(StsSurface surfaceIn)
    {
        surface = surfaceIn;
        return surface;
    }

    public StsToolbar getToolbar() { return toolbar; }
    public void setToolbar(StsToolbar tb)
    {
        toolbar = tb;
        return;
    }
}
