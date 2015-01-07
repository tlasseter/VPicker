package com.Sts.Types;

import com.Sts.DBTypes.*;

import javax.media.opengl.*;

/**
 * Created by IntelliJ IDEA.
 * User: Tom Lasseter
 * Date: Dec 9, 2008
 * Time: 10:05:55 AM
 * To change this template use File | Settings | File Templates.
 */
abstract public class StsSurfaceTexture
{
    public StsSurface surface;

    public StsSurfaceTexture(StsSurface surface)
    {
        this.surface = surface;
    }

    abstract public byte[] getTextureData();
    abstract public boolean isDisplayable();
    abstract public String getName();
    abstract public int getColorDisplayListNum(GL gl, boolean nullsFilled);
    abstract public StsColorscale getColorscale();
    abstract public float[] getHistogram();
//    abstract public void createColorTLUT(GL gl, boolean nullsFilled);
//    abstract public FloatBuffer getComputeColormapBuffer(boolean nullsFilled);
    public float getDataMin() { return 0.0f; }
    public float getDataMax() { return 0.0f; }
    public String toString() { return getName(); }
    public void selected() { }
}

