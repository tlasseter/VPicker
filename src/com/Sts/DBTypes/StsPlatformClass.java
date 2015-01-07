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

import javax.media.opengl.*;
import java.util.*;

public class StsPlatformClass extends StsClass implements StsSerializable, StsClassDisplayable
{
    /** Index into the spectrum used for platform color.
     * Initialized to the current number of platforms minus one.
     * Incremented as platforms are added.
     */
    transient int colorIndex = 0;

    /** width and height of platform display in pixels */
     static final int platformDisplaySize = 40;

    public StsPlatformClass()
    {
        userName = "Drilling Platform";
    }

    public boolean dbInitialize()
    {
        colorIndex = getSize() - 1;
        return true;
    }

    public void deleteWellFromPlatform(String wellname)
    {
        Iterator platformIterator = currentModel.getObjectIterator(StsPlatform.class);
        while(platformIterator.hasNext())
        {
            StsPlatform platform = (StsPlatform)platformIterator.next();
            if(platform.isWellOnPlatform(wellname))
                platform.deleteWellfromPlatform(wellname);
        }
    }

    public StsPlatform getWellPlatform(String wellname)
    {
        Iterator platformIterator = currentModel.getObjectIterator(StsPlatform.class);
        while(platformIterator.hasNext())
        {
            StsPlatform platform = (StsPlatform)platformIterator.next();
            if(platform.isWellOnPlatform(wellname))
                return platform;
        }
        return null;
    }

    public StsColor getNextColor(StsSpectrum spectrum)
    {
        int nColors = spectrum.getNColors();
        colorIndex++;
        if(colorIndex >= nColors) colorIndex = 0;
        return new StsColor(spectrum.getColor(colorIndex));
    }

    public void displayClass(StsGLPanel3d glPanel3d)
    {
        if(!isVisible) return;

        int nPlatforms = getSize();
        if(nPlatforms == 0) return;
        GL gl = glPanel3d.getGL();
        gl.glDisable(GL.GL_LIGHTING);
        float distance = ((StsView3d)glPanel3d.getView()).distance;
        float projectSize = currentModel.getProject().getSize();
        float size = 0.01f*Math.min(distance, projectSize);
        for(int n = 0; n < nPlatforms; n++)
        {
            StsPlatform platform = (StsPlatform)getElement(n);
            platform.display(glPanel3d, size);
        }
        gl.glEnable(GL.GL_LIGHTING);
    }
}