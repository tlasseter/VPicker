package com.Sts.Actions.Wizards.FracSim;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.SubVolume.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineSingleSurfaceFracSimVolume extends StsWizardStep
{
    StsDefineSingleSurfaceFracSimVolumePanel panel;
    StsHeaderPanel header;
    String subVolumeName;
    private StsModelSurface surface;
    float topOffset, botOffset;
    byte offsetDomain;
    StsDualSurfaceSubVolume subVolume = null;

    public StsDefineSingleSurfaceFracSimVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineSingleSurfaceFracSimVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-PostStack3d Definition");
        header.setSubtitle("Defining Single Surface Sub-PostStack3d");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FracSim");
        header.setInfoText(wizardDialog,"(1) Specify the name of the sub-volume.\n" +
                           "(2) Select the surface that the sub-volume will be centered on.\n" +
                           "(3) Specify the top and bottom offsets from the selected surface.\n" +
                           "(4) Press the Next> Button to move to the sub-volume setup screen.\n");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        if(subVolume != null) return true;
        if(!((StsSubVolumeWizard)wizard).isCanceled())
        {
            subVolume = StsDualSurfaceSubVolume.constructor(subVolumeName, surface, topOffset, botOffset, offsetDomain, model.getProject().getSeismicVelocityModel());
            model.subVolumeChanged();
        }
        return subVolume != null;
    }

    public void setSurface(StsModelSurface surface)
    {
        this.surface = surface;
    }

    public void setSubVolumeName(String name)
    {
        subVolumeName = name;
    }

    public void setTopOffset(float topOffset)
    {
        this.topOffset = topOffset;
    }

    public void setBotOffset(float botOffset)
    {
        this.botOffset = botOffset;
    }

    public String getSubVolumeName()
    {
        return subVolumeName;
    }

    public StsModelSurface getSurface()
    {
        return surface;
    }

    public float getTopOffset()
    {
        return topOffset;
    }

    public float getBottomOffset()
    {
        return botOffset;
    }

}