package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineSurfaceSubVolume extends StsWizardStep
{
    StsDefineSurfaceSubVolumePanel panel;
    StsHeaderPanel header;
    String subVolumeName;
    private StsModelSurface surface;
    float topOffset, botOffset;
    byte offsetDomain;
    StsDualSurfaceSubVolume subVolume = null;

    public StsDefineSurfaceSubVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineSurfaceSubVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-Volume Definition");
        header.setSubtitle("Defining Single Surface Sub-Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SubVolume");
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
        StsSubVolumeWizard svWizard = (StsSubVolumeWizard)wizard;
        if(!svWizard.isCanceled())
        {
            subVolume = StsDualSurfaceSubVolume.constructor(svWizard.getSubVolumeName(), surface, topOffset, botOffset, offsetDomain, model.getProject().getSeismicVelocityModel());
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
    
    public String getOffsetDomainString()
    {
        return StsParameters.TD_STRINGS[offsetDomain];
    }

    public void setOffsetDomainString(String domainString)
    {
    	offsetDomain = StsParameters.getZDomainFromString(domainString);
    }
}
