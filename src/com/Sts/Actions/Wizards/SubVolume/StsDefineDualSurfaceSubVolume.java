package com.Sts.Actions.Wizards.SubVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.*;
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

public class StsDefineDualSurfaceSubVolume extends StsWizardStep
{
    StsDefineDualSurfaceSubVolumePanel panel;
    StsHeaderPanel header;
    String subVolumeName;
    StsModelSurface[] surfaces = new StsModelSurface[2];
    float topOffset, botOffset;
    byte offsetDomain;
    StsDualSurfaceSubVolume dualSurfaceSubVolume = null;

    public StsDefineDualSurfaceSubVolume(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineDualSurfaceSubVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 500));
        header.setTitle("Sub-Volume Definition");
        header.setSubtitle("Defining Dual Surface Sub-Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SubVolume");
        header.setInfoText(wizardDialog,"(1) Specify the name of the sub-volume.\n" +
                           "(2) Select the top surface to constrain the sub-volume.\n" +
                           "(3) Specify the offset from the top surface.\n" +
                           "    **** positive is below, negative is above. ****\n" +
                           "(4) Select the bottom surface to constrain the sub-volume.\n" +
                           "(5) Specify the offset from the bottom surface.\n" +
                           "(6) Press the Next> Button to move to the sub-volume setup screen.\n");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        StsSubVolumeWizard svWizard = (StsSubVolumeWizard)wizard;    	
        if(!svWizard.isCanceled())
        {
            if(dualSurfaceSubVolume == null)
                dualSurfaceSubVolume = StsDualSurfaceSubVolume.constructor(svWizard.getSubVolumeName(), surfaces, topOffset, botOffset, offsetDomain, model.getProject().getSeismicVelocityModel());
            model.subVolumeChanged();
        }
        else
        {
            if(dualSurfaceSubVolume != null)
                dualSurfaceSubVolume.delete();
            dualSurfaceSubVolume = null;
        }
        return dualSurfaceSubVolume != null;
    }

    public void setTopSurface(StsModelSurface topSurface)
    {
        if(surfaces[0] == topSurface) return;
        surfaces[0] = topSurface;
        panel.checkSetTopSurfaceSelected();
        checkSurfaces();
    }

    public void setBotSurface(StsModelSurface botSurface)
    {
        if(surfaces[1] == botSurface) return;
        surfaces[1] = botSurface;
        panel.checkSetBotSurfaceSelected();
        checkSurfaces();
    }

    private void checkSurfaces()
    {
        if((surfaces[0] == null) || (surfaces[1] == null)) return;
        if(surfaces[0] != surfaces[1]) return;
        new StsMessage(wizard.frame, StsMessage.WARNING, "The top surface must be different than the bottom surface.");
    }

    public StsModelSurface getTopSurface() { return surfaces[0]; }
    public StsModelSurface getBotSurface() { return surfaces[1]; }

    public void setTopOffset(float topOffset) { this.topOffset = topOffset; }
    public void setBotOffset(float botOffset) { this.botOffset = botOffset; }
    public float getTopOffset() { return topOffset; }
    public float getBottomOffset() { return botOffset; }

    public String getOffsetDomainString()
    {
        return StsParameters.TD_STRINGS[offsetDomain];
    }

    public void setOffsetDomainString(String domainString)
    {
    	offsetDomain = StsParameters.getZDomainFromString(domainString);
    }

    public void setSubVolumeName(String name) { this.subVolumeName = name; }
	public String getSubVolumeName() { return subVolumeName; }
}
