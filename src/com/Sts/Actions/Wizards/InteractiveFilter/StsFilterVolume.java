package com.Sts.Actions.Wizards.InteractiveFilter;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsFilterVolume extends StsWizardStep
{
    StsFilterVolumePanel panel;
    StsHeaderPanel header;
	StsFilterVirtualVolume virtualVolume = null;

    public StsFilterVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsFilterVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Filter Definition");
        header.setSubtitle("Defining Filter");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FilterVolume");                
        header.setInfoText(wizardDialog,"(1) Select the type of filter to apply.\n" +
                           "(2) Define the filter parameters.\n" +
                           "(3) Press the Finish>> Button to apply filter to current slices.\n");
    }

    public boolean start()
	{
        panel.initialize();
        wizard.enableFinish();
        disableNext();
        return true;
    }

    public boolean end()
    {
        if(virtualVolume != null)
           return true;
        StsSeismicBoundingBox[] volumes = new StsSeismicBoundingBox[] { StsPreStackLineSetClass.currentProjectPreStackLineSet.velocityModel };
    /*
        if(panel.getFilterType() == StsFilterVirtualVolume.CONVOLUTION)
            virtualVolume = new StsFilterVirtualVolume(volumes, "filterVolume", StsFilterVirtualVolume.CONVOLUTION,
                                                          panel.getKernel());
        else
            virtualVolume = new StsFilterVirtualVolume(volumes, "filterVolume", StsFilterVirtualVolume.RANK,
                                                          panel.getSubType(), panel.getFilterRadius());
    */
        return true;
    }

}
