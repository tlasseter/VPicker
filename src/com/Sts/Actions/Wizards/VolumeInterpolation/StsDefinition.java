package com.Sts.Actions.Wizards.VolumeInterpolation;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
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

public class StsDefinition extends StsWizardStep
{
    StsDefinitionPanel panel;
    StsHeaderPanel header;

    public StsDefinition(StsWizard wizard)
	{
        super(wizard);
        panel = new StsDefinitionPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Volume Interpolation");
        header.setSubtitle("Volume Selection & Interpolation Type");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VolumeInterpolation");
        header.setInfoText(wizardDialog,"(1) \n");
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
        String name = null;
        StsSeismicBoundingBox volume = getSelectedVolume();
        if(volume == null)
            return false;

        name = ((StsVolumeInterpolationWizard)wizard).getVolumeName();
        // Create new volume on disk
        
        
        // Load new volume into project
        
        
        return true;
    }

    public StsSeismicBoundingBox getSelectedVolume()
    {
    	return panel.getSelectedVolume();
    }

}
