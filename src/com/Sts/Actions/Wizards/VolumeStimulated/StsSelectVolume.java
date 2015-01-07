package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.VirtualVolume.StsMathVVolumePanel;
import com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.DBTypes.StsMathVirtualVolume;
import com.Sts.DBTypes.StsSeismicVolume;
import com.Sts.DBTypes.StsVirtualVolume;
import com.Sts.SeismicAttributes.StsVirtualVolumeConstructor;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.UI.StsMessage;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsSelectVolume extends StsWizardStep
{
    StsSelectVolumePanel panel;
    StsHeaderPanel header;
	StsVirtualVolume virtualVolume = null;
	StsSeismicVolume[] volumes = null;

    public StsSelectVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsSelectVolumePanel((StsVolumeStimulatedWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Volume Stimulated");
        header.setSubtitle("Compute the Stimulated Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#StimulatedVolume");
        header.setInfoText(wizardDialog,"(1) Select the volume to compute on.\n" +
                           "(2) Press the Run button to compute the volume stimulated.\n" +
                           "(3) Press the Finish> Button to dismiss the screen.\n");
    }

    public boolean start()
	{
        enableFinish();
        wizard.disableNext();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
    	return true;
    }

    public StsSeismicVolume getVolume()
    {
        return panel.getVolume();
    }


    public void setAnalysisMessage(String msg)
    {
    	panel.setMessage(msg);
    }    
}