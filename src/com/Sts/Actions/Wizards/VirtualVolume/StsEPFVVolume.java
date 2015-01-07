package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
import com.Sts.Interfaces.*;
import com.Sts.UI.Progress.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsEPFVVolume extends StsWizardStep
{
    StsEPFVVolumePanel panel;
    StsHeaderPanel header;
    StsEPFVirtualVolume virtualVolume = null;
    String volumeName = null;
    StsProgressPanel ppanel = null;

    public StsEPFVVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsEPFVVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Defining Filtered Virtual Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Select a volume to apply the smoothing filter.\n" +
                           "(2) Select the type of filter to apply.\n" +
                           "(3) Define the filter parameters.\n" +
                           "(4) Press the Finish>> Button to create the volume.\n");
    }

    public boolean start()
	{
        panel.initialize();
        ppanel = null;
        return true;
    }

    public boolean end()
    {
    	return true;
    }

    public boolean buildVolume(StsProgressPanel ppanel)
    {
    	this.ppanel = ppanel;
		if(virtualVolume != null) return true;
        panel.applyFilter();
        return true;
    }

    public void applyFilter(StsSeismicVolume volume, StsVolumeFilterFace filter, String name)
    {
    	panel.initializeVolumeName();
//        name = ((StsVirtualVolumeWizard)wizard).getVolumeName();
        if(ppanel != null) ppanel.appendLine("Constructing filter virtual volume:" + name);
        if(name.equals("VVName") || (name == null) || (name.length() <= 0))
        	name = panel.getVolumeName();
        virtualVolume = new StsEPFVirtualVolume(volume, filter, name);
        if (virtualVolume == null)
        {
        	if(ppanel != null) ppanel.appendLine("Failed to construct blended virtual volume.\n");
            return;
        }
 //       if(((StsVirtualVolumeWizard)wizard).isFloatVolume())
 //       {
 //       	System.out.println("Not available yet.");
 //       }
    	virtualVolume.computeHistogram();
        if(ppanel != null) ppanel.appendLine("Successfully created filter virtual volume:" + name);
    }

    public void setVolumeName(String name) { volumeName = name; }
}
