package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;
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

public class StsCrossplotVVolume extends StsWizardStep
{
    StsCrossplotVVolumePanel panel;
    StsHeaderPanel header;
	StsCrossplotVirtualVolume virtualVolume = null;

    public StsCrossplotVVolume(StsWizard wizard)
	{
        super(wizard);
        panel = new StsCrossplotVVolumePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Virtual Volume Definition");
        header.setSubtitle("Defining Crossplot Virtual Volume");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#VirtualVolume");
        header.setInfoText(wizardDialog,"(1) Select the volume and crossplot to operate on.\n" +
                           "(2) Select whether the merge inclusively or exclusively.\n" +
                           "(3) Press the Finish> Button to create the volume.\n" +
                           "     ***** A virtual volume combobox will be added to the toolbar panel to control the viewing *****\n");
    }

    public boolean start()
	{
        panel.initialize();
        return true;
    }

    public boolean end()
    {
    	return true;
    }
    
    public boolean buildVolume(StsProgressPanel ppanel)
    {
		if(virtualVolume != null) return true;

        String name = null;
        StsSeismicVolume volume = (StsSeismicVolume) panel.getSelectedXplotSeismicVolume();
        StsCrossplot xplot = panel.getSelectedCrossplot();
        boolean inclusive = panel.isInclusive();
        name = ((StsVirtualVolumeWizard)wizard).getVolumeName();
        if(name.equals("VVName") || (name == null) || (name.length() <= 0))
        {
            if(inclusive)
                name = "Incl_" + volume.getName() + "_" + xplot.getName();
            else
                name = "Excl_" + volume.getName() + "_" + xplot.getName();
        }
        ppanel.appendLine("Constructing crossplot virtual volume:" + name);        
        virtualVolume = new StsCrossplotVirtualVolume(volume, xplot, name, inclusive);
        if (virtualVolume == null)
        {
        	ppanel.appendLine("Failed to construct crossplot virtual volume.\n");
            return false;
        }
        if(((StsVirtualVolumeWizard)wizard).isFloatVolume())
        {
        	System.out.println("Not available yet.");
        } 
    	virtualVolume.computeHistogram();        
        ppanel.appendLine("Successfully created crossplot virtual volume:" + name);        
        return true;
    }
}
