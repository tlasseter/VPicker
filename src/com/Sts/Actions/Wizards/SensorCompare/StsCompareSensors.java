package com.Sts.Actions.Wizards.SensorCompare;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsCompareSensors extends StsWizardStep
{
	StsCompareSensorsPanel panel;
    StsHeaderPanel header;

    public StsCompareSensors(StsWizard wizard)
    {
        super(wizard);
        panel = new StsCompareSensorsPanel((StsSensorCompareWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Compare Sensors");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setInfoText(wizardDialog,"(1) Select comparison method (distance, amplitude, xy distance or amplitude).\n" +
        		"(2) Optionally select a threshold attribute used as if condition to comparison method\n" +
        		"(3) Specify the threhold value.\n" +
        		"   **** Ex: Find closest event that has an amplitude of at least ##.# ****\n" +
        		"(4) Press the Run Analysis Button to find related events within the selected sensors\n" +
        		"(5) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorCompare");
    }

    public boolean start()
    {
        panel.initialize();
        wizard.enableFinish();
        return true;
    }

    public boolean end()
    {
        return true;
    }
    
    public void setAnalysisMessage(String msg)
    {
    	panel.setMessage(msg);
    }
}
