package com.Sts.Actions.Wizards.SensorCorrelation;

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

public class StsCorrelateSensors extends StsWizardStep
{
	StsCorrelateSensorsPanel panel;
    StsHeaderPanel header;

    public StsCorrelateSensors(StsWizard wizard)
    {
        super(wizard);
        panel = new StsCorrelateSensorsPanel((StsSensorCorrelationWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Correlate Sensors");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setInfoText(wizardDialog,"(1) Select comparison method (average, min, max, sum).\n" +
        		"(2) Define the range of data to analyze.\n" +
        		"(3) Press the Run Analysis Button to find related events within the selected sensors\n" +
        		"(4) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SensorCorrelation");
    }

    public boolean start()
    {
        wizard.enableFinish();
        panel.initialize();
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