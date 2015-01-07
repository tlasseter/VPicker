package com.Sts.Actions.Wizards.ClusterAnalysis;

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

public class StsDefineClustering extends StsWizardStep
{
	StsDefineClusteringPanel panel;
    StsHeaderPanel header;

    public StsDefineClustering(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineClusteringPanel((StsClusterAnalysisWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Analyze Clusters");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#ClusterAnalysis");
        header.setInfoText(wizardDialog,"(1) Specify the minimum distance between related events.\n" +
        		" **** enter -1 to exclude distance from analysis. ****\n" +
        		"(2) Specify whether distances are computed in 2D or 3D space.\n" +
        		"(3) Specify the minimum time between events.\n" +
        		" **** enter -1 to exclude time from analysis. ****\n" +
        		"(4) Specify the minimum number of events reqiured in a cluster.\n" +
        		"(5) Optionally, specify the minimum magnitude that corresponds to the provided time and distance.\n" +
        		"(6) Force a linear relationship between time and amplitude (longer time -> greater amplitude).\n" +
        		"(7) Force a linear relationship between distance and amplitude (greater distance -> greater amplitude).\n" +
        		"(8) Press the Run Analysis Button to compute and display results.\n" +
        		"(9) Press the Export View Button to export the events in view to a file.\n" +
        		"(10) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/DefineClustering.html");
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
