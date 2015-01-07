package com.Sts.Actions.Wizards.SurfaceCurvature;

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

public class StsProcessSurface extends StsWizardStep
{
	StsProcessSurfacePanel panel;
    StsHeaderPanel header;

    public StsProcessSurface(StsWizard wizard)
    {
        super(wizard);
        panel = new StsProcessSurfacePanel((StsSurfaceCurvatureWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Calculate Curvature");
        header.setSubtitle("Define Criteria & Analyze Results");
        header.setInfoText(wizardDialog,"(1) Select Curvature attribute to calculate.\n" +
        		"(2) Select analysis window size.\n" +
        		"(3) Press the Run Analysis Button to compute and display results.\n" +
        		"(4) Press the Save To Model Button to create a persistent texture surface.\n" +
        		//"(5) Press the Export View Button to export the events in view to a file.\n" +
        		"(5) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#SurfaceCurvature");
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
    
    public void updateColorscale()
    {
    	panel.updateColorscale();
    }
}
