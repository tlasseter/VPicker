
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.DBTypes.*;

import java.awt.*;

public class StsProcessCube extends StsWizardStep implements Runnable
{
    public StsProcessCubePanel panel;
    public StsHeaderPanel header;

    public StsProcessCube(StsWizard wizard)
    {
        super(wizard);
        header = new StsHeaderPanel();
        panel = new StsProcessCubePanel((StsAnalogueCubeWizard)wizard, this);
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(450, 400));
        header.setTitle("Aspect Energy Analogue Cube Analysis");
        header.setSubtitle("Process Analogue Cube");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AnalogueCube");
        header.setLogo("AspectLogo.gif");
        header.setInfoText(wizardDialog,"(1) Set the correlation optimization method.\n" +
                           "  **** Running Average continues calculation if average exceeds threshold.\n" +
                           "  **** Spiral Average continues calculation if average for spiral exceeds threshold.\n" +
                           "  **** Spiral Max continues calculation if one trace in spiral exceeds threshold. \n" +
                           "(2) Set the correlation threshold to determine when to abandon a subvolume correlation\n" +
                           "  **** The analysis is done radially out from the center of the subvolume. \n" +
                           "       If a spiral is completed with no one correlation meeting the correlation \n" +
                           "       threshold for the selected method the remaining correlations are abandoned \n" +
                           "       and an estimate correlation is computed for the sample. \n" +
                           "(3) Press the Preview Button to enter preview mode.\n" +
                           "  **** All slices in view will be processed and each time a slice is moved\n" +
                           "       it will be automatically processed.\n" +
                           "(4) Press the Execute Button to run the processing.\n" +
                           "  **** This process could take awhile, please be patient.\n" +
                           "(5) Press the Cancel Button to cancel the operation and not process the cube.\n" +
                           "(6) Pressing the Finish Button before the processing is complete will \n" +
                           "    cause processing to be run in background\n");
        header.setLink("http://www.s2ssystems.com/marketing/aspect/AspectAnalogueCube.html");
    }

    public boolean start()
    {
		StsSubVolume subVolume = ((StsAnalogueCubeWizard)wizard).getSubVolume();
		if(subVolume != null) subVolume.setIsApplied(true);
//        panel.initializePanel();
        panel.dryRunBtn.setSelected(false);
        return true;
    }

    public boolean end()
    {
        panel.dryRunBtn.setSelected(false);
        ((StsAnalogueCubeWizard)wizard).endDataPlaneCalculation();

        // If the foreground process was executed ask the user if they want to run again in background
        if(panel.wasRunForeground())
        {
//            if(!StsMessage.questionValue(wizard.frame,"Already run in foreground, want to run again in background?"))
 //               return true;
        }

        try
        {
  //          Runtime rt = Runtime.getRuntime();
  //          Process prcs = rt.exec(executionCmd);
        }
        catch(Exception e)
        {
            System.err.println("StsCoherenceVolumeProcess:end() Error executing command" + e);
        }

        return true;
    }

    public float getCorrelationThreshold()
    {
        return panel.getCorrelationThreshold();
    }

    public byte getOptimizationMethod()
    {
        return panel.getOptimizationMethod();
    }

}
