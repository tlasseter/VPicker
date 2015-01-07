package com.Sts.Actions.Wizards.MicroseismicPreStack;

import com.Sts.Actions.Wizards.PostStack.*;
import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Interfaces.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineTimeLimits extends StsWizardStep
{
	StsDefineTimeLimitsPanel panel;
    StsHeaderPanel header;

    public StsDefineTimeLimits(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineTimeLimitsPanel((StsMicroPreStackWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Load Microseismic PreStack Data");
        header.setSubtitle("Define Time Limits");
        header.setInfoText(wizardDialog,"(1) Specify the starting time for all selected gathers.\n" +
        		"(2) Specify the duration of each gather.\n" +
        		"(3) Press the Finish Button, when done.\n");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MicroseismicPreStack");                                
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
    
    public void analyzeGrid()
    {
         final StsPostStackAnalyzer analyzer = ((StsSeismicWizard)wizard).getAnalyzer(panel.progressPanel, panel.statusTablePanel);
         StsProgressRunnable progressRunnable = new StsProgressRunnable()
         {
              public void cancel()
             {
                 analyzer.cancelProcess();
             }
             public void run()
             {
                 if(!analyzer.analyzeGrid()) return;
             }
         };
         StsToolkit.runRunnable(progressRunnable);
    }  
}
