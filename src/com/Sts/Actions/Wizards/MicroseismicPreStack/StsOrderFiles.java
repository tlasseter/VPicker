
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2001
//Author:       Tom Lasseter
//Company:      4D Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.MicroseismicPreStack;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.Types.*;
import com.Sts.UI.*;

public class StsOrderFiles extends StsWizardStep
{
    StsOrderedListPanel panel;
    StsHeaderPanel header;
    StsSeismicBoundingBox[] boxes = null;
    int[] selected = null;
    StsMicroPreStackWizard wizard;

    public StsOrderFiles(StsWizard wizard)
    {
        super(wizard, new StsOrderedListPanel(), null, new StsHeaderPanel());
        this.wizard = (StsMicroPreStackWizard)wizard;
        panel = (StsOrderedListPanel) getContainer();
        header = (StsHeaderPanel) getHdrContainer();
        header.setTitle("Microseismic PreStack Load");
        header.setSubtitle("Time Order Files");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MicroseismicPreStack");
        header.setInfoText(wizardDialog,"(1) Select files and use up and down arrows to organize into chronological order.\n" +
                                  "   ***** Oldest at Top ***** \n" +
                                  "(2) Press the Next>> Button.");
    }

    public boolean start()
    {
        if(model == null) return true;

    	if(boxes == null)
        {
			boxes = (StsSeismicBoundingBox[])wizard.getVolumes();
            panel.setItems(boxes);
        }
        return boxes != null;
    }
    
    public Object[] getOrderedVolumes()
    {
    	return panel.getItems();
    }

    public boolean end()
    {
        return true;
    }
}
