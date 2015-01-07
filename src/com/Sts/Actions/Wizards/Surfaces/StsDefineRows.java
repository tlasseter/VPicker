
//Title:        S2S: Seismic-to-simulation
//Version:
//Copyright:    Copyright (c) 2002
//Author:       Stuart A. Jackson
//Company:      S2S Systems LLC
//Description:  Web-enabled Integrated Interpretation System

package com.Sts.Actions.Wizards.Surfaces;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

public class StsDefineRows extends StsWizardStep
{
    StsDefineRowsPanel panel;
    StsHeaderPanel header;

    public StsDefineRows(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineRowsPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
 //       panel.setPreferredSize(new Dimension(600, 500));
        header.setTitle("Surface Definition");
        header.setSubtitle("Define the Selected File Rows");
        header.setInfoText(wizardDialog,"(1) Increase the number of header rows.\n" +
                           " **** Header rows will disappear as number is increased. ****\n" +
                           "(2) Once all header rows have disappeared, press the Next>> Button.");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Surfaces");
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
}

