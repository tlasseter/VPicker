package com.Sts.Actions.Wizards.Fracture;

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

public class StsDefineFractureSet extends StsWizardStep
{
    StsDefineFractureSetPanel panel;
    StsHeaderPanel header;

    public StsDefineFractureSet(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineFractureSetPanel((StsFractureWizard)wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Fracture Definition");
        header.setSubtitle("Defining Fracture Properties");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#Fracture");
        header.setInfoText(wizardDialog,"(1) Specify the name of the fracture set\n" +
        		"(2) Select the desired color. All fractures will be semi-transparent\n" +
        		"(3) Select an existing fracture set to use as a guide for the new set.\n" +
        		" *****  This selection will use the same Z planes from selected set. \n" +
        		" *****  to define the new set. Additional edges can be added if desired.\n" +
        		"(4) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/DefineFracture.html");
    }

    public boolean start()
    {
        wizard.enableNext();
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        return true;
    }
}
