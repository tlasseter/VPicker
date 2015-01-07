package com.Sts.Actions.Wizards.AnalogueCube;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author unascribed
 * @version 1.1
 */

public class StsDefineType extends StsWizardStep
{
    StsDefineTypePanel panel;
    StsHeaderPanel header;

    public StsDefineType(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineTypePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(450, 400));
        header.setTitle("Aspect Energy Analogue Cube Analysis");
        header.setSubtitle("Define Analysis Type");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AnalogueCube");
        header.setLogo("AspectLogo.gif");
        header.setInfoText(wizardDialog,"(1) Select the type of volume to produce.\n" +
                           "      Real - the data from the source and target volumes will be analyzed directly\n" +
                           "      Complex - the data from the source and target volumes will be transformed \n" +
                           "                using a Hilbert transform to produce an imaginary volume and both\n" +
                           "                the real and imaginary volumes will be analyzed.\n" +
                           "(2) Press the Next>> Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/aspect/AspectAnalogueCube.html");
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

    public byte getType()
    {
        return panel.getType();
    }
/*
    public boolean persistVolumes()
    {
        return panel.persistVolumes();
    }
*/
}
