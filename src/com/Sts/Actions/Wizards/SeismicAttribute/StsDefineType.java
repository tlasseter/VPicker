package com.Sts.Actions.Wizards.SeismicAttribute;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.UI.*;

import java.awt.*;
import java.util.*;

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
        panel.setPreferredSize(new Dimension(400, 300));
        header.setTitle("Seismic Attribute Calculation");
        header.setSubtitle("Select Attribute to Calculate");
        header.setInfoText(wizardDialog,"(1) Select the type of volume to produce.\n" +
                                        "(2) Press the Next>> Button.\n");
        header.setLink("http://localhost/marketing/s2ssystems/SeismicCalculation.html");
    }

    public boolean start()
    {
        panel.initialize();
        return true;
    }

    public boolean end()
    {
        if(getAttributes().size() == 0)
        {
            new StsMessage(model.win3d, StsMessage.WARNING, "At least one attribute checkbox must be selected.");
            return false;
        }
        return true;
    }

    public ArrayList getAttributes()
    {
        return panel.getAttributes();
    }
}
