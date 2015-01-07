package com.Sts.Actions.Wizards.PreStackExport;

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

public class StsExportType extends StsWizardStep
{
	StsExportTypePanel panel;
    StsHeaderPanel header;

    public StsExportType(StsWizard wizard)
	{
        super(wizard);
        panel = new StsExportTypePanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 200));
        header.setTitle("Export Pre-Stack Data");
        header.setSubtitle("Select Type of Data to Export");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/PreStackDataExport.html");
        header.setInfoText(wizardDialog,"(1) Select the type of pre-stack data you wish to export.\n" +
                           " (2) Select the dimension of data.\n" +
                           "     *** Disabled if project only contains 2D or 3D but not both.\n" +
                           "(3) Press the Next > Button to proceed to the data selection step.\n");
    }

    public boolean start()
    {
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public byte getExportType()
    {
        return panel.getExportType();
    }
    public byte getExportDimension()
    {
        return panel.getExportDimension();
    }

}
