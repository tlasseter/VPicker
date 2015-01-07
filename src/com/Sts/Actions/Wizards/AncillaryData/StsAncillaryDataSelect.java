package com.Sts.Actions.Wizards.AncillaryData;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.WizardHeaders.*;

import java.awt.*;


/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2002</p>
 *
 * <p>Company: </p>
 *
 * @author not attributable
 * @version 1.0
 */
public class StsAncillaryDataSelect extends StsWizardStep
{
    StsAncillaryDataSelectPanel panel;
    StsHeaderPanel header;

    public StsAncillaryDataSelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsAncillaryDataSelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 400));
        header.setTitle("Ancillary Data File Selection");
        header.setSubtitle("Selecting one or more ancillary data files");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#AncillaryData");
        header.addInfoLine(wizardDialog, "(1) Navigate to the desired directory using the Dir button.");
        header.addInfoLine(wizardDialog, "     ***** All available files in the selected directory will be placed in the left list.");
        header.addInfoLine(wizardDialog, "     ***** The default location is the project directory but any directory can be selected");
        header.addInfoLine(wizardDialog, "(2) Select the desired files from the left list and place them in the right list");
        header.addInfoLine(wizardDialog, "    using the provided controls between the lists.");
        header.addInfoLine(wizardDialog, "(3) Once all object selections are complete, press the Next>> Button.");
    }

    public String getSelectedDirectory() {  return panel.getCurrentDirectory(); }

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
