package com.Sts.Actions.Wizards.HandVelocity;

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
public class StsHandVelocitySelect extends StsWizardStep
{
    StsHandVelocitySelectPanel panel;
    StsHeaderPanel header;

    public StsHandVelocitySelect(StsWizard wizard)
    {
        super(wizard);
        panel = new StsHandVelocitySelectPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 400));
        header.setTitle("Hand Velocity (HANDVEL) Selection");
        header.setSubtitle("<html> Select one or more hand picked velocity. <br> Available files have .handvel or .dat suffix");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#HandVel");                
        header.addInfoLine(wizardDialog, "(1) Navigate to the desired file directory using the Dir button.");
        header.addInfoLine(wizardDialog, "     ***** All available objects in the selected directory will be placed in the left list.");
        header.addInfoLine(wizardDialog, "     ***** The default location is the project directory but any directory can be selected");
        header.addInfoLine(wizardDialog, "(2) Select the hand velocity (HANDVEL) files from the left list and place them in the right list");
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
