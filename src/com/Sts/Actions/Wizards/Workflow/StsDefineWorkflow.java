package com.Sts.Actions.Wizards.Workflow;

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

public class StsDefineWorkflow extends StsWizardStep
{
    StsDefineWorkflowPanel panel;
    StsHeaderPanel header;

    public StsDefineWorkflow(StsWizard wizard)
    {
        super(wizard);
        panel = new StsDefineWorkflowPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        panel.setPreferredSize(new Dimension(500, 300));
        header.setTitle("Workflow Definition");
        header.setSubtitle("Selecting Workflow");
        header.setInfoText(wizardDialog,"(1) Select the desired workflow.\n" +
                           "     **** A description of each workflow is provided upon selection.\n" +
                           "(2) Press the Finish Button.\n");
        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/Workflow.html");
    }

    public boolean start()
    {
        wizard.dialog.setTitle("Selecting Workflow");
        return true;
    }

    public boolean end()
    {
        return true;
    }

    public byte getType() { return panel.getType(); }
}
