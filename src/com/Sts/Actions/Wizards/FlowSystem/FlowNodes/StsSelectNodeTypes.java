package com.Sts.Actions.Wizards.FlowSystem.FlowNodes;

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

public class StsSelectNodeTypes extends StsWizardStep
{
    StsSelectNodeTypesPanel panel;
    StsHeaderPanel header;

    public StsSelectNodeTypes(StsWizard wizard)
	{
        super(wizard);
        panel = new StsSelectNodeTypesPanel(wizard, this);
        header = new StsHeaderPanel();
        setPanels(panel, header);
        header.setTitle("Flow Node Definition");
        header.setSubtitle("Select Node Types");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#FlowNodes");
        header.setInfoText(wizardDialog,"(1) Select two types of nodes to relate.\n" +
                           "(3) Press the Next > Button to proceed to the definition screen.\n");
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
