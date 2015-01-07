package com.Sts.Actions.Wizards.FracSim;

import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineFracSimVolumePanel extends StsJPanel
{
    private StsFracSimWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;
    StsComboBoxFieldBean typesBean;
    StsStringFieldBean nameBean;

    public StsDefineFracSimVolumePanel(StsFracSimWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
            constructPanel();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
    }

    void constructPanel() throws Exception
    {
        typesBean = new StsComboBoxFieldBean(wizard, "type", "Model type", StsFracSimWizard.fracSimTypes);
        nameBean = new StsStringFieldBean(wizard, "name", true, "Name");
        add(typesBean);
        add(nameBean);
    }
}