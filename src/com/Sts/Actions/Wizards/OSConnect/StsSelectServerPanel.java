package com.Sts.Actions.Wizards.OSConnect;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;


/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectServerPanel extends StsJPanel
{
    private StsOSConnectWizard wizard;
    private StsSelectServer wizardStep;

    public StsSelectServerPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsOSConnectWizard)wizard;
        this.wizardStep = (StsSelectServer)wizardStep;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
    	;
    }

    void jbInit() throws Exception
    {
    	;
    }
}
