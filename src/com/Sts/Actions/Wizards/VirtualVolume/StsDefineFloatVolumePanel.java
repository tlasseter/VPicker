package com.Sts.Actions.Wizards.VirtualVolume;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;
import com.Sts.Utilities.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineFloatVolumePanel extends StsJPanel
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;
    private StsFloatFieldBean userNullBean = new StsFloatFieldBean();
    private float userNull = StsParameters.nullValue;
    
    public StsDefineFloatVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        try
        {
        	constructBeans();
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void constructBeans()
    {
    	userNullBean.initialize(this, "userNull", true, "Null Value:");
    	return;
    }

    void jbInit() throws Exception
    {
    	gbc.fill = gbc.HORIZONTAL;
    	addEndRow(userNullBean);
    }

    public void setUserNull(float nullVal)
    {
    	userNull = nullVal;
    }
    public float getUserNull()
    {
    	return userNull;
    }
}
