package com.Sts.Actions.Wizards.CopyVelocity;

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

public class StsCopyVelocityPanel extends StsFieldBeanPanel
{
    private StsCopyVelocityWizard wizard;
    private StsCopyVelocity wizardStep;

    transient protected StsGroupBox outputVelocityGroupBox = new StsGroupBox("Output Velocity Model");
    transient protected StsStringFieldBean outputVelocityNameBean;

    transient protected String outputVelocityName = null;

    public StsCopyVelocityPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsCopyVelocityWizard)wizard;
        this.wizardStep = (StsCopyVelocity)wizardStep;
        try
        {
            constructBeans();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        constructPanel();        
    }
 
    private void constructBeans()
    {
        outputVelocityNameBean = new StsStringFieldBean(this, "outputVelocityName", "Output name");
    }

    public void constructPanel()
    {
        outputVelocityGroupBox.add( outputVelocityNameBean);
        add( outputVelocityGroupBox);
    }

    public String getOutputVelocityName() {return outputVelocityName;}
    public void setOutputVelocityName( String outputVelocityName) {this.outputVelocityName = outputVelocityName;}

}
