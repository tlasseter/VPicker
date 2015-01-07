package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.UI.Beans.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public abstract class StsDefineVelocityAnalysisPanel extends StsFieldBeanPanel
{
	private StsVelocityAnalysisWizard wizard;
	StsPreStackLineSet volume;
	StsPreStackVelocityModel velocityModel;

	boolean initialized = false;

    protected StsGroupBox defineGroupBox = new StsGroupBox("Analysis Increments");
    protected StsFloatFieldBean analysisRowNumStartBean;
	protected StsFloatFieldBean analysisColNumStartBean;
	protected StsFloatFieldBean analysisRowNumIncBean;
	protected StsFloatFieldBean analysisColNumIncBean;
    protected StsIntFieldBean analysisColIncBean;
	protected StsIntFieldBean traceThresholdBean;

	public StsDefineVelocityAnalysisPanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = (StsVelocityAnalysisWizard)wizard;
		try
		{
			//volume = StsPreStackLineSetClass.currentProjectPreStackLineSet;//this.wizard.getPreStackVolume();
			//velocityModel = volume.getVelocityModel(); //this.wizard.getVelocityModel();
			//volume = this.wizard.getPreStackVolume();
			//velocityModel =this.wizard.getVelocityModel();
          //  buildBeans();
			//jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	protected abstract void buildBeans();

	public void initialize()
	{
		volume = this.wizard.getPreStackVolume();
		velocityModel =this.wizard.getVelocityModel();
		
        buildBeans();
        
		defineGroupBox.setPanelObject(volume);
		defineGroupBox.gbc.fill = gbc.HORIZONTAL;
		
		if (!initialized)
		{
			jbInit();
			initialized = true;
		}
		
		gbc.fill = gbc.HORIZONTAL;
		add(defineGroupBox);
	}

	abstract void jbInit(); //throws Exception;
}
