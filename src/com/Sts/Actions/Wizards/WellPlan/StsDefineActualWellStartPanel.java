package com.Sts.Actions.Wizards.WellPlan;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2005</p>e
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version 1.0
 */

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

public class StsDefineActualWellStartPanel extends StsFieldBeanPanel
{
	private StsDefineWellPlanWizard wizard;
	private StsWellPlan plannedWell;
	private StsWellPlan planSet;
	private StsWell actualWell;

	private StsPoint[] points = null;
	StsProject project = null;

	StsStringFieldBean xyUnitsBean = new StsStringFieldBean("Horizontal units:");
	StsStringFieldBean zUnitsBean = new StsStringFieldBean("Vertical units:");

	StsGroupBox wellBox = new StsGroupBox("Well start description");
	StsFloatFieldBean wellXBean = new StsFloatFieldBean(StsWellPlan.class, "wellX", true, "Well start X:");
	StsFloatFieldBean wellYBean = new StsFloatFieldBean(StsWellPlan.class, "wellY", true, "Well start Y:");
	StsFloatFieldBean wellZBean = new StsFloatFieldBean(StsWellPlan.class, "wellZ", true, "Well start Z:");
	StsFloatFieldBean wellTBean = new StsFloatFieldBean(StsWellPlan.class, "wellT", false, "Well Start Time:");
	StsFloatFieldBean wellRateBean = new StsFloatFieldBean(StsWellPlan.class, "wellRate", true, "Well start rate (deg/100):");

	public StsDefineActualWellStartPanel(StsWizard wizard, StsWell actualWell)
	{
		this.wizard = (StsDefineWellPlanWizard) wizard;
		this.actualWell = actualWell;
		project = wizard.getModel().getProject();
		try
		{
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	public void initialize(StsWellPlan planSet)
	{
		this.planSet = planSet;
		wellBox.setPanelObject(planSet); // connect beans to this object
		xyUnitsBean.setValue(project.getXyUnitString());
		zUnitsBean.setValue(project.getDepthUnitString());
	}

	void jbInit() throws Exception
	{
//        xyUnitsBean.classInitialize("Horizontal units:");
//        zUnitsBean.classInitialize("Vertical units:");

		wellBox.add(xyUnitsBean);
		wellBox.add(zUnitsBean);
		wellBox.add(wellXBean);
		wellBox.add(wellYBean);
		wellBox.add(wellZBean);
		wellBox.add(wellTBean);
		wellBox.add(wellRateBean);
		add(wellBox);
	}
}
