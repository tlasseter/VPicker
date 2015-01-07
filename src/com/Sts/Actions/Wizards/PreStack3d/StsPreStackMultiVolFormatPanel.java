package com.Sts.Actions.Wizards.PreStack3d;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;

import java.awt.*;

/**
 * <p>Title: </p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2006</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author John Friar
 * @version 1.0
 */

public class StsPreStackMultiVolFormatPanel extends StsJPanel
{
	private StsPreStackWizard wizard;
	private StsModel model;

    private StsViewHeadersPanel viewHeadersPanel = null;

    private StsGroupBox parametersGroupBox = new StsGroupBox("Multi-PreStack3d Header Parameters");
    private StsComboBoxFieldBean typeBean;
    private StsFloatFieldBean incBean;
    private StsIntFieldBean numCdpsBean;
    private StsIntFieldBean numVolsBean;
    private StsListFieldBean timeValueList;
    private StsBooleanFieldBean ignoreMultiVolumeBean;

	public StsPreStackMultiVolFormatPanel(StsWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = (StsPreStackWizard)wizard;
		model = wizard.getModel();
		viewHeadersPanel = new StsViewHeadersPanel(wizard.frame);
		try
		{
			constructBeans();
			jbInit();
		}
		catch (Exception ex)
		{
			ex.printStackTrace();
		}
	}

	private void constructBeans()
	{
            typeBean = new StsComboBoxFieldBean(StsPreStackWizard.class, "multiVolType", "PreStack3d Type:", StsSEGYFormat.multiVolumeTypes);
            typeBean.setEditable(false);
            incBean = new StsFloatFieldBean(StsPreStackWizard.class, "velocityInc", false, "Velocity Increment:");
            numCdpsBean = new StsIntFieldBean(StsPreStackWizard.class, "numCdps", false, "CDP's per PreStack3d:");
            numVolsBean = new StsIntFieldBean(StsPreStackWizard.class, "numVolumes", false, "Number of Volumes:");
            ignoreMultiVolumeBean = new StsBooleanFieldBean(StsPreStackWizard.class, "ignoreMultiVolume", "Only Process the First PreStack3d");
            timeValueList = new StsListFieldBean(StsPreStackWizard.class, "timeValuePairs", "Time-Value Pairs:", null);
	}

	private void jbInit() throws Exception
	{
		this.setLayout(new GridBagLayout());

        parametersGroupBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        parametersGroupBox.gbc.anchor = gbc.NORTH;
        parametersGroupBox.add(typeBean);
        parametersGroupBox.add(incBean);
        parametersGroupBox.add(numCdpsBean);
        parametersGroupBox.add(numVolsBean);
        parametersGroupBox.add(timeValueList);

        addEndRow(ignoreMultiVolumeBean);
		gbc.fill = GridBagConstraints.HORIZONTAL;
        addEndRow(parametersGroupBox);
        gbc.anchor = gbc.SOUTH;
        addEndRow(viewHeadersPanel);
	}

	public void initialize()
	{
        parametersGroupBox.setPanelObject(wizard);
        timeValueList.setListItems(wizard.getTimeValuePairs());
	}
}
