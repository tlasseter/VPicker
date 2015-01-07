package com.Sts.Actions.Wizards.Alarms;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineWellAlarmPanel extends StsJPanel
{
    private StsAlarmsWizard wizard;
    private StsDefineWellAlarm wizardStep;

    private StsModel model = null;

    StsListFieldBean wellListBean;

    StsGroupBox alarmPanel;
    StsFloatFieldBean offsetBean;
    StsBooleanFieldBean insideBean;

    public StsDefineWellAlarmPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsAlarmsWizard)wizard;
    	this.wizardStep = (StsDefineWellAlarm)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsObject[] wells = model.getObjectList(StsWell.class);
    	if(wells == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no wells available.\n" +
    				"Please load some before attempting to define alarms.");
    		return;
    	}
        wellListBean = new StsListFieldBean(wizard, "selectedWell", "Wells:", wells);
        wellListBean.setSingleSelect();

        gbc.fill = gbc.BOTH;
		add(wellListBean);

        alarmPanel = new StsGroupBox();
        offsetBean = new StsFloatFieldBean(wizard, "Offset", true, "Radius:");
        offsetBean.setToolTipText("Specify the radius around the well.");

        insideBean = new StsBooleanFieldBean(wizard, "inside", true, "Inside:", true);
        insideBean.setToolTipText("Alarm sounds if inside radius, or outside radius.");

        alarmPanel.addEndRow(offsetBean);
        alarmPanel.addEndRow(insideBean);

        gbc.fill = gbc.HORIZONTAL;
        add(alarmPanel);
    }

    public void initialize()
    {
    	;
    }

}