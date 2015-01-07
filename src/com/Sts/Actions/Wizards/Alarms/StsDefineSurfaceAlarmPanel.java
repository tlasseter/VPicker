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

public class StsDefineSurfaceAlarmPanel extends StsJPanel
{
    private StsAlarmsWizard wizard;
    private StsDefineSurfaceAlarm wizardStep;

    private StsModel model = null;

    StsListFieldBean surfaceListBean;

    StsGroupBox alarmPanel;
    StsFloatFieldBean offsetBean;
    StsComboBoxFieldBean directionBean;

    public StsDefineSurfaceAlarmPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsAlarmsWizard)wizard;
    	this.wizardStep = (StsDefineSurfaceAlarm)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsObject[] surfaces = model.getObjectList(StsSurface.class);
    	if(surfaces == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no surfaces available.\n" +
    				"Please load some before attempting to define alarms.");
    		return;
    	}
        surfaceListBean = new StsListFieldBean(wizard, "selectedSurface", "Surfaces:", surfaces);
        surfaceListBean.setSingleSelect();

        gbc.fill = gbc.BOTH;
		add(surfaceListBean);

        alarmPanel = new StsGroupBox();
        directionBean = new StsComboBoxFieldBean(wizard, "direction", "Reference Direction:", StsSurfaceAlarm.REF_DIRECTION);
        offsetBean = new StsFloatFieldBean(wizard, "Offset", true, "Offset:");
        alarmPanel.addEndRow(directionBean);        
        alarmPanel.addEndRow(offsetBean);

        gbc.fill = gbc.HORIZONTAL;
        add(alarmPanel);
    }

    public void initialize()
    {
    	;
    }

}