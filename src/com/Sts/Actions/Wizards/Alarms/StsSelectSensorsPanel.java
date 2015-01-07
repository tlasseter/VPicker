package com.Sts.Actions.Wizards.Alarms;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.UI.Sounds.StsSound;
import com.Sts.Utilities.StsMath;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectSensorsPanel extends StsJPanel
{
    private StsAlarmsWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsListFieldBean sensorListBean;

    StsGroupBox alarmPanel;
    StsComboBoxFieldBean alarmTypeBean;
    StsComboBoxFieldBean queueTypeBean;
    StsComboBoxFieldBean soundBean;

    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsAlarmsWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsObject[] sensors = model.getObjectList(StsDynamicSensor.class);
        sensors = (StsObject[])StsMath.arrayAddArray(sensors, model.getObjectList(StsStaticSensor.class));
    	if(sensors == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no sensors available.\n" +
    				"Please load some before attempting to define alarms.");
    		return;
    	}
        sensorListBean = new StsListFieldBean(wizard, "selectedSensors", "Sensors:", sensors);

        gbc.fill = gbc.BOTH;
		add(sensorListBean);

        alarmPanel = new StsGroupBox("Define Alarm");
        alarmTypeBean = new StsComboBoxFieldBean(wizard, "alarmType", "Alarm Type:", StsAlarm.ALARM_TYPE_STRINGS);
        queueTypeBean = new StsComboBoxFieldBean(wizard, "queueType", "Queue Type:", StsAlarm.QUEUE_TYPE_STRINGS);
        soundBean = new StsComboBoxFieldBean(wizard, "soundFile", "Sound File:", StsSound.sounds);

        alarmPanel.addEndRow(alarmTypeBean);
        alarmPanel.addEndRow(queueTypeBean);
        alarmPanel.addEndRow(soundBean);

        gbc.fill = gbc.HORIZONTAL;
        add(alarmPanel);
    }

    public void initialize()
    {
    	;
    }

    public Object[] getSelectedSensors()
    {
    	return sensorListBean.getSelectedObjects();
    }
}