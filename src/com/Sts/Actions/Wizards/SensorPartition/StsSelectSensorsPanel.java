package com.Sts.Actions.Wizards.SensorPartition;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.DBTypes.StsObject;
import com.Sts.DBTypes.StsSensorClass;
import com.Sts.MVC.StsModel;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsListFieldBean;

import javax.swing.*;

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
    private StsSensorPartitionWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsListFieldBean sensorListBean;
	JScrollPane sensorScrollPane = new JScrollPane();

    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsSensorPartitionWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsSensorClass sensorClass = (StsSensorClass)model.getStsClass("com.Sts.DBTypes.StsSensor");
    	StsObject[] sensors = sensorClass.getSensors();

        sensorListBean = new StsListFieldBean(wizard, "selectedSensors", null, sensors);
        sensorScrollPane.getViewport().add(sensorListBean, null);

    	gbc.fill = gbc.HORIZONTAL;
    	gbc.weighty = 0.0;
        JLabel label = new JLabel("Available Sensors");
        addEndRow(label);
    	gbc.fill = gbc.BOTH;
    	gbc.weighty = 1.0;
		addEndRow(sensorScrollPane);
    }

    public void initialize()
    {
    	sensorListBean.setSingleSelect();
        sensorListBean.setSelectedIndex(0);
    }

    public Object[] getSelectedSensors()
    {
    	return sensorListBean.getSelectedObjects();
    }
}