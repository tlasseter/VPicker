package com.Sts.Actions.Wizards.ClusterAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

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
    private StsClusterAnalysisWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsListFieldBean sensorListBean;
	JScrollPane sensorScrollPane = new JScrollPane();    
	
    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsClusterAnalysisWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsDynamicSensorClass sensorClass = (StsDynamicSensorClass)model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
    	StsObject[] sensors = sensorClass.getSensors();
    	if(sensors == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no sensors available to analyze.\n" + 
    				"Please load some before attempting to run this wizard again.");
    		return;
    	}
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
    	;
    }

    public Object[] getSelectedSensors()
    {
    	return sensorListBean.getSelectedObjects();
    }
}
