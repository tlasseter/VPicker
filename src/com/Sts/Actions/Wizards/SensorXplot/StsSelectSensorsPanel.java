package com.Sts.Actions.Wizards.SensorXplot;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;

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
    private StsSensorXplotWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsComboBoxFieldBean primeSensorBean;
    StsGroupBox primeBox = new StsGroupBox("Select Sensor to Crossplot");
    
    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsSensorXplotWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        if(!buildPanel())
        	wizard.cancel();
    }

    public boolean buildPanel()
    {
    	StsDynamicSensorClass sensorClass = (StsDynamicSensorClass)model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
    	StsObject[] sensors = sensorClass.getSensors();
        primeBox.gbc.fill = gbc.HORIZONTAL;
        primeSensorBean = new StsComboBoxFieldBean(wizard, "primeSensor", "Primary Sensor:", sensors);
        primeBox.addEndRow(primeSensorBean);
    	addEndRow(primeBox);
    	return true;
    }

    public void initialize()
    {
    	;
    }
}