package com.Sts.Actions.Wizards.MicroseismicCorrelation;

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

public class StsSelectSensorsPanel extends StsJPanel
{
    private StsMicroseismicCorrelationWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsListFieldBean sensorListBean;

    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsMicroseismicCorrelationWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsDynamicSensor.class);
    	if(sensors == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no dynamic sensors available to analyze.\n" +
    				"Please load some before attempting to correlate microseismic events with volume");
    		return;
    	}
        sensorListBean = new StsListFieldBean(wizard, "selectedSensors", "Sensors:", sensors);
    	gbc.fill = gbc.BOTH;
		add(sensorListBean);
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