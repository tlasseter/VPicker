package com.Sts.Actions.Wizards.GriddedSensorAtts;

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
    private StsGriddedSensorAttsWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsListFieldBean sensorListBean;
    private JScrollPane listsScrollPane = new JScrollPane();
    private StsJPanel sensorPanel = new StsJPanel();

    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsGriddedSensorAttsWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsObject[] sensors = (StsObject[])model.getObjectList(StsDynamicSensor.class);
    	if(sensors == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no sensors available to analyze.\n" + 
    				"Please load some before attempting to construct a sensor virtual volume.");
    		return;
    	}
        sensorListBean = new StsListFieldBean(wizard, "selectedSensors", null, sensors);
        sensorPanel.gbc.fill = gbc.BOTH;
        sensorPanel.add(sensorListBean);

        gbc.fill = gbc.BOTH;
        add(listsScrollPane);
        listsScrollPane.getViewport().add(sensorPanel, null);

		//add(sensorListBean);
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
