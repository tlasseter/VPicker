package com.Sts.Actions.Wizards.FractureTrack;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;

import javax.swing.*;
import java.awt.*;

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
    private StsFractureTrackWizard wizard;
    private StsSelectSensors wizardStep;

    private StsModel model = null;

    StsListFieldBean sensorListBean;
    StsComboBoxFieldBean wellBean;
	JScrollPane sensorScrollPane = new JScrollPane();    
	
    public StsSelectSensorsPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsFractureTrackWizard)wizard;
    	this.wizardStep = (StsSelectSensors)wizardStep;
    	model  = wizard.getModel();
        buildPanel();
    }

    public void buildPanel()
    {
    	StsDynamicSensorClass sensorClass = (StsDynamicSensorClass)model.getStsClass("com.Sts.DBTypes.StsDynamicSensor");
    	StsObject[] sensors = sensorClass.getSensors();
    	StsObject[] wells = model.getObjectList(StsWell.class);
    	if(sensors == null)
    	{
    		new StsMessage(wizard.frame,StsMessage.ERROR, "There are no sensors available to analyze.\n" + 
    				"Please load some before attempting to run this wizard again.");
    		return;
    	}
        sensorListBean = new StsListFieldBean(wizard, "selectedSensors", null, sensors);
        sensorScrollPane.getViewport().add(sensorListBean, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
                , GridBagConstraints.WEST, GridBagConstraints.BOTH, new Insets(1, 1, 1, 1), 0, 0));
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0;
        wellBean = new StsComboBoxFieldBean(wizard, "selectedWell", "Well:", wells);
    	addEndRow(wellBean);
        gbc.gridwidth = 2;
        JLabel label = new JLabel("Available Sensors");
        addEndRow(label);  
        gbc.fill = gbc.BOTH;
        gbc.gridwidth = 2;
        gbc.weighty = 1.0;        
		addEndRow(sensorScrollPane);
		wizard.rebuild();
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
