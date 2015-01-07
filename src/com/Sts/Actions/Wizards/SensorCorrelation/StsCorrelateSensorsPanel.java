package com.Sts.Actions.Wizards.SensorCorrelation;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Actions.Wizards.SensorCompare.StsSensorCompareWizard;
import com.Sts.Actions.Wizards.SensorCompare.StsCompareSensors;

import javax.swing.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsCorrelateSensorsPanel extends StsJPanel
{
    private StsSensorCorrelationWizard wizard;
    private StsCorrelateSensors wizardStep;

    private StsModel model = null;

    StsGroupBox criteriaBox = new StsGroupBox("Define Correlation Criteria");
    StsComboBoxFieldBean methodBean = new StsComboBoxFieldBean();
    StsIntFieldBean rangeTimeBean;

	StsGroupBox analyzeBox = new StsGroupBox("Analyze");
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	StsButton analyzeButton;

    public StsCorrelateSensorsPanel(StsSensorCorrelationWizard wizard, StsCorrelateSensors wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void initialize()
    {
        methodBean.setSelectedItem(wizard.getCorrelationMethod());
    }

    public void buildPanel()
    {
    	methodBean.initialize(wizard, "correlationMethod", "Correlation Method:", wizard.methodList);
    	rangeTimeBean = new StsIntFieldBean(wizard, "timeRange", 1, 1000, "Time Range(sec):", true);
        rangeTimeBean.setStep(1);
    	rangeTimeBean.setToolTipText("Specify the time range in seconds to analyze.");

        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.gbc.anchor = gbc.WEST;
        criteriaBox.addEndRow(methodBean);
        criteriaBox.addEndRow(rangeTimeBean);

        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = StsJPanel.addInsets();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);

        StsJPanel btnPanel = StsJPanel.addInsets();
        //runButton = new StsButton("Run Analysis", "Run comparison with specified criteria.", wizard, "analyze");
        analyzeButton = new StsButton("Save Attribute", "Save the values as an attribute within the respective sensors", wizard, "saveAttributes");

        btnPanel.gbc.fill = gbc.NONE;
        //btnPanel.addToRow(runButton);
        btnPanel.addEndRow(analyzeButton);

        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);
        gbc.fill = gbc.HORIZONTAL;
        add(criteriaBox);
        add(analyzeBox);
    }

    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }

}