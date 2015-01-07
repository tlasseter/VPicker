package com.Sts.Actions.Wizards.SensorCompare;

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

public class StsCompareSensorsPanel extends StsJPanel
{
    private StsSensorCompareWizard wizard;
    private StsCompareSensors wizardStep;
 
    private StsModel model = null;
    
    StsGroupBox criteriaBox = new StsGroupBox("Define Comparison Criteria");
    StsGroupBox thresholdBox = new StsGroupBox();
    StsGroupBox typeBox = new StsGroupBox();

    StsBooleanFieldBean interactiveModeBean;
    ButtonGroup criteriaGrp = new ButtonGroup();    
    JRadioButton timeBtn = new JRadioButton("Time");
    JRadioButton distanceBtn = new JRadioButton("Distance");
    JRadioButton xyDistanceBtn = new JRadioButton("XY Distance");
    JRadioButton amplitudeBtn = new JRadioButton("Amplitude");
    StsComboBoxFieldBean thresholdAttributeBean = new StsComboBoxFieldBean();
    StsFloatFieldBean thresholdBean;
    
	StsGroupBox analyzeBox = new StsGroupBox("Analyze");
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	StsButton clusterButton;
	
    public StsCompareSensorsPanel(StsSensorCompareWizard wizard, StsCompareSensors wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }
    
    public void buildPanel()
    { 
    	thresholdAttributeBean.initialize(wizard, "thresholdAttribute", "Threshold Attribute:", wizard.attList);
    	thresholdBean = new StsFloatFieldBean(wizard, "threshold", 1, 1000, "Threshold:", true);
    	thresholdBean.setToolTipText("Specify the analysis threshold.");
        
        criteriaGrp.add(timeBtn);
        criteriaGrp.add(distanceBtn);
        criteriaGrp.add(xyDistanceBtn);
        criteriaGrp.add(amplitudeBtn);
        
        typeBox.gbc.fill = gbc.HORIZONTAL;
        typeBox.gbc.anchor = gbc.WEST;
        typeBox.addToRow(timeBtn);
        typeBox.addEndRow(distanceBtn);
        typeBox.addToRow(xyDistanceBtn);
        typeBox.addEndRow(amplitudeBtn);
        thresholdBox.gbc.fill = gbc.HORIZONTAL;
        thresholdBox.addEndRow(thresholdAttributeBean);
        thresholdBox.addEndRow(thresholdBean);

        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.gbc.anchor = gbc.WEST;
        criteriaBox.addEndRow(typeBox);        
        criteriaBox.addEndRow(thresholdBox);

        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = StsJPanel.addInsets();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        
        StsJPanel btnPanel = StsJPanel.addInsets();
        interactiveModeBean = new StsBooleanFieldBean(wizard, "interactiveMode", "Interactive Mode", false);
        interactiveModeBean.setToolTipText("Interactively pick sensor points and run comparison.");
        runButton = new StsButton("Run Analysis", "Run comparison with specified criteria.", wizard, "analyze");
        clusterButton = new StsButton("Save Clusters", "Save the current cluster values as an attribute within the respective sensors", wizard, "saveClusters");
       
        btnPanel.gbc.fill = gbc.NONE;
        //btnPanel.addToRow(interactiveModeBean);
        btnPanel.addToRow(runButton);
        //btnPanel.addEndRow(clusterButton);
        
        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);
        gbc.fill = gbc.HORIZONTAL;
        add(criteriaBox);
        add(analyzeBox);
    }
    
    public void reinitializeThreshold()
    {
        if(thresholdAttributeBean.getSelectedIndex() == wizard.AMPL)
        {
        	float[] minMax = wizard.getAmplitudeRange(wizard.getPrimeSensor());
        	thresholdBean.setValueAndRangeFixStep(minMax[1], minMax[0], minMax[1], (minMax[1] - minMax[0])/100.0f);
        }
        else if(thresholdAttributeBean.getSelectedIndex() == wizard.DIST)
        {
        	thresholdBean.setValueAndRangeFixStep(10000, 0, 10000, 10);
        }
        else if(thresholdAttributeBean.getSelectedIndex() == wizard.XYDIST)
        {
        	thresholdBean.setValueAndRangeFixStep(10000, 0, 10000, 10);
        }
        else if(thresholdAttributeBean.getSelectedIndex() == wizard.TIME)
        {
        	long max = (long)wizard.getPrimeSensor().getTimeMax() - (long)wizard.getPrimeSensor().getTimeMin();
        	thresholdBean.setValueAndRangeFixStep(max, 0, max, (int)(max/1000));
        }
    }
    
    public void initialize()
    {
    	timeBtn.setSelected(true);
    	thresholdAttributeBean.setSelectedIndex(wizard.NONE);   	
    }

    public byte getComparisonType()
    {
    	if(timeBtn.isSelected())
    		return wizard.TIME;
    	else if (distanceBtn.isSelected())
    		return wizard.DIST;
    	else if(xyDistanceBtn.isSelected())
    		return wizard.XYDIST;
    	else if(amplitudeBtn.isSelected())
    		return wizard.AMPL;
    	return wizard.NONE;
    }
    
    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }

}