package com.Sts.Actions.Wizards.ProximityAnalysis;

import com.Sts.DBTypes.StsObject;
import com.Sts.DBTypes.StsSensor;
import com.Sts.DBTypes.StsSensorClass;
import com.Sts.DBTypes.StsSensorQualifyFilter;
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

public class StsDefineProximityPanel extends StsJPanel
{
    private StsProximityAnalysisWizard wizard;
    private StsDefineProximity wizardStep;
 
    private StsModel model = null;

    StsGroupBox criteriaBox = new StsGroupBox("Define Proximity Criteria");
    StsFloatFieldBean distanceBean;
	StsBooleanFieldBean ignoreVertical;
	StsBooleanFieldBean insideLimits;	
	StsGroupBox analyzeBox = new StsGroupBox("Analyze Proximity");
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	StsButton exportButton;
	StsButton clusterButton;
    StsButton saveFilterButton;
	
    public StsDefineProximityPanel(StsProximityAnalysisWizard wizard, StsDefineProximity wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
        distanceBean = new StsFloatFieldBean(wizard, "distanceLimit", true, "Distance Limit:");
        distanceBean.setToolTipText("Specify the maximum distance from a well to accept an event.");
        ignoreVertical = new StsBooleanFieldBean(wizard, "ignoreVertical", "Compute distance in only XY", false);
        ignoreVertical.setToolTipText("Ignore the vertical component when computing distance.");
        insideLimits = new StsBooleanFieldBean(wizard, "insideLimits", "Accept events inside distance limit", false);
        insideLimits.setToolTipText("Accept events inside or outside specified distance limit.");
        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.add(distanceBean);
        criteriaBox.add(ignoreVertical);
        criteriaBox.add(insideLimits);

        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = StsJPanel.addInsets();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        
        StsJPanel btnPanel = StsJPanel.addInsets();
        runButton = new StsButton("Run Analysis", "Execute Proximity Analysis with specified criteria.", wizard, "analyze");
        exportButton = new StsButton("Export View", "Export the sensor events currently in view.", wizard, "exportView");        
        clusterButton = new StsButton("Save Clusters", "Save the current cluster values as a an attribute within the respective sensors", wizard, "saveClusters");        
        saveFilterButton = new StsButton("Save Filter", "Save the current proximity values as a filter to be selectively applied to sensors", wizard, "saveAsFilter");
        
        btnPanel.gbc.fill = gbc.NONE;
        btnPanel.addToRow(runButton);
        btnPanel.addToRow(exportButton);
        btnPanel.addToRow(clusterButton);
        btnPanel.addEndRow(saveFilterButton);
        
        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);
        gbc.fill = gbc.HORIZONTAL;
        add(criteriaBox);
        add(analyzeBox);
    }

    public void initialize()
    {
    	;
    }
    
    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }

}