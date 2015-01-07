package com.Sts.Actions.Wizards.ClusterAnalysis;

import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineClusteringPanel extends StsJPanel
{
    private StsClusterAnalysisWizard wizard;
    private StsDefineClustering wizardStep;
 
    private StsModel model = null;

    StsProgressPanel progressPanel;
    StsGroupBox criteriaBox = new StsGroupBox("Define Cluster Criteria");
    StsFloatFieldBean distanceBean;
	StsBooleanFieldBean ignoreVertical;    
    StsIntFieldBean timeBean;
    StsIntFieldBean minimumSizeBean;
    
	StsGroupBox relateBox = new StsGroupBox("Define Amplitude Relationships");
	StsFloatFieldBean amplitudeBean;
	StsBooleanFieldBean relateTime;
	StsBooleanFieldBean relateDistance;

	StsGroupBox analyzeBox = new StsGroupBox("Analyze Clusters");
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	StsButton exportButton;
	StsButton clusterButton;	
	
    public StsDefineClusteringPanel(StsClusterAnalysisWizard wizard, StsDefineClustering wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
    	progressPanel = StsProgressPanel.constructorWithCancelButton();
    	
        minimumSizeBean = new StsIntFieldBean(wizard, "minimumClusterSize", 2, 10, "Minimum Cluster Size:", true);
        minimumSizeBean.setToolTipText("Specify the minimum number of events within a cluster.");
        distanceBean = new StsFloatFieldBean(wizard, "distanceLimit", true, "Distance Limit:");
        distanceBean.setToolTipText("Specify the maximum distance between related events. -1 indicates no distance limit");
        ignoreVertical = new StsBooleanFieldBean(wizard, "ignoreVertical", "Compute distance in only XY", false);
        ignoreVertical.setToolTipText("Ignore the vertical component when computing distance.");
        timeBean = new StsIntFieldBean(wizard,"timeLimit", true, "Time Limit(ms):");
        timeBean.setToolTipText("Specify the maximum time in milliseconds between related events. -1 indicates no distance limit.");
        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.gbc.gridwidth = 2;
        criteriaBox.add(distanceBean);
        criteriaBox.gbc.gridwidth = 1;
        criteriaBox.add(ignoreVertical);
        criteriaBox.gbc.gridwidth = 2;
        criteriaBox.add(timeBean);
        criteriaBox.add(minimumSizeBean);
        
        amplitudeBean = new StsFloatFieldBean(wizard, "amplitude", true, "Amplitude:");
        amplitudeBean.setToolTipText("Specify the amplitude that corresponds with the supplied time and distanc limits.");
        relateTime = new StsBooleanFieldBean(wizard, "relateTime", "Linearly relate time with amplitude", false);
        relateTime.setToolTipText("Check to linearly relate a change in amplitude with a change in time criteria");
        relateDistance = new StsBooleanFieldBean(wizard, "relateDistance", "Linearly relate distance with amplitude", false);
        relateTime.setToolTipText("Check to linearly relate a change in amplitude with a change in distance criteria");        
        relateBox.gbc.fill = gbc.HORIZONTAL;
        relateBox.gbc.gridwidth = 1;
        relateBox.add(amplitudeBean);
        relateBox.gbc.gridwidth = 2;
        relateBox.add(relateTime);
        relateBox.add(relateDistance);
        
        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = StsJPanel.addInsets();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        
        StsJPanel btnPanel = StsJPanel.addInsets();
        runButton = new StsButton("Run Analysis", "Execute Cluster Analysis with specified criteria.", wizard, "analyze", progressPanel);
        exportButton = new StsButton("Export View", "Export the sensor events currently in view.", wizard, "exportView");        
        clusterButton = new StsButton("Save Clusters", "Save the current cluster values as a an attribute within the respective sensors", wizard, "saveClusters");        

        btnPanel.gbc.fill = gbc.NONE;
        btnPanel.addToRow(runButton);
        btnPanel.addToRow(exportButton);
        btnPanel.addEndRow(clusterButton);
        
        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(progressPanel);
        analyzeBox.addEndRow(btnPanel);
        gbc.fill = gbc.HORIZONTAL;
        add(criteriaBox);
        add(relateBox);
        add(analyzeBox);
    }

    public void initialize()
    {

    }
    
    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }

}