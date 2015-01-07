package com.Sts.Actions.Wizards.FractureTrack;

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

public class StsDefineTrackerPanel extends StsJPanel
{
    private StsFractureTrackWizard wizard;
    private StsDefineTracker wizardStep;
 
    private StsModel model = null;

    StsGroupBox criteriaBox = new StsGroupBox("Define Tracking Criteria");
    StsRangeSliderBean distanceBean;
    StsRangeSliderBean timeBean;
    StsRangeSliderBean azimuthBean;
    StsBooleanFieldBean transparentNotTrack;
    
	StsGroupBox analyzeBox = new StsGroupBox("Track");
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;
	StsButton exportButton;
	StsButton clearButton;
	
    public StsDefineTrackerPanel(StsFractureTrackWizard wizard, StsDefineTracker wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
		JLabel dlabel = new JLabel("Distance (m/ft):");    	
        distanceBean = new StsRangeSliderBean(0, 1000, 0, StsRangeSliderBean.HORIZONTAL);
        distanceBean.setHighValue(250);
        distanceBean.setToolTipText("Specify the maximum distance from a well to accept an event.");
		
        JLabel tlabel = new JLabel("Time (seconds):");
        timeBean = new StsRangeSliderBean(0, 1000, 0, StsRangeSliderBean.HORIZONTAL);
        timeBean.setHighValue(250);
        timeBean.setToolTipText("Specify the maximum time to look for event in fracture.");
        
        JLabel alabel = new JLabel("Azimuth (degrees):");
        azimuthBean = new StsRangeSliderBean(0, 360, 0, StsRangeSliderBean.HORIZONTAL);
        azimuthBean.setHighValue(45);
        azimuthBean.setToolTipText("Specify the azimuth corridor to search.");
        
        transparentNotTrack = new StsBooleanFieldBean(wizard, "showNotTracked", "Do not show events that are not tracked", false);
        transparentNotTrack.setToolTipText("Do not show events that were not included in tracking.");
		
        criteriaBox.gbc.fill = gbc.NONE;
		criteriaBox.gbc.anchor = gbc.EAST;
		criteriaBox.addToRow(dlabel);
		criteriaBox.gbc.fill = gbc.HORIZONTAL;
		criteriaBox.gbc.anchor = gbc.WEST;    		
		criteriaBox.addEndRow(distanceBean);
		
        criteriaBox.gbc.fill = gbc.NONE;
		criteriaBox.gbc.anchor = gbc.EAST;
		criteriaBox.addToRow(tlabel);
		criteriaBox.gbc.fill = gbc.HORIZONTAL;
		criteriaBox.gbc.anchor = gbc.WEST;    		
		criteriaBox.addEndRow(timeBean);
		
        criteriaBox.gbc.fill = gbc.NONE;
		criteriaBox.gbc.anchor = gbc.EAST;
		criteriaBox.addToRow(alabel);
		criteriaBox.gbc.fill = gbc.HORIZONTAL;
		criteriaBox.gbc.anchor = gbc.WEST;    		
		criteriaBox.addEndRow(azimuthBean);
		
		//criteriaBox.gbc.gridwidth = 2;
		//criteriaBox.addEndRow(transparentNotTrack);
		
        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = new StsJPanel();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        
        StsJPanel btnPanel = new StsJPanel();
        // runButton = new StsButton("Run Analysis", "Execute Proximity Analysis with specified criteria.", wizard, "analyze");
        exportButton = new StsButton("Export View", "Export the sensor events currently in view.", wizard, "exportView");        
        clearButton = new StsButton("Clear Tracked", "Clear the currently tracked events.", wizard, "clearTracking");        
        
        btnPanel.gbc.fill = gbc.NONE;
        //btnPanel.addToRow(runButton);
        btnPanel.addToRow(clearButton);
        btnPanel.addEndRow(exportButton);        
        
        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);
        gbc.fill = gbc.HORIZONTAL;
        add(criteriaBox);
        add(transparentNotTrack);
        add(analyzeBox);
    }

    public void initialize()
    {
    	;
    }
    public int[] getTimeRange()
    {
    	return new int[] {timeBean.getLowValue(), timeBean.getHighValue() };
    }
    public int[] getAzimuthRange()
    {
    	return new int[] {azimuthBean.getLowValue(), azimuthBean.getHighValue() };
    }    
    public int[] getDistanceRange()
    {
    	return new int[] {distanceBean.getLowValue(), distanceBean.getHighValue() };
    }
    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }

}
