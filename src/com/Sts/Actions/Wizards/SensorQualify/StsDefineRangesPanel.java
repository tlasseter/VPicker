package com.Sts.Actions.Wizards.SensorQualify;

import com.Sts.Actions.Import.StsSensorImport;
import com.Sts.Actions.Import.StsSensorKeywordIO;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.StsMath;

import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsDefineRangesPanel extends StsJPanel implements ChangeListener
{
    private StsSensorQualifyWizard wizard;
    private StsDefineRanges wizardStep;
 
    private StsObject[] commonCurves = null;
    private StsModel model = null;
    StsGroupBox criteriaBox = new StsGroupBox("Define Qualification Criteria");
    StsRangeSliderBean[] rangeBeans = null;
    StsButton[] rangeEditBtns = null;

    int selectedIndex = 0;
    
	StsGroupBox analyzeBox = new StsGroupBox("Results");
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton exportButton;
	StsButton clusterButton;
    StsButton saveFilterButton;
	
	StsProgressPanel progressPanel;
	
    public StsDefineRangesPanel(StsSensorQualifyWizard wizard, StsDefineRanges wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {	
    	progressPanel = StsProgressPanel.constructorWithCancelButton(); 

        criteriaBox.gbc.fill = gbc.HORIZONTAL;
        criteriaBox.gbc.gridwidth = 2;
        
        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Adjust range sliders to desired settings");
        StsJPanel msgPanel = new StsJPanel();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);
        
        StsJPanel btnPanel = new StsJPanel();
        exportButton = new StsButton("Export View", "Export the sensor events currently in view.", wizard, "exportView", progressPanel);        
        clusterButton = new StsButton("Save Clusters", "Save the current cluster values as a an attribute within the respective sensors", wizard, "saveClusters");        
        saveFilterButton = new StsButton("Save Filter", "Save the current attribute ranges as a filter to be selectively applied to sensors", this, "saveAsFilter");

        btnPanel.gbc.fill = gbc.NONE;
        btnPanel.addToRow(exportButton);
        btnPanel.addToRow(clusterButton);
        btnPanel.addEndRow(saveFilterButton);
        
        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(btnPanel);
        gbc.fill = gbc.BOTH;
        gbc.weighty = 1.0;
        add(criteriaBox);
        gbc.fill = gbc.HORIZONTAL;
        gbc.weighty = 0.0f;
        add(analyzeBox);
    }

    public void initialize()
    {
    	criteriaBox.removeAll();
    	
    	commonCurves = wizard.getCommonAttributes();
    	rangeBeans = new StsRangeSliderBean[commonCurves.length];
        rangeEditBtns = new StsButton[commonCurves.length];
    	for(int i=0; i<commonCurves.length; i++)
    	{
    		StsTimeCurve curve = (StsTimeCurve)commonCurves[i];
    		float min = getMinValue(curve.getName());
    		float max = getMaxValue(curve.getName());
    		rangeBeans[i] = new StsRangeSliderBean(min, max, min, true, StsRangeSliderBean.HORIZONTAL);
    		rangeBeans[i].setThumbColor(new Color(30,140,0,120));
    		rangeBeans[i].addChangeListener(this);

            rangeEditBtns[i] = new StsButton("Edit Bounds", "Set range outside data limits.", this, "editRange" + i);

    		JLabel label = new JLabel(curve.getName());
    		criteriaBox.gbc.fill = gbc.NONE;
    		criteriaBox.gbc.anchor = gbc.EAST;
    		criteriaBox.addToRow(label);
    		criteriaBox.gbc.fill = gbc.HORIZONTAL;
    		criteriaBox.gbc.anchor = gbc.WEST;    		
    		criteriaBox.addToRow(rangeBeans[i]);

            criteriaBox.gbc.fill = gbc.NONE;
            criteriaBox.gbc.anchor = gbc.EAST;
            criteriaBox.addEndRow(rangeEditBtns[i]);
    	}
    	wizard.rebuild();
    }
    
    private float getMinValue(String curveName)
    {
    	Object[] sensors = (Object[])wizard.getSelectedSensors();
    	if(sensors.length == 0) 
    		return 0.0f;
        float min = ((StsSensor)sensors[0]).getTimeCurve(curveName).getMinValue() + (float)((StsSensor)sensors[0]).getTimeCurve(curveName).getValueVector().getOrigin();
    	for(int i=1; i<sensors.length; i++)
    	{
            StsTimeCurve curve = ((StsSensor)sensors[i]).getTimeCurve(curveName);
            float value = curve.getMinValue() + (float)curve.getValueVector().getOrigin();
    		if(value < min)
    			min = value;
    	}
    	return min;	
    }
    private float getMaxValue(String curveName)
    {
    	Object[] sensors = (Object[])wizard.getSelectedSensors();
    	if(sensors.length == 0) 
    		return 0.0f;
    	float max = ((StsSensor)sensors[0]).getTimeCurve(curveName).getMaxValue() + (float)((StsSensor)sensors[0]).getTimeCurve(curveName).getValueVector().getOrigin();
    	for(int i=1; i<sensors.length; i++)
    	{
            StsTimeCurve curve = ((StsSensor)sensors[i]).getTimeCurve(curveName);
            float value = curve.getMaxValue() + (float)curve.getValueVector().getOrigin();
    		if(value > max)
    			max = value;
    	}
    	return max;	
    }    
    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }

    public void stateChanged(ChangeEvent ev) 
    {
    	Object source = ev.getSource();
    	Object[] sensors = (Object[])wizard.getSelectedSensors();
    	for(int j=0; j<sensors.length; j++)
    		((StsSensor)sensors[j]).setClusters(null);
    	
    	for(int i=0; i<rangeBeans.length; i++)
    	{
    	    for(int j=0; j<sensors.length; j++)
    	    {
                float lowVal = rangeBeans[i].getFloatLowValue();
                float highVal = rangeBeans[i].getFloatHighValue();
    	    	((StsSensor)sensors[j]).setAttributeLimits(((StsTimeCurve)commonCurves[i]), lowVal, highVal);
    	    	wizard.getModel().viewObjectRepaint(this, (StsSensor)sensors[j]);
    	    } 
    	}
    }

    public boolean saveAsFilter()
    {
        int numRanges = 0;
        StsSensorQualifyFilter filter = null;

        for(int i=0; i<commonCurves.length; i++)
    	{
            String curveName = commonCurves[i].getName();
            float lowVal = rangeBeans[i].getFloatLowValue();
            float highVal = rangeBeans[i].getFloatHighValue();
            if((lowVal <= rangeBeans[i].getFloatMinimum()) && (highVal >= rangeBeans[i].getFloatMaximum()))
                  continue;

            numRanges++;
            if(filter == null)
                 filter = new StsSensorQualifyFilter("QualifyFilter");

            if(!filter.addRange(curveName, lowVal, highVal))
                numRanges--;
    	}

        if(numRanges > 0)
        {
            msgBean.setText("Qualify filter has been successfully saved.");
            return true;
        }
        else
        {
            msgBean.setText("Failed to create filter. Check ranges and try again.");
            return false;
        }
    }

    public float getMinimum()
    {
        return rangeBeans[selectedIndex].getFloatMinimum();
    }
    public void setMinimum(float min)
    {
        rangeBeans[selectedIndex].setFloatMinimum(min);
    }
    public float getMaximum()
    {
        return rangeBeans[selectedIndex].getFloatMaximum();
    }
    public void setMaximum(float max)
    {
        rangeBeans[selectedIndex].setFloatMaximum(max);
    }
    public void editRange(int index)
    {
        selectedIndex = index;

        // Dialog to get user range
        StsJPanel panel = new StsJPanel();
        StsFloatFieldBean minBean = new StsFloatFieldBean(this, "minimum", true, "Minimum:", false);
        StsFloatFieldBean maxBean = new StsFloatFieldBean(this, "maximum", true, "Maximum:", false);

        panel.addToRow(minBean);
        panel.addEndRow(maxBean);

        JDialog dialog = new JDialog(wizard.frame, commonCurves[index].getName() + " Range Override", true);
        dialog.setSize(400, 150);
        dialog.getContentPane().add(panel);
        dialog.setVisible(true);

        // Reset ranges and UI
        wizard.rebuild();
    }

    public void editRange0() { editRange(0); }
    public void editRange1() { editRange(1); }
    public void editRange2() { editRange(2); }
    public void editRange3() { editRange(3); }
    public void editRange4() { editRange(4); }
    public void editRange5() { editRange(5); }
    public void editRange6() { editRange(6); }
    public void editRange7() { editRange(7); }
    public void editRange8() { editRange(8); }
    public void editRange9() { editRange(9); }
    public void editRange10() { editRange(10); }
    public void editRange11() { editRange(11); }
    public void editRange12() { editRange(12); }
    public void editRange13() { editRange(13); }
    public void editRange14() { editRange(14); }
    public void editRange15() { editRange(15); }
    public void editRange16() { editRange(16); }
    public void editRange17() { editRange(17); }
    public void editRange18() { editRange(18); }
    public void editRange19() { editRange(19); }
}
