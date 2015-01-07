package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.Types.PreStack.*;
import com.Sts.UI.Beans.*;

public class StsDefineVelocityAnalysisPanel2d extends StsDefineVelocityAnalysisPanel 
{
	transient StsComboBoxFieldBean lineSelectBean = null;
	transient StsPreStackLine2d selectedLine;
	transient StsPreStackLine[] lines;
	
	public StsDefineVelocityAnalysisPanel2d(StsWizard wizard, StsWizardStep wizardStep) 
	{
		super(wizard, wizardStep);
	}
	
	protected void buildBeans() 
	{
	    volume.setTraceThreshold(2);  //start off 2D w/ default is 2 instead of 10 (many 2D lines are 6 fold or less)
		lines = volume.lines;
		selectedLine = (StsPreStackLine2d) lines[0];
		lineSelectBean  = new StsComboBoxFieldBean(this, "selectedLine", "2D Data Set:" , lines);
		
		analysisColNumStartBean = new StsFloatFieldBean(selectedLine, "analysisColStart", true, "CDP Start:", true);
    	analysisColNumStartBean.setToolTipText("Starting CDP number for velocity analysis");

        analysisColIncBean = new StsIntFieldBean(selectedLine, "analysisColInc", true, "CDP Increment:", true);
        analysisColIncBean.setToolTipText("CDP increment for velocity analysis");
        
        traceThresholdBean = new StsIntFieldBean(volume, "traceThreshold", true, "Minimum Traces Required:", true);
        traceThresholdBean.setToolTipText("Minimum number of traces required in a gather for velocity analysis");
	}
	
	public StsPreStackLine2d getSelectedLine() {
		return selectedLine;
	}

	public void setSelectedLine(StsPreStackLine2d selectedLine) {
		this.selectedLine = selectedLine;
		analysisColNumStartBean.setBeanObject(selectedLine);
		analysisColIncBean.setBeanObject(selectedLine);
	}

	protected void jbInit() //throws Exception
	{
		analysisColNumStartBean.setValueAndRangeFixStep(selectedLine.getAnalysisColStart(), selectedLine.colNumMin, selectedLine.colNumMax, 1);
        analysisColIncBean.setValueAndRangeFixStep(selectedLine.getAnalysisColInc(), 1, selectedLine.nCols-1, 1);
        
        defineGroupBox.addEndRow(traceThresholdBean);
        defineGroupBox.addEndRow(lineSelectBean);
		defineGroupBox.addEndRow(analysisColNumStartBean);
        defineGroupBox.addEndRow(analysisColIncBean);
	}

}
