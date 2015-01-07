package com.Sts.Actions.Wizards.VelocityAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.UI.Beans.*;

public class StsDefineVelocityAnalysisPanel3d extends StsDefineVelocityAnalysisPanel 
{

	public StsDefineVelocityAnalysisPanel3d(StsWizard wizard, StsWizardStep wizardStep) 
	{
		super(wizard, wizardStep);
	}
	
	protected void buildBeans()
	{
            analysisRowNumStartBean = new StsFloatFieldBean(volume, "analysisRowNumStart", true, "Inline Start:", true);
            analysisRowNumStartBean.setToolTipText("Starting inline number for velocity analysis");
            analysisColNumStartBean = new StsFloatFieldBean(volume, "analysisColNumStart", true, "Crossline Start:", true);
            analysisColNumStartBean.setToolTipText("Starting crossline number for velocity analysis");

            analysisRowNumIncBean = new StsFloatFieldBean(volume, "analysisRowNumInc", true, "Inline Increment:", true);
            analysisRowNumIncBean.setToolTipText("Inline increment for velocity analysis");
            analysisColNumIncBean = new StsFloatFieldBean(volume, "analysisColNumInc", true, "Crossline Increment:", true);
            analysisColNumIncBean.setToolTipText("Crossline increment for velocity analysis");

            traceThresholdBean = new StsIntFieldBean(volume, "traceThreshold", true, "Minimum Traces Required:", true);
            traceThresholdBean.setToolTipText("Minimum number of traces required in a gather for velocity analysis");
	}
	
	protected void jbInit() //throws Exception
	{
		analysisRowNumStartBean.setValueAndRangeFixStep(volume.getAnalysisRowNumStart(), volume.getRowNumMin(), volume.getRowNumMax(), volume.getRowNumInc());
		analysisRowNumIncBean.setValueAndRangeFixStep(volume.getAnalysisRowNumInc(), volume.getRowNumInc(), volume.getRowNumMax() - volume.getRowNumMin(), volume.getRowNumInc());
		analysisColNumStartBean.setValueAndRangeFixStep(volume.getAnalysisColNumStart(), volume.getColNumMin(), volume.getColNumMax(), volume.getColNumInc());
		analysisColNumIncBean.setValueAndRangeFixStep(volume.getAnalysisColNumInc(), volume.getColNumInc(), volume.getColNumMax() - volume.getColNumMin(), volume.getColNumInc());

		defineGroupBox.addEndRow(traceThresholdBean);
		defineGroupBox.addEndRow(analysisRowNumStartBean);
		defineGroupBox.addEndRow(analysisRowNumIncBean);
		defineGroupBox.addEndRow(analysisColNumStartBean);
		defineGroupBox.addEndRow(analysisColNumIncBean);
	}

}
