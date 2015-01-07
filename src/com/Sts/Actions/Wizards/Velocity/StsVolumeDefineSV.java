package com.Sts.Actions.Wizards.Velocity;

import java.awt.Dimension;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.WizardHeaders.StsHeaderPanel;
import com.Sts.DBTypes.StsObjectList;
import com.Sts.DBTypes.StsWell;
import com.Sts.Actions.Import.StsWellKeywordIO;
import com.Sts.DBTypes.StsLogCurve;
import com.Sts.Utilities.StsParameters;
import com.Sts.UI.StsMessage;
import com.Sts.DBTypes.StsLogVector;

/**
 * <p>Title: jS2S development</p>
 *
 * <p>Description: </p>
 *
 * <p>Copyright: Copyright (c) 2004</p>
 *
 * <p>Company: S2S Systems LLC</p>
 *
 * @author J.B.West
 * @version c51c
 */
public class StsVolumeDefineSV extends StsWizardStep
{
	StsVolumeDefineSVPanel panel;
	StsHeaderPanel header, info;

	/**
	 * StsEditVelocity creates or modifies the velocity model used for
	 * time-to-depth conversion.
	 *
	 * @param wizard StsSeisVelWizard of which this is a step.
	 */
	public StsVolumeDefineSV(StsWizard wizard)
	{
		super(wizard);
		panel = new StsVolumeDefineSVPanel(this);
		header = new StsHeaderPanel();
		setPanels(panel, header);
        panel.setPreferredSize(new Dimension(400, 400));
		header.setTitle("Velocity PostStack3d Definition");
		header.setSubtitle("Define Selected PostStack3d");
//        header.setLink("http://www.s2ssystems.com/marketing/s2ssystems/WellLoad.html");
		header.setInfoText(wizardDialog,"(1) Specify type (average or instantaneous) of velocity in volume.\n"
                    + "(2) Review the minimum and maximum velocity detected in the volume.\n"
                    + "(3) Select the appropriate velocity units for volume data.\n"
                    + "(4) Specify whether one-way or two-way travel time.\n"
                    + "(5) Press Next> Button to proceed to model construction.");
	}
	public boolean start()
    {
        panel.initializePanel();
		setTimeDepthDatum();
        return true;
    }
	private void setTimeDepthDatum()
	{
		float timeDatum = 0.0f;
		panel.topTimeDatumBean.setValue(timeDatum);
		((StsSeisVelWizard)wizard).setTopTimeDatum(timeDatum);
		StsObjectList wells = model.getInstances(StsWell.class);
		int nWells = wells.getSize();
		float minDepthDatum = StsParameters.largeFloat;
		float maxDepthDatum = -StsParameters.largeFloat;

		for(int n = 0; n < nWells; n++)
		{
			StsWell well = (StsWell)wells.getElement(n);
			StsLogCurve tdCurve = well.getLastLogCurveOfType(StsWellKeywordIO.TIME);
			if(tdCurve == null) continue;
			StsLogVector timeVector = tdCurve.getValueVector();
			float indexF = timeVector.getIndexF(0.0f);
			if(indexF == StsParameters.nullValue) continue;
			StsLogVector depthVector = tdCurve.getDepthVector();
			float depth = depthVector.getValue(indexF);
			minDepthDatum = Math.min(minDepthDatum, depth);
			maxDepthDatum = Math.max(maxDepthDatum, depth);
		}
		float depthDatum;
		if(minDepthDatum == StsParameters.largeFloat)
		  depthDatum = 0.0f;
		else
		{
			if(maxDepthDatum - minDepthDatum > 1.0f)
			{
				StsMessage.printMessage("Depth datum values computed from all td curves range from " + minDepthDatum + " to " + maxDepthDatum + ".\n" +
							   "Set Top depth datum manually.");
			}
			depthDatum = (maxDepthDatum + minDepthDatum)/2;
		}
		panel.topDepthDatumBean.setValue(depthDatum);
		((StsSeisVelWizard)wizard).setTopDepthDatum(depthDatum);
   }

    public boolean end()
    {
        return true;
    }
}
