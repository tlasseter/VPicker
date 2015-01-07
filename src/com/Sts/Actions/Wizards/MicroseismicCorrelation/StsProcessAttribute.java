package com.Sts.Actions.Wizards.MicroseismicCorrelation;

import com.Sts.Actions.Wizards.*;
import com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard;
import com.Sts.Actions.Wizards.WizardHeaders.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.StsCursor3d;
import com.Sts.UI.Progress.*;
import com.Sts.Utilities.*;
import com.Sts.DBTypes.StsSeismicVolume;
import com.Sts.DBTypes.StsDynamicSensor;

public class StsProcessAttribute extends StsWizardStep implements Runnable
{
	public StsProgressPanel panel;
    private StsHeaderPanel header;
	private boolean canceled = false;

    public StsProcessAttribute(StsMicroseismicCorrelationWizard wizard)
    {
		super(wizard);
		this.wizard = wizard;
		panel = StsProgressPanel.constructor(false,10,40);
		header = new StsHeaderPanel();
		setPanels(panel, header);
		header.setTitle("Compute Correlation Attributes");
		header.setSubtitle("Correlating Volume to Microseismic Events");
        header.setLink("http://www.s2ssystems.com/Protected/Marketing/AppLinks.html#MicroseismicCorrelation");                                        
		header.setInfoText(wizardDialog, "(1) Once complete, press the Finish Button to dismiss the screen");
    }

	public boolean start()
	{
		run();
		return true;
	}

    public void run()
    {
		if(canceled)
		{
			success = false;
			return;
		}
		StsModel model = wizard.getModel();
		StsProject project = model.getProject();
        String name = ((StsMicroseismicCorrelationWizard)wizard).getAttributeName();
		try
		{
			disablePrevious();
			model.disableDisplay();
			panel.appendLine("Creating sensor attribute:" + name);

            // Compute Attribute
            StsSeismicVolume volume = ((StsMicroseismicCorrelationWizard)wizard).getSelectedVolume();
            if(!project.hasVelocityModel() && (volume.getZDomain() != StsProject.TD_DEPTH))
            {
                statusArea.setText("Selected volume is not available in depth domain: " + name);
                panel.appendLine("Selected volume is not available in depth domain.\n");
                panel.appendLine("New attribute {" + name + ") computation unsuccessfully.");
                success = false;
                panel.finished();
                return;
            }
            boolean isInputFloats = true;
            if(!volume.setupReadRowFloatBlocks())
            {
                volume.setupRowByteBlocks();
                isInputFloats = false;
            }
            Object[] sensors = (Object[])((StsMicroseismicCorrelationWizard)wizard).getSelectedSensors();
            for(int i=0; i<sensors.length; i++)
            {
                StsDynamicSensor sensor = (StsDynamicSensor)sensors[i];
                float[][] xyzs = sensor.getXYZVectors();
                float[] values = new float[sensor.getNumValues()];
                for(int j=0; j<values.length; j++)
                {
                    float row = volume.getRowCoor(xyzs[0][j]);
                    float col = volume.getColCoor(xyzs[1][j]);
                    float slice = volume.getSliceCoor(xyzs[2][j]);
                    values[j] = computeValue(volume, row, col, slice, isInputFloats);
                }
                sensor.setAttribute(values, name);
                sensor.saveAttribute();
            }
            statusArea.setText("Processing successfully completed for:" + name);
            panel.appendLine("New attribute " + name + " successfully created.");
            success = true;
            panel.finished();
		}
		catch(Exception e)
		{
            StsException.outputWarningException(this, "run", e);
            statusArea.setText("Processing failed for:" + name);
			success = false;
			return;
		}
        finally
        {
            wizard.completeLoading(success);
        }
    }

    // Compute value - assumes 3x3x3 analysis size
    // TODO: make size user specified.
    private float computeValue(StsSeismicVolume volume, float rowF, float colF, float sliceF, boolean floatVolume)
    {
        float[] values = new float[27];
        int row = Math.round(rowF);
        int col = Math.round(colF);
        int slice = Math.round(volume.getSliceCoor(sliceF));
        int cnt = 0;
        for(int i=row-1; i<=row+1; i++)
        {
            for(int j=col-1; j<=col+1; j++)
            {
                for(int k=slice-1; k<=slice+1; k++)
                {
                    if(floatVolume)
                        values[cnt++] = volume.getRowFloatBlockValue(row,col,slice);
                    else
                        values[cnt++] = volume.getValue(row, col, slice);
                }
            }
        }

        int op = ((StsMicroseismicCorrelationWizard)wizard).getOperator();
        switch(op)
        {
            case StsComputeSetupPanel.AVG:
                return StsMath.average(values, StsParameters.nullValue);
            case StsComputeSetupPanel.MAX:
                return StsMath.maxExcludeNull(values);
            case StsComputeSetupPanel.MIN:
                return StsMath.minExcludeNull(values);
            case StsComputeSetupPanel.VAL:
                return values[13];
            default:
                break;
        }
        return wizard.getModel().getProject().getLogNull();
    }

	public boolean end()
	{
		return true;
	}
}