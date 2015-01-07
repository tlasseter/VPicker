package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.SeismicAttributes.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import java.awt.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsProcessVolumePanel extends StsJPanel
{
    private StsVolumeCurvatureWizard wizard;
    private StsProcessVolume wizardStep;

    public int filterSize = 3;
    public int interpSize = 5;
    public float curvatureMin;
    public float curvatureMax;
    public boolean runAllPatches = true;
    public byte attributeType = StsSurfaceCurvatureAttribute.CURVPos;

    StsGroupBox curvatureBox = new StsGroupBox("Create Curvature Attribute");
    StsGroupBox editRangeBox = new StsGroupBox("Edit Colorscale");
	StsGroupBox outputVolumeBox = new StsGroupBox("Output Volume");
	StsEditableColorscaleFieldBean colorscaleBean = null;
    StsIntFieldBean filterSizeBean;
    StsBooleanFieldBean runAllBean;
    StsIntFieldBean interpSizeBean;
    StsComboBoxFieldBean curvAttributeBean;
    StsFloatFieldBean curvatureMinBean;
    StsFloatFieldBean curvatureMaxBean;
	StsButton runButton;
	StsButton volumeButton;
	StsJPanel mainPane = new StsJPanel();
	StsJPanel leftPane = new StsJPanel();
	StsJPanel rightPane = new StsJPanel();
    public StsProgressPanel progressPanel = StsProgressPanel.constructorWithCancelButton();

    public StsProcessVolumePanel(StsVolumeCurvatureWizard wizard, StsProcessVolume wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
        buildPanel();
    }

    public void buildPanel()
    {
        curvatureBox.gbc.fill = GridBagConstraints.HORIZONTAL;

        curvAttributeBean = new StsComboBoxFieldBean(this, "attrString","Attribute:", StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS);
    	curvAttributeBean.setToolTipText("Select Curvature Attribute.");
        curvatureBox.add(curvAttributeBean);

        filterSizeBean = new StsIntFieldBean(this, "filterSize", true, "Filter Size:", true);
        filterSizeBean.setValueAndRangeFixStep(3, 3, 21, 2);
        filterSizeBean.setToolTipText("Specify anaylsis window size.");
        curvatureBox.add(filterSizeBean);

        runAllBean = new StsBooleanFieldBean(this, "runAllPatches", "Run all patches");
        runAllBean.setToolTipText("Run curvature on all instead of selected patches.");
        curvatureBox.add(runAllBean);
        runButton = new StsButton("Run Analysis", "Push button to start analysis.", this, "analyze");
        curvatureBox.add(runButton);

        editRangeBox.gbc.fill = GridBagConstraints.HORIZONTAL;
        curvatureMinBean = new StsFloatFieldBean(wizard, "curvatureMin", true, "Min: ");
        curvatureMaxBean = new StsFloatFieldBean(wizard, "curvatureMax", true, "Max: ");
        editRangeBox.addToRow(curvatureMinBean);
        editRangeBox.addEndRow(curvatureMaxBean);

        outputVolumeBox.gbc.fill = GridBagConstraints.HORIZONTAL;

        interpSizeBean = new StsIntFieldBean(this, "interpSize", true, "Max vertical interpolation (samples):", true);
        interpSizeBean.setValueAndRangeFixStep(3, 0, 100, 1);
        interpSizeBean.setToolTipText("Specify trace Interpolation size.");
        outputVolumeBox.add(interpSizeBean);
        volumeButton = new StsButton("Output Volume", "Push button to create attribute volume.", this, "processVolume");
        outputVolumeBox.add(volumeButton);

        leftPane.gbc.fill = GridBagConstraints.HORIZONTAL;
        leftPane.gbc.anchor = GridBagConstraints.WEST;
        leftPane.addEndRow(curvatureBox);
        leftPane.addEndRow(editRangeBox);
        leftPane.addEndRow(outputVolumeBox);

        mainPane.addToRow(leftPane);

        colorscaleBean = new StsEditableColorscaleFieldBean(StsSurface.class, "colorscale");
        rightPane.add(colorscaleBean);

        mainPane.addEndRow(rightPane);
        addEndRow(mainPane);

        gbc.fill = GridBagConstraints.HORIZONTAL;
        addEndRow(progressPanel);
    }

    public void runUpdateColorscale()
    {
        Runnable runnable = new Runnable()
        {
            public void run()
            {
                updateColorscale();
            }
        };
        StsToolkit.runLaterOnEventThread(runnable);
    }

    public void updateColorscale()
    {
        StsException.systemDebug(this, "updateColorscale", "calculateHistogram(histogramSamples, nHistogramSamples");
        StsPatchVolume patchVolume = wizard.getPatchVolume();
        if(patchVolume == null) return;
        StsColorscale colorscale = patchVolume.getColorscale();
        if(colorscale == null) return;

        colorscaleBean.setHistogram(patchVolume.getHistogram());
        colorscaleBean.setValueObject(colorscale);
        float[] range = colorscale.getRange();
        curvatureMinBean.setValue(range[0]);
        curvatureMaxBean.setValue(range[1]);
        // colorscaleBean.repaint();
    }

//    public void setMessage(String msg)
//    {
//    	msgBean.setText(msg);
//    	msgPanel.setVisible(true);
//    	wizard.rebuild();
//    }
    
    public void setFilterSize(int size) { filterSize = size; }
    public int getFilterSize() { return filterSize; }
    public void setInterpSize(int size) { interpSize = size; }
    public int getInterpSize() { return interpSize; }

    public String getAttrString()
    {
    	return StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS[attributeType];
    }

    public void setAttrString(String string)
    {
        attributeType = StsParameters.getStringMatchByteIndex(StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS, string);
    }
    
    public byte getAttributeType()
    {
    	return attributeType;
    }

    public boolean analyze()
    {
        String attributeMessage = "Computing " + StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS[attributeType];
//        wizardStep.setAnalysisMessage(attributeMessage + " on volume (" + wizard.getSelectedVolume().getName() + ") with a "
//                + filterSize + "x" + filterSize + " filter.");
//        progressPanel.appendLine(attributeMessage + " on volume (" + wizard.getSelectedVolume().getName() + ") with a "
//                + filterSize + "x" + filterSize + " filter.");
        wizard.rebuild();
        Main.logUsage();
        runPatchCurvature();
        return true;
    }

    public boolean saveToModel()
    {
        return true;
    }

    public void runPatchCurvature()
    {
        Runnable runPatchCurvature = new Runnable()
        {
            public void run()
            {
            	StsPatchVolume patchVolume = wizard.getPatchVolume();
            	if (patchVolume == null)
            	{
            		//TODO
            		return;
            	}
                patchVolume.runCurvature(progressPanel, filterSize, attributeType, runAllPatches);
                runUpdateColorscale();

                //wizardStep.enableFinish();
		        wizardStep.enableNext();
            }
        };
        StsToolkit.runRunnable(runPatchCurvature);
//        saveButton.setEnabled(true);

    }
    public void processVolume()
	{
        Runnable runCreateVolume = (new Runnable()
		{
			public void run()
			{
                StsPatchVolume volume = wizard.getPatchVolume();
                StsCurvatureVolumeConstructor.createInterpolatedVolume(wizard.model, volume, progressPanel, attributeType, interpSize);
                // wizard.enableFinish();
                wizard.model.win3dDisplay();
			}
		});

		Thread createVolumeThread = new Thread(runCreateVolume);
		createVolumeThread.start();
	}

    public void setRunAllPatches(boolean runAll)
    {
        this.runAllPatches = runAll;
    }

    public boolean getRunAllPatches()
    {
        return runAllPatches;
    }
}