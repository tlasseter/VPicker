package com.Sts.Actions.Wizards.SurfaceCurvature;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.MVC.View3d.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 *
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsPatchPickPanel extends StsFieldBeanPanel // implements ActionListener
{
    StsVolumeCurvatureWizard wizard;
    private StsPatchPick wizardStep;

    StsSeismicGridEditPanel gridEditPanel = new StsSeismicGridEditPanel();
    //volume patch output info
    // StsGroupBox operationsBox = new StsGroupBox("Volume Patch Setup");
	StsGroupBox filterBox = new StsGroupBox("Patch filter parameters");
	StsGroupBox debugBox = new StsGroupBox("Debug parameters");
    //control parameters
    StsGroupBox pickBox = new StsGroupBox("Volume Patch Parameters");
    StsComboBoxFieldBean typeListBean = new StsComboBoxFieldBean();
    // StsIntFieldBean windowSizeBean = new StsIntFieldBean();
    // StsIntFieldBean testWindowSizeBean = new StsIntFieldBean();
    StsFloatFieldBean maxPickDifBean = new StsFloatFieldBean();
    StsIntFieldBean minPatchSizeBean = new StsIntFieldBean();
    StsIntFieldBean corWindowBean = new StsIntFieldBean();
    // StsFloatFieldBean minCorrelationBean = new StsFloatFieldBean();
    StsFloatFieldBean minAmpFractionBean = new StsFloatFieldBean();
    StsFloatFieldBean maxStretchBean = new StsFloatFieldBean();

    StsGroupBox processBox = new StsGroupBox("AutoPick Process");
    public StsButton runButton;
    StsBooleanFieldBean isIterativeBean = new StsBooleanFieldBean();
    StsFloatFieldBean autoCorMaxBean = new StsFloatFieldBean();
    StsFloatFieldBean autoCorMinBean = new StsFloatFieldBean();
    StsFloatFieldBean autoCorIncBean = new StsFloatFieldBean();
    StsFloatFieldBean minCorrelBean = new StsFloatFieldBean();
    StsFloatFieldBean minAmpRatioBean = new StsFloatFieldBean();

    StsBooleanFieldBean filterBean = new StsBooleanFieldBean();
	StsIntFieldBean boxFilterWidthBean = new StsIntFieldBean();

	StsBooleanFieldBean useFalseTypesBean = new StsBooleanFieldBean();
	StsBooleanFieldBean checkBackMatchBean = new StsBooleanFieldBean();
	// StsBooleanFieldBean cycleSkipOnlyBean = new StsBooleanFieldBean();

	StsBooleanFieldBean debugBean = new StsBooleanFieldBean();
	StsIntFieldBean patchIdBean = new StsIntFieldBean();
	StsIntFieldBean pointRowBean = new StsIntFieldBean();
	StsIntFieldBean pointColBean = new StsIntFieldBean();
	StsIntFieldBean pointSliceBean = new StsIntFieldBean();

    JLabel statusLbl = new JLabel("Pick Status...");
    // private float minCorrelation = defaultMinCorrel;
    public int corWavelength = defaultCorWavelength;
    public float maxPickDif = defaultPickDifWavelength;
    public int minPatchSize = defaultMinPatch;
    public float minAmpFraction = defaultMinAmpFraction;
    // private int windowSize = defaultWindow;
    // private int testWindowSize = defaultWindow;
    public float minAmpRatio = defaultMinAmpRatio;
    public float maxStretch = defaultMaxStretchWavelength;

    /** indicates picking is being run iteratively from autoCorMax to autoCorMin; otherwise run once at minCorrel */
    public boolean isIterative = true;
    /** manual picking minimum acceptable cross-correlation */
    public float manualCorMin = defaultMinCorrel;
    /** iterative picking operation max correl */
    public float autoCorMax = defaultAutoCorMax;
    /** iterative picking operation max correl */
    public float autoCorMin = defaultAutoCorMin;
    /** iterative picking operation max correl */
    public float autoCorInc = defaultAutoCorInc;
	public boolean filter = true;
	public int boxFilterWidth = 1;
	public boolean useFalseTypes = true;
	/** when closing a local box ij box, don't allow same patch at two different points on closing trace */
	public boolean checkBackMatch = true;
	//public boolean cycleSkipOnly = false;
	public boolean debug = false;
	public int patchId = NO_DEBUG;
	public int pointRow = NO_DEBUG;
	public int pointCol = NO_DEBUG;
	public int pointSlice = NO_DEBUG;

    public static final String[] stopCriteriaNames = new String[]{"Stop", "Replace", "Stop if same Z"};

    public byte pickType = StsPatchVolume.PICK_ALL;
    public String pickTypeName = StsPatchVolume.pickTypeNames[StsPatchVolume.PICK_ALL];

    public StsProgressBar progressBar = StsProgressBar.constructorWithInsetsAndCancel();

    static final String autoButtonTip = "Runs all picks from max to min correl.";
    static final String runButtonTip = "After selecting seed, this runs picking.";
    static final String testButtonTip = "Using same seed with new parameters, this reruns picking.";
    static final String deletePatchButtonTip = "Delete the selected seed and patch.";
    static final String deletePicksButtonTip = "Delete the patch for this seed.";
    static final String noCorrelString = "Correlation Range: not available.";
    static final String noPickDifString = "Pick Dif Range: not available.";
    static final String noWaveLengthString = "WaveLength Range: not available.";
	static final String filterString  = "Smooth patches before drawing";
	static final String boxFilterWidthString = "Width of box filter for smoothing (3x3 has width of 1, 5x5 has width of 2).";
	static final String useFalseTypesTip = "Allow false types to correlate with true types (e.g., false Max and Max.";
	static final String checkBackMatchTip = "After initial match, see if it matches back; if not, reject match.";
    static final String minAmpRatioTip = "Amplitude ratio of two points (smaller/larger) must be >= current iterative correlation.";
	//static final String cycleSkipOnlyTip = "DEBUG: only run the cycle skip check.";
    static public final float defaultAutoCorMax = 0.9f;
    static public final float defaultAutoCorMin = 0.5f;
    static public final float defaultAutoCorInc = 0.1f;
    static public final float defaultMinCorrel = 0.8f;
    static public final float defaultMinAmpRatio = 0.4f;

    static public final String displayPropertyNone = "None";
    static public final String displayPropertyPatchColor = "Patch Colors";
    static public final String displayPropertyCorrelCoefs = "Correlation Coefficients";

    // static public final float defaultMinCorrel = 0.80f;
    static public final float defaultPickDifWavelength = 0.25f;
    static public final int defaultCorWavelength = 0;
    static public final int defaultMinPatch = 1;
    static public final float defaultMaxStretchWavelength = 0.25f;
    // static public final int defaultWindow = 21;
    static public final float defaultMinAmpFraction = 0.01f;

	static final int NO_DEBUG = StsPatchGrid.NO_DEBUG;

    public StsPatchPickPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = (StsVolumeCurvatureWizard) wizard;
        this.wizardStep = (StsPatchPick) wizardStep;
        try
        {
            typeListBean.initialize(this, "pickTypeName", "Pick Type:                   ", StsPatchVolume.pickTypeNames);
            // windowSizeBean.initialize(this, "windowSize", true, "Window Size:                 ");
            // testWindowSizeBean.initialize(this, "testWindowSize", true, "Test Window Size:  ");
            maxPickDifBean.initialize(this, "maxPickDif", true, "Max Dif between Trace Picks: ");
            maxPickDifBean.setToolTipText("Max dif allowed between picks (wavelengths)");
            corWindowBean.initialize(this, "corWavelength", true, "Correlation window size wavelengths:     ");
            corWindowBean.setToolTipText("Correlation window size in wavelengths (0 is half-wave).");
            minAmpFractionBean.initialize(this, "minAmpFraction", true, "Min amplitude fraction:   ");
            minAmpFractionBean.setToolTipText("Volume max amp times this factor to be correlated.");
            maxStretchBean.initialize(this, "maxStretch", true, "Max window stretch multiplier:   ");
            maxStretchBean.setToolTipText("Max stretch allowed between correlation windows");
            minPatchSizeBean.initialize(this, "minPatchSize", true, "Min allowable patch size:    ");
            minPatchSizeBean.setToolTipText("Minimum number of points allowed for a patch");
            runButton = new StsButton("Run", runButtonTip, this, "pickVolume");

            isIterativeBean.initialize(this, "isIterative", "Iterative Processing");
            //isIterativeBean.setSelected(true);
            autoCorMaxBean.initialize(this, "autoCorMax", true, "Start Correlation: ");
            autoCorMaxBean.setValueAndRangeFixStep(defaultAutoCorMax, 0.0f, 1.0f, 0.1f);
            autoCorMinBean.initialize(this, "autoCorMin", true, "End Correlation: ");
            autoCorMinBean.setValueAndRangeFixStep(defaultAutoCorMin, 0.0f, 1.0f, 0.1f);
            autoCorIncBean.initialize(this, "autoCorInc", true, "Correlation Increment: ");
            autoCorIncBean.setValueAndRange(defaultAutoCorInc, 0.01f, 0.2f);
            minCorrelBean.initialize(this, "manualCorMin", true, "Min Correlation:  ");
            minCorrelBean.setToolTipText("Minimum correlation criteria.");
            minAmpRatioBean.initialize(this, "minAmpRatio", true, "Min Amp Ratio (< 1):  ");
            minAmpRatioBean.setToolTipText(minAmpRatioTip);
			filterBean.initialize(this, "filter", "Smooth box filter");
			filterBean.setToolTipText(filterString);
			boxFilterWidthBean.initialize(this, "boxFilterWidth", true, "Box filter width");
			boxFilterWidthBean.setToolTipText(boxFilterWidthString);

			useFalseTypesBean.initialize(this, "useFalseTypes", "Use false types");
			useFalseTypesBean.setToolTipText(useFalseTypesTip);
			checkBackMatchBean.initialize(this, "checkBackMatch", "Check back match");
			checkBackMatchBean.setToolTipText(checkBackMatchTip);
			//cycleSkipOnlyBean.initialize(this, "cycleSkipOnly", "Cycle skip check only");
			//cycleSkipOnlyBean.setToolTipText(cycleSkipOnlyTip);

			debugBean.initialize(this, "debug", "Debug");
			patchIdBean.initialize(this, "patchId", true, "Patch ID");
			pointRowBean.initialize(this, "pointRow", true, "Point row");
			pointColBean.initialize(this, "pointCol", true, "Point col");
			pointSliceBean.initialize(this, "pointSlice", true, "Point slice");
            buildPanel();
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    public void initialize()
    {
        try
        {
            StsSeismicVolume seismicVolume = wizard.getSelectedVolume();
            gridEditPanel.initialize(seismicVolume);
        }
        catch (Exception e)
        {
            e.printStackTrace();
        }
    }

    private void buildPanel() throws Exception
    {
        this.gbc.fill = GridBagConstraints.BOTH;
        pickBox.gbc.fill = GridBagConstraints.BOTH;
        pickBox.add(typeListBean);
        // pickBox.add(windowSizeBean);
        // pickBox.add(corWindowBean);
        pickBox.add(minAmpFractionBean);
        // pickBox.add(maxPickDifBean);
//		pickBox.add(maxStretchBean);
        pickBox.add(minPatchSizeBean);
		pickBox.add(useFalseTypesBean);
		pickBox.add(checkBackMatchBean);
		//pickBox.add(cycleSkipOnlyBean);
        StsJPanel operationsPanel = new StsJPanel();
        operationsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
        operationsPanel.add(gridEditPanel);
        operationsPanel.add(pickBox);

        addEndRow(operationsPanel);

        processBox.gbc.fill = GridBagConstraints.BOTH;
        processBox.add(isIterativeBean);
        processBox.add(autoCorMaxBean);
        processBox.add(autoCorMinBean);
        processBox.add(autoCorIncBean);
        processBox.add(minCorrelBean);
        processBox.add(minAmpRatioBean);
		minCorrelBean.setEditable(false);
        processBox.gbc.gridwidth = 2;
        processBox.add(runButton);

        addEndRow(processBox);

		filterBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		filterBox.add(filterBean);
		filterBox.add(boxFilterWidthBean);

		addEndRow(filterBox);

		debugBox.gbc.fill = GridBagConstraints.HORIZONTAL;
		debugBox.add(debugBean);
		debugBox.add(patchIdBean);
		debugBox.add(pointRowBean);
		debugBox.add(pointColBean);
		debugBox.add(pointSliceBean);

		addEndRow(debugBox);

        addEndRow(progressBar);

        wizard.rebuild();
    }

    public void setStatusLabel(String label)
    {
        statusLbl.setText(label);
        repaint();
        return;
    }

    public void togglePickBeans(boolean b)
    {
        typeListBean.setEditable(b);
    }

    public void pickVolume()
    {
//      String attributeMessage = "Computing " + StsSurfaceCurvatureAttribute.CURV_ATTRIBUTE_STRINGS[attributeType];
//      wizardStep.setAnalysisMessage(attributeMessage + " on volume (" + wizard.getSelectedVolume().getName() + ") with a "
//              + filterSize + "x" + filterSize + " filter.");
        Main.logUsage();

        runPickVolume();
    }

/*
  public void testParameters()
  {
      StsPatchVolume patchVolume = wizard.getPatchVolume();
      if(patchVolume.getTestParameters(windowSize))
          setRangeLabels();
      else
          setRangeLabelsNull();


  }
*/

    public void runPickVolume()
    {
        Runnable runPickVolume = new Runnable()
        {
            public void run()
            {
                StsPatchVolume patchVolume = wizard.getPatchVolume();
                patchVolume.setCroppedBoundingBox(gridEditPanel.croppedBoundingBox);
                patchVolume.constructPatchVolume(StsPatchPickPanel.this);
                wizardStep.enableNext();
                wizardStep.enableFinish();
            }
        };
        StsToolkit.runRunnable(runPickVolume);
    }

    /*
        public void computeTestCorrelRange()
        {
            StsPatchVolume patchVolume = wizard.getPatchVolume();
            if(patchVolume.computeTestCorrelRange())
                setRangeLabels();
            else
                setRangeLabelsNull();

        }

        private void setRangeLabels()
        {
            StsPatchVolume patchVolume= wizard.getPatchVolume();
            corRangeLabel.setText("Correlation Range: " +  patchVolume.correlTestMin + " to " + patchVolume.correlTestMax);
            pickDifRangeLabel.setText("Pick Dif Range.  Avg: " +  patchVolume.avgTestPickDif + " Max: " + patchVolume.maxTestPickDif);
            waveLengthRangeLabel.setText("Wavelength Range.  Min: " +  patchVolume.minWaveLength + " Max: " + patchVolume.maxWaveLength);

        }

        private void setRangeLabelsNull()
        {
            corRangeLabel.setText(noCorrelString);
            pickDifRangeLabel.setText(noPickDifString);
            waveLengthRangeLabel.setText(noWaveLengthString);
        }
    */
    public byte getPickType()
    { return pickType; }

    /*
        public float getMinCorrelation()
        {
            return minCorrelation;
        }

        public void setMinCorrelation(float corMin)
        {
            this.minCorrelation = corMin;
        }
    */
    public int getCorWavelength()
    {
        return corWavelength;
    }

    public void setMinAmpRatio(float value)
    {
        minAmpRatio = value;
    }

    public float getMinAmpRatio()
    {
        return minAmpRatio;
    }

    public void setCorWavelength(int cor)
    {
        this.corWavelength = cor;
    }

    public void setMinPatchSize(int patchSize)
    {
        this.minPatchSize = patchSize;
    }

    public int getMinPatchSize() { return minPatchSize; }

    /*
        public void setWindowSize(int size)
        {
            this.windowSize = size;
        }

        public int getWindowSize() { return windowSize; }

        public void setTestWindowSize(int size)
        {
            this.testWindowSize = size;
        }

        public int getTestWindowSize() { return testWindowSize; }
    */
    public void setMaxPickDif(float maxPickDif)
    {
        this.maxPickDif = maxPickDif;
    }

    public float getMaxPickDif()
    {
        return maxPickDif;
    }

    public void setMaxStretch(float maxStretch)
    {
        this.maxStretch = maxStretch;
    }

    public float getMaxStretch()
    {
        return maxStretch;
    }

    public String getPickTypeName()
    {
        return StsPatchVolume.pickTypeNames[pickType];

    }

    public void setPickTypeName(String string)
    {
        pickType = StsParameters.getStringMatchByteIndex(StsPatchVolume.pickTypeNames, string);
    }

//    public StsPatchVolume getPatchVolume()
//    {
//          return patchVolume;
//    }
//    
//    public void setPatchVolume(StsPatchVolume volume)
//    {
//          patchVolume = volume;
//    }

    static public void main(String[] args)
    {
        StsModel model = new StsModel();
        StsActionManager actionManager = new StsActionManager(model);
        model.mainWindowActionManager = actionManager;
        StsVolumeCurvatureWizard wizard = new StsVolumeCurvatureWizard(actionManager);
        StsPatchPickPanel panel = new StsPatchPickPanel(wizard, new StsPatchPick(wizard));
        StsToolkit.createDialog(panel, true, 200, 400);
    }

    public float getAutoCorMax()
    {
        return autoCorMax;
    }

    public void setAutoCorMax(float autoCorMax)
    {
        this.autoCorMax = autoCorMax;
    }

    public float getAutoCorMin()
    {
        return autoCorMin;
    }

    public void setAutoCorMin(float autoCorMin)
    {
        this.autoCorMin = autoCorMin;
    }

    public float getAutoCorInc()
    {
        return autoCorInc;
    }

    public void setAutoCorInc(float autoCorInc)
    {
        this.autoCorInc = autoCorInc;
    }

    public void setIsIterative(boolean isIterative)
    {
        if (this.isIterative == isIterative) return;
        this.isIterative = isIterative;

        if (!isIterative)
        {
            autoCorMaxBean.setEditable(false);
            autoCorMinBean.setEditable(false);
            autoCorIncBean.setEditable(false);
            minCorrelBean.setEditable(true);
        }
        else
        {
            autoCorMaxBean.setEditable(true);
            autoCorMinBean.setEditable(true);
            autoCorIncBean.setEditable(true);
            minCorrelBean.setEditable(false);
        }
    }

    public boolean getIsIterative() { return isIterative; }

    public float getManualCorMin()
    {
        return manualCorMin;
    }

    public void setManualCorMin(float manualCorMin)
    {
        this.manualCorMin = manualCorMin;
    }

	/** if true, use falseTypes (false max, min, +/- zero crossings */
	public boolean getUseFalseTypes()
	{
		return useFalseTypes;
	}

	public void setUseFalseTypes(boolean useFalseTypes)
	{
		this.useFalseTypes = useFalseTypes;
	}

	/** when closing a local box ij box, don't allow same patch at two different points on closing trace */
	public boolean getCheckBackMatch()
	{
		return checkBackMatch;
	}

	public void setCheckBackMatch(boolean checkBackMatch)
	{
		this.checkBackMatch = checkBackMatch;
	}

	public boolean getDebug()
	{
		return debug;
	}

	public void setDebug(boolean debug)
	{
		this.debug = debug;
	}

	public int getPatchId()
	{
		return patchId;
	}

	public void setPatchId(int patchId)
	{
		this.patchId = patchId;
	}

	public int getPointRow()
	{
		return pointRow;
	}

	public void setPointRow(int pointRow)
	{
		this.pointRow = pointRow;
	}

	public int getPointCol()
	{
		return pointCol;
	}

	public void setPointCol(int pointCol)
	{
		this.pointCol = pointCol;
	}

	public int getPointSlice()
	{
		return pointSlice;
	}

	public void setPointSlice(int pointSlice)
	{
		this.pointSlice = pointSlice;
	}

	/** smooth patch surfaces with box filter */
	public boolean getFilter()
	{
		return filter;
	}

	public void setFilter(boolean filter)
	{
		this.filter = filter;
		wizard.getPatchVolume().setFilter(filter);
		wizard.model.repaintViews(StsView3d.class);
	}

	/** box filter width (3x3 filter has width of 1, 5x5 has width of 2 */
	public int getBoxFilterWidth()
	{
		return boxFilterWidth;
	}

	public void setBoxFilterWidth(int boxFilterWidth)
	{
		this.boxFilterWidth = boxFilterWidth;
		wizard.getPatchVolume().setBoxFilterWidth(boxFilterWidth);
		wizard.model.repaintViews(StsView3d.class);
	}

    public float getMinAmpFraction()
    {
        return minAmpFraction;
    }

    public void setMinAmpFraction(float minAmpFraction)
    {
        this.minAmpFraction = minAmpFraction;
    }

    /** For debugging: only run the cycle skip check */
	//public boolean isCycleSkipOnly()
	//{
	//	return cycleSkipOnly;
	//}

	//public void setCycleSkipOnly(boolean cycleSkipOnly)
	//{
	//	this.cycleSkipOnly = cycleSkipOnly;
	//}
}
