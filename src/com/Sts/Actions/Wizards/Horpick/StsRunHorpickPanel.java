package com.Sts.Actions.Wizards.Horpick;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author TJLasseter
 * @version beta 1.0
 */

public class StsRunHorpickPanel extends StsFieldBeanPanel // implements ActionListener
{
    StsHorpickWizard wizard;

    StsGroupBox traceBox = new StsGroupBox("Selected Trace");

    StsGroupBox operationsBox = new StsGroupBox("Seed Operations");
//    public StsToggleButton selectButton;
//    public StsButton runButton;
    public StsButton deletePatchButton;
    public StsButton deletePicksButton;

    StsGroupBox processBox = new StsGroupBox("AutoPick Process");
    public StsButton runButton;
    StsBooleanFieldBean isIterativeBean = new StsBooleanFieldBean();
//    StsCheckbox iterativeChk = new StsCheckbox();
    JLabel statusLbl = new JLabel("Pick Status...");
    StsFloatFieldBean autoCorMaxBean = new StsFloatFieldBean();
    StsFloatFieldBean autoCorMinBean = new StsFloatFieldBean();
    StsFloatFieldBean autoCorIncBean = new StsFloatFieldBean();
    StsFloatFieldBean minCorrelBean = new StsFloatFieldBean();

    StsGroupBox pickBox = new StsGroupBox("Seed Parameters");
    StsGroupBox statusBox = new StsGroupBox("Pick Status");
    StsColorItemComboBoxFieldBean patchColorItemList = new StsColorItemComboBoxFieldBean();
    StsFloatFieldBean pickZBean = new StsFloatFieldBean();
    StsFloatFieldBean inlineBean = new StsFloatFieldBean();
    StsFloatFieldBean xlineBean = new StsFloatFieldBean();
    StsComboBoxFieldBean typeListBean = new StsComboBoxFieldBean();
    StsIntFieldBean windowSizeBean = new StsIntFieldBean();
    StsFloatFieldBean maxPickDifBean = new StsFloatFieldBean();
    JLabel ilinePickLabel = new JLabel(noILineString);
    JLabel xlinePickLabel = new JLabel(noXLineString);
    JLabel zPickLabel = new JLabel(noZString);
    JLabel corRangeLabel = new JLabel(noCorrelString);
    JLabel pickDifRangeLabel = new JLabel(noPickDifString);
    JLabel waveLengthRangeLabel = new JLabel(noWaveLengthString);

    static final String autoButtonTip = "Runs all picks from max to min correl.";
    static final String runButtonTip = "After selecting seed, this runs picking.";
    static final String rerunButtonTip = "Using same seed with new parameters, this reruns picking.";
    static final String deletePatchButtonTip = "Delete the selected seed and patch.";
    static final String deletePicksButtonTip = "Delete the patch for this seed.";
    static final String noCorrelString = "Correlation Range: not available.";
    static final String noPickDifString = "Pick Dif Range: not available.";
    static final String noWaveLengthString = "WaveLength Range: not available.";
    static final String noILineString = "Iline: not available.";
    static final String noXLineString = "Xline: not available.";
    static final String noZString = "Time or depth: not available.";
    public static final byte PICK_MAX = 0;
    public static final byte PICK_MIN = 1;
    public static final byte PICK_ZERO_PLUS = 2;
    public static final byte PICK_ZERO_MINUS = 3;

    public static final String[] pickTypeNames = new String[] { "Maximum", "Minimum", "Zero-crossing, up.", "Zero-crosing, down." };
    public static final String[] stopCriteriaNames = new String[] { "Stop", "Replace", "Stop if same Z" };
    public static final StsColorListItem nullColorListItem = StsColorListItem.nullColorListItem();

    protected byte pickType = PICK_MAX;

    public StsRunHorpickPanel(StsWizard wizard)
    {
        this.wizard = (StsHorpickWizard)wizard;
		try
		{
//			selectButton = new StsToggleButton("Select Seed", "Select seed point on vertical seismic cursor in 2d or 3d.", this, "setSelectSeed");
			deletePatchButton = new StsButton("Delete Seed & Patch", deletePatchButtonTip, wizard, "deletePatch");
			deletePicksButton = new StsButton("Delete Patch", deletePicksButtonTip, wizard, "deletePicks");

			isIterativeBean.classInitialize(StsHorpick.class, "isIterative", "Iterative Processing");
			isIterativeBean.setSelected(false);
			autoCorMaxBean.classInitialize(StsHorpick.class, "autoCorMax", true, "Start Correlation: ");
			autoCorMaxBean.setValueAndRangeFixStep(StsHorpick.defaultAutoCorMax, 0.0f, 1.0f, 0.01f);
			autoCorMinBean.classInitialize(StsHorpick.class, "autoCorMin", true, "End Correlation: ");
			autoCorMinBean.setValueAndRangeFixStep(StsHorpick.defaultAutoCorMin, 0.0f, 1.0f, 0.01f);
			autoCorIncBean.classInitialize(StsHorpick.class, "autoCorInc", true, "Correlation Increment: ");
			autoCorIncBean.setValueAndRange(StsHorpick.defaultAutoCorInc, 0.001f, 0.1f);
			minCorrelBean.classInitialize(StsHorpick.class, "manualCorMin", true, "Min Correlation:  ");
//            minCorrelBean.setValueAndRange(StsHorpick.defaultMinCorrel, 0.0f, 1.0f);
			minCorrelBean.setToolTipText("Minimum correlation criteria.");

//			Color[] patchColors = wizard.getModel().getSpectrum("Basic").getColors();
			patchColorItemList.initialize(wizard, "pickPatchColorItem", "Patch:", new StsColorListItem[] { nullColorListItem } );
			pickZBean.initialize(wizard, "pickZ", true, "T or Z:  ");
			inlineBean.initialize(wizard, "pickInline", true, "Inline:  ");
			xlineBean.initialize(wizard, "pickXline", true, "Crossline:  ");
			typeListBean.initialize(this, "pickTypeName", "Pick Type:  ", pickTypeNames);
			windowSizeBean.initialize(wizard, "windowSize", true, "Window Size: ");
			maxPickDifBean.initialize(wizard, "maxPickDif", true, "Max Dif between Trace Picks: ");
			maxPickDifBean.setValueAndRange(1.0f, 0.0f, 10.0f);

			runButton = new StsButton("Run", runButtonTip, this, "runPickSurface");

			jbInit();
		}
		catch(Exception e)
		{
			e.printStackTrace();
        }
    }

	public void initialize()
	{
        try
        {
//			selectButton.setEnabled(true);
//			selectButton.setSelected(false);
			StsHorpick horpick = wizard.getHorpick();
			isIterativeBean.setBeanObject(horpick);
			autoCorMaxBean.setBeanObject(horpick);
			autoCorMinBean.setBeanObject(horpick);
			autoCorIncBean.setBeanObject(horpick);
			minCorrelBean.setBeanObject(horpick);
//			minCorrelBean.setValueFromPanelObject(horpick);
			patchColorItemList.setBeanObject(wizard);
			windowSizeBean.setBeanObject(wizard);
			maxPickDifBean.setBeanObject(wizard);
			int nItems = initializeComboBox();
			deletePatchButton.setEnabled(nItems > 0);
			deletePicksButton.setEnabled(nItems > 0);
			setIsIterative(horpick.getIsIterative());
        }
        catch(Exception e)
        {
            e.printStackTrace();
        }
    }

    private void jbInit() throws Exception
    {
		this.gbc.fill = GridBagConstraints.BOTH;

	    operationsBox.gbc.fill = GridBagConstraints.BOTH;
//		operationsBox.addToRow(selectButton);
		operationsBox.addToRow(deletePatchButton);
		operationsBox.addEndRow(deletePicksButton);
		operationsBox.gbc.gridwidth = 3;
        operationsBox.add(ilinePickLabel);
        operationsBox.add(xlinePickLabel);
        operationsBox.add(zPickLabel);
        operationsBox.add(corRangeLabel);
		operationsBox.add(pickDifRangeLabel);
		operationsBox.add(waveLengthRangeLabel);

	    processBox.gbc.fill = GridBagConstraints.BOTH;
	    processBox.add(isIterativeBean);
		processBox.add(autoCorMaxBean);
		processBox.add(autoCorMinBean);
		processBox.add(autoCorIncBean);
		processBox.add(minCorrelBean);
		processBox.gbc.gridwidth = 2;
		processBox.add(runButton);

	    pickBox.gbc.fill = GridBagConstraints.BOTH;
		pickBox.add(patchColorItemList);
		pickBox.add(pickZBean);
		pickBox.add(inlineBean);
		pickBox.add(xlineBean);
		pickBox.add(typeListBean);
		pickBox.add(windowSizeBean);
		pickBox.add(maxPickDifBean);

		statusBox.gbc.fill = GridBagConstraints.BOTH;
        statusBox.add(statusLbl);

		StsJPanel operationsPanel = new StsJPanel();
		operationsPanel.gbc.fill = GridBagConstraints.HORIZONTAL;
		operationsPanel.add(operationsBox);
		operationsPanel.add(pickBox);
		operationsPanel.add(processBox);
		operationsPanel.add(statusBox);

	    addToRow(traceBox);
	    addEndRow(operationsPanel);

		wizard.rebuild();
    }

    public void setIsIterative(boolean isIterative)
    {
        if(!isIterative)
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

    public void setStatusLabel(String label_)
    {
        final String label = label_;
        StsToolkit.runLaterOnEventThread
        (
            new Runnable()
            {
                public void run()
                {
                    statusLbl.setText(label);
                    repaint();
                }
            }
        );
    }

    public void togglePickBeans(boolean b)
    {
        pickZBean.setEditable(b);
        inlineBean.setEditable(b);
        xlineBean.setEditable(b);
        typeListBean.setEditable(b);
        windowSizeBean.setEditable(b);
        maxPickDifBean.setEditable(b);
        minCorrelBean.setEditable(b);
    }

    public void runPickSurface()
    {
//        selectSeedButton.setSelected(false);
        wizard.runPickSurface();
    }
/*
    public void setSelectSeed()
    {
        boolean selected = selectSeedButton.isSelected();
//        runButton.setEnabled(!selected);
        wizard.setSelectSeed(selected);
    }
*/
/*
    public void selectSeed()
    {
        boolean selected = selectButton.isSelected();
        wizard.setSelectSeed(selected);
//        runButton.setEnabled(!selected);
    }
*/
    public void seedSelected(StsPickPatch patch)
    {
        computeTestCorrelRange();
        addPatchToComboBox(patch);
//        selectButton.setSelected(false);
        deletePatchButton.setEnabled(true);
        deletePicksButton.setEnabled(true);
        setPatchSelected(patch);
        runButton.setEnabled(true);
        typeListBean.setBeanObject(this);
    }

    private int initializeComboBox()
    {
        int nPatches = 0;
        patchColorItemList.removeAllItems();
        StsObjectRefList patches = wizard.getHorpick().getPatches();
        if(patches != null)  nPatches = patches.getSize();
        StsColorListItem listItem = null;
        if(nPatches == 0)
        {
            patchColorItemList.addItem(nullColorListItem);
        }
        else
        {
            listItem = new StsColorListItem(null, StsColor.WHITE, "All", 16, 16, 0);
            patchColorItemList.addItem(listItem);
            for(int n = 0; n < nPatches; n++)
            {
                StsPickPatch patch = (StsPickPatch)patches.getElement(n);
                listItem = new StsColorListItem(patch);
                patchColorItemList.addItem(listItem);
            }
        }
        patchColorItemList.setEditable(true);
        return nPatches;
    }

    public void addPatchToComboBox(StsPickPatch patch)
    {
        // first item is a dummy and has a null object; if this is only item,
        // change it to a legitimate item
        int nItems = patchColorItemList.getListSize();
        StsColorListItem listItem;
    /*
        if(nItems == 1)
        {
            listItem = patchColorList.getItemAt(0);
            if(listItem.getObject() == null) patchColorList.removeItem(listItem);
        }
    */
        if(nItems == 1)
        {
            listItem = patchColorItemList.getItemAt(0);
            listItem.setName("All");
        }
        listItem = new StsColorListItem(patch);
        patchColorItemList.addItem(listItem);
        patchColorItemList.setSelectedItem(listItem);
    }

    public void computeTestCorrelRange()
    {
//        wizard.setRerun(false);
        StsPickPatch pickPatch = wizard.getPickPatch();
        if(pickPatch == null) return;

        if(wizard.computeTestCorrelRange())
            setRangeLabels(pickPatch);
        else
            setRangeLabelsNull();

        pickZBean.getValueFromPanelObject();
    }

    private void setRangeLabels(StsPickPatch pickPatch)
    {
        StsSeismicVolume volume = wizard.seismicVolume;
        StsGridPoint seedPoint = pickPatch.getSeedPoint();
        ilinePickLabel.setText("Iline: " + volume.getNearestRowNumFromY(seedPoint.getY()));
        xlinePickLabel.setText("Xline: " + volume.getNearestColNumFromX(seedPoint.getX()));
        zPickLabel.setText("Time or depth: " + (seedPoint.getZ()));
        corRangeLabel.setText("Correlation Range: " +  pickPatch.correlTestMin + " to " + pickPatch.correlTestMax);
        pickDifRangeLabel.setText("Pick Dif Range.  Avg: " +  pickPatch.avgTestPickDif + " Max: " + pickPatch.maxTestPickDif);
        waveLengthRangeLabel.setText("Wavelength Range.  Min: " +  pickPatch.minWaveLength + " Max: " + pickPatch.maxWaveLength);

    }

    private void setRangeLabelsNull()
    {
        corRangeLabel.setText(noCorrelString);
        pickDifRangeLabel.setText(noPickDifString);
        waveLengthRangeLabel.setText(noWaveLengthString);
    }

    // a generic method should be created in StsFieldBeanPanel
    public void setPatchSelected(StsPickPatch patch)
    {
        if(patch == null)
        {
            pickZBean.setValue(1.e30);
            inlineBean.setValue(0.0);
            xlineBean.setValue(0.0);
            windowSizeBean.setValue(0);
            maxPickDifBean.setValue(1.0);
            setRangeLabelsNull();
        }
        else
        {
            pickZBean.getValueFromPanelObject();
            inlineBean.getValueFromPanelObject();
            xlineBean.getValueFromPanelObject();
            windowSizeBean.getValueFromPanelObject();
            maxPickDifBean.getValueFromPanelObject();
            setRangeLabels(patch);
        }
    }

    public void pickCompleted()
    {
//        selectSeedButton.setSelected(false);
//        selectSeedButton.setEnabled(true);
//        wizard.setRerun(true);
    }

    public void deletePick(StsPickPatch pickPatch)
    {
        StsColorListItem colorListItem;

        if(pickPatch == null) return;
        patchColorItemList.deleteObject(pickPatch);
//        wizard.getHorpick().deletePatch(pickPatch);
        colorListItem = (StsColorListItem) patchColorItemList.getSelectedItem();
        StsPickPatch patch = (StsPickPatch)colorListItem.getObject();
        setPatchSelected(patch);

        int nPatches = wizard.getHorpick().getPatches().getSize();
        if(nPatches == 0)
        {
            colorListItem = patchColorItemList.getItemAt(0);
            colorListItem.setName("None");
            runButton.setEnabled(false);
        }
//        wizard.setRerun(true);
    }

    public void deleteAllPicks()
    {
        int nItems = patchColorItemList.getListSize();
        for(int n = nItems-1; n >= 1; n--)
            patchColorItemList.removeItemAtIndex(n);
        patchColorItemList.getItemAt(0).setName("None");

        wizard.getHorpick().deleteAllPatches();
        runButton.setEnabled(false);
    }


    public void setPickType(byte pickType)
    {
        if(this.pickType == pickType) return;
        this.pickType = pickType;
        wizard.setPickType(pickType);
        computeTestCorrelRange();
    }

    public byte getPickType() { return wizard.getPickType(); }

    public void setPickTypeName(String pickTypeName)
    {
        for(byte n = 0; n < 4; n++)
        {
            if(pickTypeName.equals(pickTypeNames[n]))
            {
                setPickType(n);
                wizard.getModel().win3dDisplay();
                break;
            }
        }
    }
    public String getPickTypeName() { return pickTypeNames[pickType]; }

	static public void main(String[] args)
	{
		StsModel model = new StsModel();
		StsActionManager actionManager = new StsActionManager(model);
        model.mainWindowActionManager = actionManager;
        StsHorpickWizard wizard = new StsHorpickWizard(actionManager);
		StsRunHorpickPanel panel = new StsRunHorpickPanel(wizard);
		StsToolkit.createDialog(panel, true, 200, 400);
	}
}
