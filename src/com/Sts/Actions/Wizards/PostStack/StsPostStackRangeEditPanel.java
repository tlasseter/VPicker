package com.Sts.Actions.Wizards.PostStack;

import com.Sts.Actions.Wizards.Seismic.*;
import com.Sts.Actions.Wizards.*;
import com.Sts.Interfaces.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.Histogram.*;
import com.Sts.UI.Progress.*;
import com.Sts.UI.*;
import com.Sts.UI.Table.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;

/**
 * <p>Title: S2S Development</p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: S2S Systems LLC</p>
 * @author T.Lasseter
 * @version 1.1
 */

public class StsPostStackRangeEditPanel extends StsFieldBeanPanel implements StsSelectRowNotifyListener, StsTableModelListener
{
    /** the volumes selected for this pass thru the wizard */
    StsSeismicBoundingBox[] selectedVolumes;
    // total range of histogram
    float totalDataMin, totalDataMax;
    // current range set by data min and max textboxes
    float dataMin, dataMax;
    private String clipString = StsSeismicBoundingBox.CLIP_NONE_STRING;
    private boolean hasUserNull = false;
    private float userNull = 0.0f;

	public StsSeismic3dRangeEditPanel cropped3dRangeBox = new StsSeismic3dRangeEditPanel();
    private StsGroupBox dataRangeGroupBox = new StsGroupBox("Data Range");
    StsSeismicVolumesHistogramPanel histogramPanel;
    private StsJPanel nullPanel = StsJPanel.addInsets();
    private StsJPanel userNullPanel = StsJPanel.addInsets();
	private StsSegyInformationPanel infoPanel = null;

    public StsProgressPanel progressPanel = StsProgressPanel.constructor(5, 50);
    private StsComboBoxFieldBean clipBean;
    private StsBooleanFieldBean hasUserNullBean;
    private StsFloatFieldBean userNullBean;
    
    private StsSeismicWizard wizard;
	private StsWizardStep wizardStep;

	JButton displaySeismicBtn   = new JButton("Compute Histogram");

    private StsTablePanelNew volumeStatusPanel = null;

	public StsPostStackRangeEditPanel(StsSeismicWizard wizard, StsWizardStep wizardStep)
	{
		this.wizard = wizard;
		this.wizardStep = wizardStep;

        try
		{
            buildTablePanel();
            constructBeans();
            buildPanel();
        }
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	public void initialize()
	{
        cropped3dRangeBox.initialize(wizard.getSegyVolumes());
        cropped3dRangeBox.initializeBeanValuesAndRanges();
        selectedVolumes = wizard.getSelectedVolumes();
        initializeHistogram( wizard.getSegyVolumes());
        volumeStatusPanel.replaceRows(wizard.getSegyVolumesList());
        volumeStatusPanel.setSelectionIndex(0);
        initializeData();

        infoPanel = wizard.getInfoPanel();
        if(infoPanel != null)
        {
            infoPanel.initialize(wizard.getFirstSelectedVolume().getSegyFormat(), wizard);
            infoPanel.addToPanel(this);
            wizard.rebuild();
        }
    }

    public void buildTablePanel()
    {
        String[] columnNames  = {"stemname", "dataMin", "dataMax", "statusString"};
        String[] columnTitles = {"Name", "Data Min", "Data Max", "Status"};
        StsEditableTableModel tableModel = new StsEditableTableModel(StsSeismicBoundingBox.class, columnNames, columnTitles, true);
        volumeStatusPanel = new StsTablePanelNew(tableModel);
        volumeStatusPanel.setLabel("Volumes");
        volumeStatusPanel.setSelectionMode( ListSelectionModel.SINGLE_SELECTION);
        volumeStatusPanel.setSize(400, 50);
        volumeStatusPanel.initialize();
    }

    private void constructBeans()
    {
        clipBean  = new StsComboBoxFieldBean(this, "clipString", "Clip:", StsSeismicBoundingBox.CLIP_STRINGS);
        hasUserNullBean = new StsBooleanFieldBean(this, "hasUserNull", false, null);
        userNullBean = new StsFloatFieldBean(this, "userNull", true, "Assign null value:");
    }

    private void buildPanel()
	{
        gbc.fill = GridBagConstraints.HORIZONTAL;
        gbc.anchor = GridBagConstraints.NORTH;

        add(cropped3dRangeBox);

        dataRangeGroupBox.gbc.fill = GridBagConstraints.BOTH;
        dataRangeGroupBox.add(volumeStatusPanel);

        dataRangeGroupBox.gbc.fill = gbc.BOTH;

        histogramPanel = new StsSeismicVolumesHistogramPanel(volumeStatusPanel);
        dataRangeGroupBox.add(histogramPanel);

        userNullPanel.gbc.fill = gbc.HORIZONTAL;
        userNullPanel.addToRow(hasUserNullBean);
        userNullPanel.addEndRow(userNullBean);

        nullPanel.addToRow(clipBean);
        nullPanel.addToRow(userNullPanel);

        dataRangeGroupBox.add(nullPanel);

        gbc.fill = GridBagConstraints.BOTH;
        add(dataRangeGroupBox);

        add(progressPanel);

        if( !wizard.cropEnabled())
            cropped3dRangeBox.disableCropOptions();

        volumeStatusPanel.addSelectRowNotifyListener(this);
        volumeStatusPanel.addTableModelListener(this);
    }

    public StsSeismicBoundingBox[] getSegyVolumes()
    {
        return wizard.getSegyVolumes();
    }

    private void initializeHistogram(StsSeismicBoundingBox[] volumes)
    {
        totalDataMin = volumes[0].dataMin;
        totalDataMax = volumes[0].dataMax;
        for(int n = 1; n < volumes.length; n++)
        {
            totalDataMin = Math.min(totalDataMin, volumes[n].dataMin);
            totalDataMax = Math.max(totalDataMax, volumes[n].dataMax);
        }
        histogramPanel.setDataRange(totalDataMin, totalDataMax);
//        reanalyzeTraces();
    }

    public float getDataMin()
    {
        return dataMin;
    }

    public void setDataMin( float dataMin)
    {
        this.dataMin = dataMin;
        if( selectedVolumes == null) return;
        for(int n = 0; n < selectedVolumes.length; n++)
            selectedVolumes[n].setDataMin(dataMin);
        volumeStatusPanel.replaceRows(wizard.getSegyVolumesList());
        histogramPanel.setClipDataMin(dataMin);
    }

    public float getDataMax()
    {
        return dataMax;
    }

    public void setDataMax(float dataMax)
    {
        this.dataMax = dataMax;
        if( selectedVolumes == null) return;
        for(int n = 0; n < selectedVolumes.length; n++)
            selectedVolumes[n].setDataMax(dataMax);
        volumeStatusPanel.replaceRows(wizard.getSegyVolumesList());
        histogramPanel.setClipDataMax(dataMax);
    }

     public void initializeData()
     {
         if(selectedVolumes == null) return;
         StsSeismicBoundingBox selectedVolume = selectedVolumes[0];
         dataMin = selectedVolumes[0].dataMin;
         dataMax = selectedVolumes[0].dataMax;
         clipString = selectedVolume.getClipString();
         clipBean.doSetValueObject(clipString);
         checkSetClipEditable();
         hasUserNull =  selectedVolume.getHasUserNull();
         hasUserNullBean.setValue(hasUserNull);
         if(hasUserNull)
         {
             userNull = selectedVolume.getUserNull();
             userNullBean.setValue(userNull);
             userNullBean.setEditable(true);
         }
         else
         {
            userNullBean.setValue(0.0);
            userNullBean.setEditable(false);
         }
         histogramPanel.updateData(selectedVolume);
     }

    public void rowsSelected( int[] selectedIndices)
    {
        wizard.setSelectedVolumes( selectedIndices);
        selectedVolumes = wizard.getSelectedVolumes();
        histogramPanel.setSelectedVolumes(selectedVolumes);
        initializeData(); 
    }

    public void removeRows( int firstRow, int lastRow)
    {
        for( int i = firstRow; i <= lastRow; i++)
            wizard.moveSegyVolumeToAvailableList( wizard.getSegyVolume(i));

        if (wizard.getSegyVolumesList().size() >  0)
            wizard.enableNext();
        else
            wizard.disableNext();
    }

    public void setClipString(String clipString)
    {
        this.clipString = clipString;
        for(int n = 0; n < selectedVolumes.length; n++)
            selectedVolumes[n].setClipString(clipString);
        checkSetClipEditable();
    }

    private void checkSetClipEditable()
    {
        boolean clip = clipString != StsSeismicBoundingBox.CLIP_NONE_STRING;
        histogramPanel.setClipEditable(clip);
    }

    public String getClipString()
    {
        return clipString;
    }

    public void setHasUserNull(boolean value)
    {
        if(hasUserNull == value) return;
        hasUserNull = value;
        userNullBean.setEnabled(hasUserNull);
        selectedVolumes[0].setHasUserNull(hasUserNull);
        if(hasUserNull)
            userNullBean.setValue(selectedVolumes[0].getUserNull());
        else
            userNullBean.setValue(0.0);
        reanalyzeTraces();
    }

    public boolean getHasUserNull()
    {
        return hasUserNull;
    }

    public void setUserNull(float value)
    {
        userNull = value;
        selectedVolumes[0].setUserNull(value);
        reanalyzeTraces();
    }

    public float getUserNull()
    {
        return userNull;
    }

    public void reanalyzeTraces()
    {

        StsProgressRunnable progressRunnable = new StsProgressRunnable()
        {
             StsPostStackAnalyzer analyzer = wizard.getAnalyzer(progressPanel, volumeStatusPanel);
             public void run()
             {
                for(int n = 0; n < selectedVolumes.length; n++)
                {
                    if(!analyzer.analyzeTraces(selectedVolumes[n])) continue;
                    analyzer.analyzeGrid(progressPanel, false, selectedVolumes[n]);
                }
                updatePanel();
             }

             public void cancel()
             {
                 progressPanel.cancel();
             }
        };
        StsToolkit.runRunnable(progressRunnable);
    }

    private void updatePanel()
    {
        initializeData();
    }
}
