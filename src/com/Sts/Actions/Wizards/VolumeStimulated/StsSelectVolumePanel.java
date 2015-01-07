package com.Sts.Actions.Wizards.VolumeStimulated;

import com.Sts.Actions.Wizards.StsWizard;
import com.Sts.Actions.Wizards.StsWizardStep;
import com.Sts.Actions.Wizards.VirtualVolume.StsVirtualVolumeWizard;
import com.Sts.DBTypes.*;
import com.Sts.MVC.StsModel;
import com.Sts.Types.StsArrayList;
import com.Sts.UI.Beans.StsComboBoxFieldBean;
import com.Sts.UI.Beans.StsGroupBox;
import com.Sts.UI.Beans.StsJPanel;
import com.Sts.UI.Beans.StsStringFieldBean;
import com.Sts.UI.Progress.StsProgressPanel;
import com.Sts.UI.StsButton;
import com.Sts.UI.StsCheckbox;

import javax.swing.*;
import java.awt.*;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.TreeSet;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsSelectVolumePanel extends StsJPanel
{
    private StsVolumeStimulatedWizard wizard;
    private StsSelectVolume wizardStep;

    private StsModel model = null;

    private StsSensorVirtualVolume[] virtualVolumes = null;

    private StsSeismicVolume selectedVolume = null;

    StsComboBoxFieldBean volumeCombo = new StsComboBoxFieldBean();

    StsGroupBox analyzeBox = new StsGroupBox("Analyze Volume");
	StsComboBoxFieldBean unitsBean = new StsComboBoxFieldBean();
	StsStringFieldBean msgBean = new StsStringFieldBean(false);
	StsButton runButton;

    StsProgressPanel progressPanel;

    public StsSelectVolumePanel(StsWizard wizard, StsWizardStep wizardStep)
    {
    	this.wizard = (StsVolumeStimulatedWizard)wizard;
    	this.wizardStep = (StsSelectVolume)wizardStep;
        try
        {
            jbInit();
        }
        catch(Exception ex)
        {
            ex.printStackTrace();
        }
    }

    public void initialize()
    {
        model = wizard.getModel();
        //
        // Add all existing seismic volumes
        //
        virtualVolumes = (StsSensorVirtualVolume[])model.getCastObjectList(StsSensorVirtualVolume.class);
        volumeCombo.initialize(this, "volume", "Volume:", virtualVolumes);
        wizard.rebuild();
    }

    void jbInit() throws Exception
    {
        progressPanel = StsProgressPanel.constructorWithCancelButton();

        addEndRow(volumeCombo);

        analyzeBox.gbc.fill = gbc.HORIZONTAL;
        msgBean.setText("Push button to start analysis");
        StsJPanel msgPanel = new StsJPanel();
        msgPanel.gbc.fill = gbc.HORIZONTAL;
        msgPanel.addEndRow(msgBean);

        StsJPanel btnPanel = new StsJPanel();
        unitsBean.initialize(wizard, "units", "Units:");
        unitsBean.setListItems(wizard.unitStrings, wizard.unitStrings[0]);
        runButton = new StsButton("Run Analysis", "Compute the volume within stimulated area.", this, "analyzeVolume", progressPanel);
        btnPanel.gbc.fill = gbc.NONE;
        btnPanel.addToRow(runButton);
        btnPanel.gbc.fill = gbc.HORIZONTAL;
        btnPanel.addEndRow(unitsBean);

        analyzeBox.addEndRow(msgPanel);
        analyzeBox.addEndRow(progressPanel);
        analyzeBox.addEndRow(btnPanel);
        add(analyzeBox);
    }

    public void analyzeVolume(StsProgressPanel panel)
    {
        progressPanel.resetProgressBar();
        progressPanel.setValue(0);
        wizard.analyzeVolume(panel);
    }

    public StsSeismicVolume getVolume()
    {
        return selectedVolume;
    }

    public void setVolume(StsSeismicVolume volume)
    {
        selectedVolume = volume;
    }

    public void setMessage(String msg)
    {
    	msgBean.setText(msg);
    }
}