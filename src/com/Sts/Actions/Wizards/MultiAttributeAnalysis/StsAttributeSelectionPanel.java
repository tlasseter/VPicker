package com.Sts.Actions.Wizards.MultiAttributeAnalysis;

import com.Sts.Actions.Wizards.*;
import com.Sts.DBTypes.*;
import com.Sts.MVC.*;
import com.Sts.Types.*;
import com.Sts.UI.Beans.*;
import com.Sts.UI.*;
import com.Sts.Utilities.*;

import javax.swing.*;
import java.awt.*;
import java.awt.event.*;

/**
 * <p>Title: </p>
 * <p>Description: </p>
 * <p>Copyright: Copyright (c) 2002</p>
 * <p>Company: </p>
 * @author unascribed
 * @version 1.0
 */

public class StsAttributeSelectionPanel extends JPanel implements ActionListener
{
    private StsWizard wizard;
    private StsWizardStep wizardStep;

    private StsModel model = null;

    private StsSeismicVolume[] seismicVolumes = null;
    private StsVirtualVolume[] virtualVolumes = null;

    private int nSeismicVolumes = 0;
    private int nVirtualVolumes = 0;

    private int selectedAzimuthIndex = 0;
    private int selectedColorIndex = 0;
    private int selectedLengthIndex = 0;

    StsCheckbox lengthThresholdChkBox = new StsCheckbox("Enable Threshold","Ignore values below the threshold value");
    StsFloatFieldBean lengthThresholdBean = new StsFloatFieldBean();
    private float lengthThreshold = StsParameters.nullValue;

    JComboBox azimuthCombo = new JComboBox();
    JLabel azimuthLabel = new JLabel();
    StsCheckbox azimuthChkBox = new StsCheckbox("Normalize","Scale values to between 0-360 degrees.");
    JComboBox colorCombo = new JComboBox();
    JLabel colorLabel = new JLabel();
    JComboBox lengthCombo = new JComboBox();
    JLabel lengthLabel = new JLabel();

    JPanel jPanel3 = new JPanel();
    GridBagLayout gridBagLayout1 = new GridBagLayout();
    GridBagLayout gridBagLayout2 = new GridBagLayout();

    public StsAttributeSelectionPanel(StsWizard wizard, StsWizardStep wizardStep)
    {
        this.wizard = wizard;
        this.wizardStep = wizardStep;
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
        seismicVolumes = ((StsMultiAttrAnalysisWizard)wizard).getSeismicVolumes();
        nSeismicVolumes = seismicVolumes.length;
        azimuthCombo.addItem("None");
        lengthCombo.addItem("None");
        colorCombo.addItem("None");
        for(int v = 0; v < nSeismicVolumes; v++)
        {
            azimuthCombo.addItem(seismicVolumes[v].getName());
            lengthCombo.addItem(seismicVolumes[v].getName());
            colorCombo.addItem(seismicVolumes[v].getName());
        }
        virtualVolumes = ((StsMultiAttrAnalysisWizard)wizard).getVirtualVolumes();
        nVirtualVolumes = virtualVolumes.length;
        for(int v = 0; v < nVirtualVolumes; v++)
        {
            azimuthCombo.addItem(virtualVolumes[v].getName());
            lengthCombo.addItem(virtualVolumes[v].getName());
            colorCombo.addItem(virtualVolumes[v].getName());
        }
        if((nSeismicVolumes + nVirtualVolumes) > 0)
        {
            setSelectedAzimuthVolume(0);
            setSelectedLengthVolume(0);
            setSelectedColorVolume(0);
        }
        lengthThresholdChkBox.setEnabled(false);

    }

    // Set the seismic combobox to the current index
    private void setSelectedLengthVolume(int volumeIndex)
    {
        lengthCombo.setSelectedIndex(volumeIndex);
    }
    // Set the seismic combobox to the current index
    private void setSelectedAzimuthVolume(int volumeIndex)
    {
        azimuthCombo.setSelectedIndex(volumeIndex);
    }
    // Set the seismic combobox to the current index
    private void setSelectedColorVolume(int volumeIndex)
    {
        colorCombo.setSelectedIndex(volumeIndex);
    }

    // Get the current seismic index
    public int getLengthVolumeIndex()
    {
        return selectedLengthIndex;
    }
    // Get the current seismic index
    public int getColorVolumeIndex()
    {
        return selectedColorIndex;
    }
    // Get the current seismic index
    public int getAzimuthVolumeIndex()
    {
        return selectedAzimuthIndex;
    }
    void jbInit() throws Exception
    {
        this.setLayout(gridBagLayout2);
        azimuthLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        azimuthLabel.setText("Azimuth PostStack3d:");
        colorLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        colorLabel.setText("Color PostStack3d:");
        lengthLabel.setFont(new java.awt.Font("Dialog", 3, 11));
        lengthLabel.setText("Length PostStack3d:");

        lengthThresholdBean = new StsFloatFieldBean(this, "lengthThreshold", true, null, true);

        azimuthCombo.setToolTipText("Select seismic volume to represent vector azimuth");
        lengthCombo.setToolTipText("Select seismic volume to represent vector length");
        colorCombo.setToolTipText("Select seismic volume to represent vector color");

        jPanel3.setBorder(BorderFactory.createEtchedBorder());
        jPanel3.setLayout(gridBagLayout1);

        jPanel3.add(azimuthLabel, new GridBagConstraints(0, 0, 1, 1, 0.0, 0.0
            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 2), 0, 0));
        jPanel3.add(azimuthCombo,     new GridBagConstraints(1, 0, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        jPanel3.add(azimuthChkBox,     new GridBagConstraints(2, 0, 2, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));

        jPanel3.add(lengthLabel, new GridBagConstraints(0, 1, 1, 1, 0.0, 0.0
            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 2), 0, 0));
        jPanel3.add(lengthCombo,     new GridBagConstraints(1, 1, 1, 1, 0.5, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));
        jPanel3.add(lengthThresholdChkBox,     new GridBagConstraints(2, 1, 1, 1, 0.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.NONE, new Insets(0, 0, 0, 0), 0, 0));
        jPanel3.add(lengthThresholdBean,     new GridBagConstraints(3, 1, 1, 1, 0.5, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        jPanel3.add(colorLabel, new GridBagConstraints(0, 2, 1, 1, 0.0, 0.0
            , GridBagConstraints.EAST, GridBagConstraints.NONE, new Insets(0, 0, 5, 2), 0, 0));
        jPanel3.add(colorCombo,     new GridBagConstraints(1, 2, 1, 1, 1.0, 0.0
            ,GridBagConstraints.WEST, GridBagConstraints.HORIZONTAL, new Insets(0, 0, 0, 0), 0, 0));

        this.add(jPanel3, new GridBagConstraints(0, 0, 1, 1, 1.0, 1.0
            , GridBagConstraints.CENTER, GridBagConstraints.BOTH, new Insets(4, 3, 5, 4), 3, -5));

        azimuthCombo.addActionListener(this);
        lengthCombo.addActionListener(this);
        colorCombo.addActionListener(this);
        azimuthChkBox.addActionListener(this);
        lengthThresholdChkBox.addActionListener(this);
    }

    public void actionPerformed(ActionEvent e)
    {
        Object source = e.getSource();

        if(source == lengthCombo)
        {
            selectedLengthIndex = lengthCombo.getSelectedIndex();
            StsSeismicVolume volume = (StsSeismicVolume)getSelectedLengthVolume();
            if(volume != null)
            {
                lengthThresholdBean.setEnabled(true);
                lengthThresholdChkBox.setEnabled(true);
                float step = (volume.getDataMin() - volume.getDataMax()) / 1000.0f;
                lengthThresholdBean.setValueAndRangeFixStep(volume.getDataMin(), volume.getDataMin(), volume.getDataMax(), step);
            }
            else
            {
                lengthThresholdBean.setEnabled(false);
                lengthThresholdChkBox.setEnabled(false);
            }
        }
        else if(source == colorCombo)
            selectedColorIndex = colorCombo.getSelectedIndex();
        else if(source == azimuthCombo)
        {
            selectedAzimuthIndex = azimuthCombo.getSelectedIndex();
            isAzimuthVolumeInRange();
        }
        else if(source == azimuthChkBox)
            isAzimuthVolumeInRange();
        else if(source == lengthThresholdChkBox)
        {
            if(lengthThresholdChkBox.isSelected())
                lengthThresholdChkBox.setEnabled(true);
            else
                lengthThresholdChkBox.setEnabled(false);
        }
    }

    public float getLengthThreshold() { return lengthThreshold; }
    public void setLengthThreshold(float threshold) { lengthThreshold = threshold; }
    public StsSeismicBoundingBox getSelectedAzimuthVolume()
    {
        int idx = selectedAzimuthIndex;
        idx = idx -1;
        if(idx == -1)
            return null;
        if(idx > seismicVolumes.length-1)
            return virtualVolumes[idx - seismicVolumes.length];
        else
            return seismicVolumes[idx];
    }

    public StsSeismicBoundingBox getSelectedColorVolume()
    {
        int idx = selectedColorIndex;
        idx = idx -1;
        if(idx == -1)
            return null;
        if(idx > seismicVolumes.length-1)
            return virtualVolumes[idx - seismicVolumes.length];
        else
            return seismicVolumes[idx];
    }

    public StsSeismicBoundingBox getSelectedLengthVolume()
    {
        int idx = selectedLengthIndex;
        idx = idx -1;
        if(idx == -1)
            return null;
        if(idx > seismicVolumes.length-1)
            return virtualVolumes[idx - seismicVolumes.length];
        else
            return seismicVolumes[idx];
    }

    public StsSeismicBoundingBox[] getSelectedVolumes()
    {
        StsSeismicBoundingBox[] volumes = new StsSeismicBoundingBox[3];  // Order dependent
        volumes[StsMultiAttributeVector.AZIMUTH] = getSelectedAzimuthVolume();
        volumes[StsMultiAttributeVector.LENGTH] = getSelectedLengthVolume();
        volumes[StsMultiAttributeVector.COLOR] = getSelectedColorVolume();
        return volumes;
    }

    private boolean isAzimuthVolumeInRange()
    {
        // Verify that the value range for the azimuth attribute is between 0-360 if user selects no normalization
        boolean normalize = getNormalizeAzimuth();
        StsSeismicBoundingBox azVolume = getSelectedAzimuthVolume();
        if(azVolume == null) return true;
        if(!normalize)
        {
            if((azVolume.getDataMax() > 360.0) && (azVolume.getDataMin() < 0.0))
            {
                if(StsYesNoDialog.questionValue(wizard.frame,"Must normalize azimuth volume since range is outside 0-360"))
                    azimuthChkBox.setSelected(true);
                else
                {
                    new StsMessage(wizard.frame, StsMessage.WARNING, "Then must select a different volume that is in range.",true);
                    selectedAzimuthIndex = 0;
                    azimuthCombo.setSelectedIndex(selectedAzimuthIndex);
                    return false;
                }
            }
        }
        return true;
    }

    public boolean getNormalizeAzimuth()
    {
        return azimuthChkBox.isSelected();
    }
}
